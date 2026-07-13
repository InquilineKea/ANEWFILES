"""
Boston 311 × MassGIS Traffic + EPA AQS Air Quality
Sources:
  - MassDOT AADT: ArcGIS cmp_network_2024_roads_traffic_inv
  - EPA AQS 2023 annual monitors: Suffolk County stations
  - 311 complaints: 311_environmental.csv
  - Census ACS 2022 income/poverty: hardcoded

Method:
  1. Fetch 1000 road segments from MassDOT ArcGIS → assign to ZIP by nearest centroid
  2. Aggregate mean AADT per ZIP (weighted by segment count)
  3. Interpolate 2 EPA stations (Kenmore Sq, Harrison Ave) to all 24 ZIPs via IDW
  4. Cross-correlate traffic + air quality with 311 complaints by ZIP

Outputs:
  massgis_aqs_scatter.png   — 6-panel scatter (AADT, PM2.5, NO2, benzene, ozone, burden)
  massgis_aqs_map.png       — dual choropleth map
  massgis_aqs_table.csv     — ZIP-level merged table
"""

import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from scipy import stats

# ── CONSTANTS ─────────────────────────────────────────────────────────────────

MASSDOT_URL = (
    "https://services1.arcgis.com/jIRgb54Jq9V3BUeD/arcgis/rest/services"
    "/cmp_network_2024_roads_traffic_inv/FeatureServer/0/query"
)

# EPA AQS 2023 annual averages for Suffolk County monitoring stations
# Station 1: Kenmore Square (25-025-0042)
# Station 2: Harrison Ave / Roxbury (25-025-0044)
EPA_STATIONS = [
    {
        "name": "Kenmore Sq",
        "lat": 42.348873,
        "lon": -71.097163,
        "pm25":    8.416,   # μg/m³ annual mean
        "no2":    21.24,    # ppb annual mean
        "so2":     0.40,    # ppb annual mean
        "ozone":   0.036,   # ppm annual mean (4th max 8-hr)
        "benzene": 0.85,    # ppbC (estimated lower — primarily traffic monitor)
        "bc":      0.35,    # μg/m³ black carbon
    },
    {
        "name": "Harrison Ave",
        "lat": 42.329500,
        "lon": -71.082600,
        "pm25":    8.10,    # μg/m³ (speciated PM2.5 FRM)
        "no2":    18.50,    # ppb
        "so2":     0.38,    # ppb
        "ozone":   0.039,   # ppm
        "benzene": 1.13,    # ppbC (NATTS toxic air monitor)
        "bc":      0.426,   # μg/m³ black carbon (aethalometer)
    },
]

ZIP_CENTROIDS = {
    "02108": (42.3601, -71.0549), "02109": (42.3601, -71.0518),
    "02110": (42.3572, -71.0504), "02111": (42.3630, -71.0596),
    "02113": (42.3658, -71.0560), "02114": (42.3614, -71.0706),
    "02115": (42.3447, -71.0875), "02116": (42.3508, -71.0757),
    "02118": (42.3361, -71.0714), "02119": (42.3199, -71.0875),
    "02120": (42.3308, -71.0931), "02121": (42.3077, -71.0864),
    "02122": (42.2993, -71.0691), "02124": (42.2935, -71.0797),
    "02125": (42.3077, -71.0704), "02126": (42.2801, -71.1003),
    "02127": (42.3341, -71.0376), "02128": (42.3700, -71.0155),
    "02129": (42.3839, -71.0576), "02130": (42.3096, -71.1148),
    "02131": (42.2906, -71.1249), "02132": (42.2786, -71.1519),
    "02134": (42.3508, -71.1102), "02135": (42.3536, -71.1502),
}

POPULATION = {
    "02108": 5575,  "02109": 3851,  "02110": 4260,  "02111": 14120,
    "02113": 9870,  "02114": 14710, "02115": 38350, "02116": 26180,
    "02118": 26420, "02119": 29840, "02120": 14340, "02121": 28100,
    "02122": 24980, "02124": 34740, "02125": 32060, "02126": 29280,
    "02127": 34150, "02128": 45770, "02129": 18410, "02130": 34750,
    "02131": 36710, "02132": 28880, "02134": 12650, "02135": 32580,
}

INCOME = {
    "02108": 170379, "02109": 162179, "02110": 126157, "02111":  72321,
    "02113": 106625, "02114": 121418, "02115":  51454, "02116": 131648,
    "02118":  74155, "02119":  37158, "02120":  51698, "02121":  42411,
    "02122":  81025, "02124":  72189, "02125":  64983, "02126":  66689,
    "02127": 145357, "02128":  84023, "02129": 146815, "02130": 114933,
    "02131": 101566, "02132": 132457, "02134":  80027, "02135":  87603,
}

POVERTY = {
    "02108": 14.7, "02109":  7.1, "02110": 14.8, "02111": 27.8,
    "02113": 11.0, "02114": 11.3, "02115": 33.9, "02116": 15.1,
    "02118": 26.2, "02119": 30.9, "02120": 38.8, "02121": 24.5,
    "02122": 14.8, "02124": 17.6, "02125": 23.5, "02126": 18.7,
    "02127": 10.5, "02128": 14.2, "02129":  8.4, "02130": 12.3,
    "02131": 10.0, "02132":  4.1, "02134": 21.4, "02135": 19.2,
}


# ── STEP 1: FETCH MASSDOT AADT ────────────────────────────────────────────────

def fetch_massdot_aadt(max_records=2000):
    """Fetch road segments with AADT from MassDOT ArcGIS FeatureServer."""
    all_features = []
    offset = 0
    batch = 1000

    while len(all_features) < max_records:
        params = {
            "where": "avg_aadt > 0",
            "outFields": "route_id,avg_aadt",
            "returnGeometry": "true",
            "geometryType": "esriGeometryPolyline",
            "outSR": "4326",
            "f": "json",
            "resultOffset": offset,
            "resultRecordCount": min(batch, max_records - len(all_features)),
        }
        try:
            r = requests.get(MASSDOT_URL, params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f"  MassDOT fetch error at offset {offset}: {e}")
            break

        features = data.get("features", [])
        if not features:
            break
        all_features.extend(features)
        print(f"  fetched {len(all_features)} road segments...")

        if not data.get("exceededTransferLimit", False):
            break
        offset += len(features)

    print(f"  total road segments: {len(all_features)}")
    return all_features


def segments_to_zip(features):
    """Assign each road segment to nearest ZIP centroid; aggregate AADT per ZIP."""
    centroids = np.array(list(ZIP_CENTROIDS.values()))  # (N, 2) lat/lon
    zip_keys  = list(ZIP_CENTROIDS.keys())

    zip_aadt = {z: [] for z in zip_keys}

    for feat in features:
        aadt = feat.get("attributes", {}).get("avg_aadt", 0)
        if not aadt:
            continue
        geom = feat.get("geometry", {})
        paths = geom.get("paths", [])
        if not paths:
            continue

        # Use midpoint of first path
        pts = paths[0]
        if not pts:
            continue
        mid = pts[len(pts) // 2]
        lon, lat = mid[0], mid[1]

        # Nearest centroid (Euclidean OK at this scale)
        dists = np.sqrt((centroids[:, 0] - lat) ** 2 + (centroids[:, 1] - lon) ** 2)
        nearest = zip_keys[int(np.argmin(dists))]
        zip_aadt[nearest].append(aadt)

    result = {}
    for z, vals in zip_aadt.items():
        result[z] = np.mean(vals) if vals else np.nan
        if vals:
            print(f"  ZIP {z}: {len(vals)} segments, mean AADT {result[z]:,.0f}")

    return result


# ── STEP 2: IDW INTERPOLATION FROM EPA STATIONS ───────────────────────────────

def idw_interpolate(station_vals, zip_lat, zip_lon, power=2):
    """Inverse distance weighting from EPA stations to a single ZIP centroid."""
    weights = []
    for st in EPA_STATIONS:
        d = np.sqrt((st["lat"] - zip_lat) ** 2 + (st["lon"] - zip_lon) ** 2)
        d = max(d, 1e-6)
        weights.append(1.0 / d ** power)

    total_w = sum(weights)
    result = {}
    for key in station_vals[0].keys():
        result[key] = sum(w * st[key] for w, st in zip(weights, station_vals)) / total_w
    return result


def interpolate_all_zips():
    """IDW-interpolate all EPA pollutant fields to every ZIP centroid."""
    station_vals = [
        {k: v for k, v in st.items() if k not in ("name", "lat", "lon")}
        for st in EPA_STATIONS
    ]
    rows = {}
    for z, (lat, lon) in ZIP_CENTROIDS.items():
        rows[z] = idw_interpolate(station_vals, lat, lon, power=2)
    return rows


# ── STEP 3: LOAD 311 COMPLAINTS ───────────────────────────────────────────────

def load_311():
    try:
        df = pd.read_csv("311_environmental.csv")
        print(f"  loaded {len(df)} 311 rows from CSV")
    except FileNotFoundError:
        print("  311_environmental.csv not found — returning empty df")
        return pd.DataFrame()

    # Try explicit zipcode column first
    zip_col = None
    for c in df.columns:
        if "zip" in c.lower():
            zip_col = c
            break

    if zip_col:
        df["zipcode"] = df[zip_col].astype(str).str.strip().str.zfill(5)
    elif "lat" in df.columns and "long" in df.columns:
        # Assign to nearest ZIP centroid by lat/long
        df["lat"]  = pd.to_numeric(df["lat"],  errors="coerce")
        df["long"] = pd.to_numeric(df["long"], errors="coerce")
        df = df.dropna(subset=["lat", "long"])
        centroids = np.array(list(ZIP_CENTROIDS.values()))
        zip_keys  = list(ZIP_CENTROIDS.keys())

        def nearest_zip(row):
            dists = np.sqrt((centroids[:, 0] - row["lat"]) ** 2
                            + (centroids[:, 1] - row["long"]) ** 2)
            return zip_keys[int(np.argmin(dists))]

        df["zipcode"] = df.apply(nearest_zip, axis=1)
        print(f"  assigned {len(df)} rows to ZIPs via lat/long nearest centroid")
    else:
        return pd.DataFrame()

    counts = df["zipcode"].value_counts().reset_index()
    counts.columns = ["zipcode", "count_311"]
    return counts


# ── STEP 4: BUILD MERGED TABLE ────────────────────────────────────────────────

def build_merged(zip_aadt, zip_aq, counts_311):
    rows = []
    for z in ZIP_CENTROIDS:
        row = {
            "zipcode":    z,
            "aadt":       zip_aadt.get(z, np.nan),
            "income":     INCOME.get(z, np.nan),
            "poverty":    POVERTY.get(z, np.nan),
            "population": POPULATION.get(z, np.nan),
        }
        aq = zip_aq.get(z, {})
        row.update({
            "pm25":    aq.get("pm25",    np.nan),
            "no2":     aq.get("no2",     np.nan),
            "ozone":   aq.get("ozone",   np.nan),
            "benzene": aq.get("benzene", np.nan),
            "bc":      aq.get("bc",      np.nan),
        })
        rows.append(row)

    df = pd.DataFrame(rows)
    if not counts_311.empty:
        df = df.merge(counts_311, on="zipcode", how="left")
    else:
        df["count_311"] = np.nan

    df["count_311"] = df["count_311"].fillna(0)
    df["per_1k"]    = df["count_311"] / (df["population"] / 1000)

    # Pollution burden: z-score composite (AADT + PM2.5 + NO2 + benzene)
    for col in ["aadt", "pm25", "no2", "benzene"]:
        m, s = df[col].mean(), df[col].std()
        df[f"z_{col}"] = (df[col] - m) / s if s > 0 else 0.0
    df["burden"] = df[["z_aadt", "z_pm25", "z_no2", "z_benzene"]].mean(axis=1)

    return df


# ── STEP 5: CORRELATION ANALYSIS ─────────────────────────────────────────────

def analyze(df):
    print("\n── Correlation with raw 311 count ──────────────────────────────")
    metrics = ["aadt", "pm25", "no2", "ozone", "benzene", "bc", "burden",
               "income", "poverty"]
    for m in metrics:
        sub = df[["count_311", m]].dropna()
        if len(sub) < 4:
            continue
        r, p = stats.pearsonr(sub[m], sub["count_311"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        print(f"  {m:12s}  r={r:+.3f}  p={p:.3f} {sig}")

    print("\n── Correlation with per-1k 311 rate ────────────────────────────")
    for m in metrics:
        sub = df[["per_1k", m]].dropna()
        if len(sub) < 4:
            continue
        r, p = stats.pearsonr(sub[m], sub["per_1k"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        print(f"  {m:12s}  r={r:+.3f}  p={p:.3f} {sig}")


# ── STEP 6: SCATTER PLOTS ─────────────────────────────────────────────────────

def build_scatter(df):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Boston 311 Complaints × Real Traffic & Air Quality\n"
                 "(MassDOT AADT 2024 + EPA AQS 2023 IDW-interpolated)",
                 fontsize=13, fontweight="bold")

    pairs = [
        ("aadt",    "Mean AADT (vehicles/day)"),
        ("pm25",    "PM2.5 (μg/m³) — IDW from 2 stations"),
        ("no2",     "NO₂ (ppb) — IDW"),
        ("benzene", "Benzene (ppbC) — IDW"),
        ("ozone",   "Ozone (ppm) — IDW"),
        ("burden",  "Pollution Burden (z-score composite)"),
    ]

    colors = cm.RdYlGn_r(np.linspace(0.15, 0.85, len(df)))
    income_sorted = df["income"].rank(pct=True).values

    for ax, (xcol, xlabel) in zip(axes.flat, pairs):
        sub = df[["zipcode", xcol, "count_311", "per_1k", "income"]].dropna()
        if sub.empty:
            ax.set_visible(False)
            continue

        sc_colors = cm.RdYlGn_r(sub["income"].rank(pct=True).values)
        ax.scatter(sub[xcol], sub["per_1k"], c=sc_colors, s=60, alpha=0.8, edgecolors="k", lw=0.4)

        for _, row in sub.iterrows():
            ax.annotate(row["zipcode"], (row[xcol], row["per_1k"]),
                        fontsize=6, ha="left", va="bottom", alpha=0.7)

        # Regression line
        if len(sub) >= 4:
            m_slope, b_int, r, p, _ = stats.linregress(sub[xcol], sub["per_1k"])
            xr = np.linspace(sub[xcol].min(), sub[xcol].max(), 100)
            ax.plot(xr, m_slope * xr + b_int, "r--", lw=1.2, alpha=0.7)
            sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
            ax.set_title(f"vs {xlabel.split('(')[0].strip()}\nr={r:+.3f} p={p:.3f}{sig}", fontsize=9)
        else:
            ax.set_title(xlabel, fontsize=9)

        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("311 complaints per 1,000 residents", fontsize=8)
        ax.tick_params(labelsize=7)

    sm = plt.cm.ScalarMappable(cmap="RdYlGn_r", norm=mcolors.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.flat[-1], shrink=0.7)
    cbar.set_label("Median income (low→high)", fontsize=8)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["Low", "Med", "High"])

    plt.tight_layout()
    plt.savefig("massgis_aqs_scatter.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved massgis_aqs_scatter.png")


# ── STEP 7: DUAL CHOROPLETH MAP ───────────────────────────────────────────────

def build_map(df):
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle("Boston ZIPs: Traffic AADT (left) vs 311 per 1k (right)\n"
                 "Point size = PM2.5 (IDW), colour = AADT / complaint rate",
                 fontsize=12, fontweight="bold")

    lats = [c[0] for c in ZIP_CENTROIDS.values()]
    lons = [c[1] for c in ZIP_CENTROIDS.values()]

    for ax, val_col, title, cmap_name in [
        (axes[0], "aadt",    "Mean AADT (vehicles/day)", "YlOrRd"),
        (axes[1], "per_1k",  "311 complaints / 1k residents", "Blues"),
    ]:
        sub = df[["zipcode", val_col, "pm25"]].dropna()
        val_arr  = sub[val_col].values
        pm25_arr = sub["pm25"].values

        norm  = mcolors.Normalize(vmin=val_arr.min(), vmax=val_arr.max())
        cmap  = plt.get_cmap(cmap_name)
        colors = [cmap(norm(v)) for v in val_arr]

        # Bubble size proportional to PM2.5
        pm25_norm = (pm25_arr - pm25_arr.min()) / (np.ptp(pm25_arr) + 1e-6)
        sizes = 60 + 300 * pm25_norm

        zip_lat = [ZIP_CENTROIDS[z][0] for z in sub["zipcode"]]
        zip_lon = [ZIP_CENTROIDS[z][1] for z in sub["zipcode"]]

        ax.scatter(zip_lon, zip_lat, c=colors, s=sizes, alpha=0.85, edgecolors="k", lw=0.5)

        for i, row in sub.iterrows():
            ax.annotate(row["zipcode"],
                        (ZIP_CENTROIDS[row["zipcode"]][1], ZIP_CENTROIDS[row["zipcode"]][0]),
                        fontsize=6.5, ha="center", va="center", color="white", fontweight="bold")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, shrink=0.6, label=title)

        # Station markers
        for st in EPA_STATIONS:
            ax.plot(st["lon"], st["lat"], "D", color="cyan", markersize=8,
                    markeredgecolor="navy", lw=1.5, zorder=5)
            ax.annotate(st["name"], (st["lon"], st["lat"]),
                        fontsize=7, color="navy", ha="left",
                        xytext=(4, 4), textcoords="offset points")

        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title, fontsize=10)
        ax.set_facecolor("#e8f4f8")

    # Legend for bubble size
    for ax in axes:
        for sz, label in [(60, "Low PM2.5"), (200, "Mid"), (360, "High PM2.5")]:
            ax.scatter([], [], c="gray", alpha=0.6, s=sz, label=label)
        ax.legend(scatterpoints=1, frameon=True, labelspacing=0.8,
                  title="Bubble = PM2.5", fontsize=7, title_fontsize=7,
                  loc="lower left")

    plt.tight_layout()
    plt.savefig("massgis_aqs_map.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved massgis_aqs_map.png")


# ── MAIN ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=== Step 1: Fetch MassDOT AADT road segments ===")
    features = fetch_massdot_aadt(max_records=2000)

    print("\n=== Step 2: Assign segments → ZIP, aggregate AADT ===")
    zip_aadt = segments_to_zip(features)

    print("\n=== Step 3: IDW-interpolate EPA AQS stations → all ZIPs ===")
    zip_aq = interpolate_all_zips()
    sample = list(zip_aq.items())[:3]
    for z, aq in sample:
        print(f"  {z}: PM2.5={aq['pm25']:.3f} NO2={aq['no2']:.2f} "
              f"benz={aq['benzene']:.3f} ozone={aq['ozone']:.4f}")

    print("\n=== Step 4: Load 311 complaints ===")
    counts_311 = load_311()

    print("\n=== Step 5: Build merged table ===")
    df = build_merged(zip_aadt, zip_aq, counts_311)
    df.to_csv("massgis_aqs_table.csv", index=False)
    print(f"  saved massgis_aqs_table.csv ({len(df)} ZIPs)")
    print(df[["zipcode", "aadt", "pm25", "no2", "benzene", "burden",
               "count_311", "per_1k"]].to_string(index=False))

    print("\n=== Step 6: Correlation analysis ===")
    analyze(df)

    print("\n=== Step 7: Build scatter plots ===")
    build_scatter(df)

    print("\n=== Step 8: Build map ===")
    build_map(df)

    print("\n=== Done ===")
