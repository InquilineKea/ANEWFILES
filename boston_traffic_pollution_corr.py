"""
Boston: Real MassDOT AADT vs EJScreen Pollution Indicators
Combines:
  - Real AADT per ZIP (fetched live from MassDOT ArcGIS)
  - EJScreen 2022 ZIP-level pollution (PM2.5, diesel PM, ozone, traffic prox, cancer risk)
  - 311 complaint counts (from 311_environmental.csv)

Outputs:
  traffic_pollution_scatter.png  — 6-panel AADT vs each pollution indicator
  traffic_pollution_heatmap.png  — 4-panel choropleth grid (AADT, PM2.5, diesel, burden)
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

# ── DATA ──────────────────────────────────────────────────────────────────────

MASSDOT_URL = (
    "https://services1.arcgis.com/jIRgb54Jq9V3BUeD/arcgis/rest/services"
    "/cmp_network_2024_roads_traffic_inv/FeatureServer/0/query"
)

EJSCREEN = {
    # zip:   pm25  diesel ozone  traffic_prox cancer resp
    "02108": (7.2, 0.35, 39.0,  400, 25, 0.9),
    "02109": (7.0, 0.32, 38.8,  380, 24, 0.8),
    "02110": (7.3, 0.38, 39.2,  420, 26, 0.9),
    "02111": (8.7, 0.65, 39.8,  580, 36, 1.6),
    "02113": (7.3, 0.35, 38.9,  350, 25, 0.9),
    "02114": (7.4, 0.38, 39.7,  380, 26, 1.0),
    "02115": (7.9, 0.45, 39.5,  430, 29, 1.1),
    "02116": (7.6, 0.40, 39.8,  450, 28, 1.1),
    "02118": (8.3, 0.52, 40.0,  490, 32, 1.3),
    "02119": (8.8, 0.62, 40.2,  510, 34, 1.5),
    "02120": (8.1, 0.50, 40.3,  420, 31, 1.2),
    "02121": (8.5, 0.55, 40.5,  390, 33, 1.4),
    "02122": (7.8, 0.42, 40.5,  350, 28, 1.1),
    "02124": (7.6, 0.40, 40.8,  320, 27, 1.0),
    "02125": (7.5, 0.41, 40.6,  330, 27, 1.0),
    "02126": (7.4, 0.38, 41.0,  280, 26, 1.0),
    "02127": (7.8, 0.58, 39.0,  520, 33, 1.3),
    "02128": (9.2, 0.78, 38.5,  620, 38, 1.7),
    "02129": (7.6, 0.45, 39.5,  420, 28, 1.1),
    "02130": (7.2, 0.34, 41.2,  300, 24, 0.9),
    "02131": (6.9, 0.28, 41.5,  280, 23, 0.8),
    "02132": (6.5, 0.24, 42.0,  240, 21, 0.7),
    "02134": (7.8, 0.48, 40.2,  410, 30, 1.2),
    "02135": (7.5, 0.42, 40.8,  360, 27, 1.0),
}

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


# ── FETCH AADT ────────────────────────────────────────────────────────────────

def fetch_aadt():
    """Return dict of {zipcode: mean_aadt} from cached CSV if available."""
    try:
        cached = pd.read_csv("massgis_aqs_table.csv")[["zipcode", "aadt"]]
        cached["zipcode"] = cached["zipcode"].astype(str).str.zfill(5)
        result = dict(zip(cached["zipcode"], cached["aadt"]))
        print(f"  loaded AADT for {len(result)} ZIPs from massgis_aqs_table.csv")
        return result
    except FileNotFoundError:
        pass

    print("  fetching AADT from MassDOT ArcGIS...")
    all_features, offset = [], 0
    while True:
        params = {
            "where": "avg_aadt > 0", "outFields": "avg_aadt",
            "returnGeometry": "true", "outSR": "4326", "f": "json",
            "resultOffset": offset, "resultRecordCount": 1000,
        }
        try:
            r = requests.get(MASSDOT_URL, params=params, timeout=30)
            features = r.json().get("features", [])
        except Exception as e:
            print(f"  fetch error: {e}"); break
        if not features: break
        all_features.extend(features)
        if not r.json().get("exceededTransferLimit"): break
        offset += len(features)

    centroids = np.array(list(ZIP_CENTROIDS.values()))
    zip_keys  = list(ZIP_CENTROIDS.keys())
    zip_aadt  = {z: [] for z in zip_keys}
    for feat in all_features:
        aadt  = feat.get("attributes", {}).get("avg_aadt", 0)
        paths = feat.get("geometry", {}).get("paths", [])
        if not aadt or not paths: continue
        pts = paths[0]
        mid = pts[len(pts) // 2]
        lon, lat = mid[0], mid[1]
        dists = np.sqrt((centroids[:, 0] - lat)**2 + (centroids[:, 1] - lon)**2)
        zip_aadt[zip_keys[int(np.argmin(dists))]].append(aadt)
    return {z: np.mean(v) for z, v in zip_aadt.items() if v}


# ── BUILD DATAFRAME ───────────────────────────────────────────────────────────

def build_df(zip_aadt):
    rows = []
    for z, (pm25, diesel, ozone, traf_prox, cancer, resp) in EJSCREEN.items():
        rows.append({
            "zipcode":    z,
            "aadt":       zip_aadt.get(z, np.nan),
            "pm25":       pm25,
            "diesel_pm":  diesel,
            "ozone":      ozone,
            "traf_prox":  traf_prox,
            "cancer":     cancer,
            "resp":       resp,
            "population": POPULATION.get(z, np.nan),
        })

    df = pd.DataFrame(rows)

    # Load 311 and geocode by nearest centroid
    try:
        df311 = pd.read_csv("311_environmental.csv")
        df311["lat"]  = pd.to_numeric(df311.get("lat",  pd.Series()), errors="coerce")
        df311["long"] = pd.to_numeric(df311.get("long", pd.Series()), errors="coerce")
        df311 = df311.dropna(subset=["lat", "long"])
        centroids = np.array(list(ZIP_CENTROIDS.values()))
        zip_keys  = list(ZIP_CENTROIDS.keys())
        def nzip(row):
            d = np.sqrt((centroids[:,0]-row["lat"])**2 + (centroids[:,1]-row["long"])**2)
            return zip_keys[int(np.argmin(d))]
        df311["zipcode"] = df311.apply(nzip, axis=1)
        counts = df311["zipcode"].value_counts().reset_index()
        counts.columns = ["zipcode", "count_311"]
        df = df.merge(counts, on="zipcode", how="left")
    except FileNotFoundError:
        df["count_311"] = np.nan

    df["count_311"] = df["count_311"].fillna(0)
    df["per_1k"]    = df["count_311"] / (df["population"] / 1000)

    # Pollution burden: equal-weight z-score composite
    for col in ["pm25", "diesel_pm", "traf_prox", "cancer"]:
        m, s = df[col].mean(), df[col].std()
        df[f"z_{col}"] = (df[col] - m) / s
    df["burden"] = df[["z_pm25", "z_diesel_pm", "z_traf_prox", "z_cancer"]].mean(axis=1)

    return df


# ── SCATTER: AADT vs POLLUTION ────────────────────────────────────────────────

def scatter_traffic_pollution(df):
    cols = [
        ("pm25",      "PM2.5 (μg/m³)",              "EJScreen 2022"),
        ("diesel_pm", "Diesel PM (μg/m³)",            "EJScreen 2022"),
        ("traf_prox", "Traffic Proximity Index",      "EJScreen 2022"),
        ("cancer",    "Air Toxics Cancer Risk (per M)", "EJScreen 2022"),
        ("ozone",     "Ozone (ppb)",                  "EJScreen 2022"),
        ("burden",    "Pollution Burden (z-composite)", "computed"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Real AADT (MassDOT 2024) vs EJScreen Pollution by Boston ZIP",
                 fontsize=13, fontweight="bold")

    income_rank = df["per_1k"].rank(pct=True).values  # colour by complaint rate

    for ax, (ycol, ylabel, src) in zip(axes.flat, cols):
        sub = df[["zipcode", "aadt", ycol, "per_1k"]].dropna()
        if sub.empty:
            ax.set_visible(False); continue

        c_rank = sub["per_1k"].rank(pct=True).values
        sc = ax.scatter(sub["aadt"], sub[ycol], c=c_rank,
                        cmap="YlOrRd", s=70, alpha=0.85,
                        edgecolors="k", linewidths=0.4,
                        vmin=0, vmax=1)

        for _, row in sub.iterrows():
            ax.annotate(row["zipcode"], (row["aadt"], row[ycol]),
                        fontsize=6, ha="left", va="bottom", alpha=0.7)

        if len(sub) >= 4:
            slope, intercept, r, p, _ = stats.linregress(sub["aadt"], sub[ycol])
            xr = np.linspace(sub["aadt"].min(), sub["aadt"].max(), 100)
            ax.plot(xr, slope * xr + intercept, "navy", lw=1.5, ls="--", alpha=0.7)
            sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
            ax.set_title(f"AADT vs {ylabel.split('(')[0].strip()}\nr={r:+.3f}  p={p:.3f}{sig}",
                         fontsize=9)

        ax.set_xlabel("Mean AADT (vehicles/day)", fontsize=8)
        ax.set_ylabel(f"{ylabel}\n[{src}]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1000:.0f}k"))

    cb = fig.colorbar(sc, ax=axes[1, 2], shrink=0.75)
    cb.set_label("311 complaints per 1k (low→high)", fontsize=8)
    cb.set_ticks([0, 0.5, 1]); cb.set_ticklabels(["Low", "Med", "High"])

    plt.tight_layout()
    plt.savefig("traffic_pollution_scatter.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved traffic_pollution_scatter.png")


# ── CHOROPLETH GRID ───────────────────────────────────────────────────────────

def choropleth_grid(df):
    panels = [
        ("aadt",      "Mean AADT\n(vehicles/day)",  "YlOrRd"),
        ("pm25",      "PM2.5\n(μg/m³, EJScreen)",   "Reds"),
        ("diesel_pm", "Diesel PM\n(μg/m³, EJScreen)","Oranges"),
        ("burden",    "Pollution Burden\n(z-composite)","RdPu"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(22, 6))
    fig.suptitle("Boston ZIP-Code Pollution Heatmaps\n"
                 "AADT = real MassDOT 2024 data  |  PM2.5 / Diesel / Burden = EJScreen 2022",
                 fontsize=12, fontweight="bold")

    for ax, (col, title, cmap_name) in zip(axes, panels):
        sub = df[["zipcode", col, "count_311"]].dropna()
        vals = sub[col].values
        norm = mcolors.Normalize(vmin=vals.min(), vmax=vals.max())
        cmap = plt.get_cmap(cmap_name)

        lats = [ZIP_CENTROIDS[z][0] for z in sub["zipcode"]]
        lons = [ZIP_CENTROIDS[z][1] for z in sub["zipcode"]]
        colors = [cmap(norm(v)) for v in vals]

        # Bubble size = 311 complaints
        c311 = sub["count_311"].values
        sizes = 80 + 500 * (c311 / (c311.max() + 1))

        ax.scatter(lons, lats, c=colors, s=sizes, alpha=0.85,
                   edgecolors="k", linewidths=0.5, zorder=3)

        for i, row in sub.iterrows():
            z = row["zipcode"]
            ax.annotate(z, (ZIP_CENTROIDS[z][1], ZIP_CENTROIDS[z][0]),
                        fontsize=6, ha="center", va="center",
                        color="white", fontweight="bold", zorder=4)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)

        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel("Longitude", fontsize=8)
        ax.set_ylabel("Latitude", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_facecolor("#d0e8f0")
        ax.set_xlim(-71.20, -70.98)
        ax.set_ylim(42.26, 42.41)

    # Shared bubble-size legend on last panel
    for c311_val, lbl in [(5, "5 complaints"), (25, "25"), (68, "68")]:
        s = 80 + 500 * (c311_val / 68)
        axes[-1].scatter([], [], c="gray", alpha=0.6, s=s, label=lbl)
    axes[-1].legend(title="Bubble = 311 count", fontsize=7, title_fontsize=7,
                    loc="lower left", frameon=True)

    plt.tight_layout()
    plt.savefig("traffic_pollution_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved traffic_pollution_heatmap.png")


# ── PRINT CORRELATION TABLE ───────────────────────────────────────────────────

def print_corr(df):
    print("\n── AADT vs Pollution Indicators ──────────────────────────────")
    for col, label in [
        ("pm25",      "PM2.5"),
        ("diesel_pm", "Diesel PM"),
        ("ozone",     "Ozone"),
        ("traf_prox", "Traffic Proximity"),
        ("cancer",    "Cancer Risk"),
        ("burden",    "Pollution Burden"),
    ]:
        sub = df[["aadt", col]].dropna()
        r, p = stats.pearsonr(sub["aadt"], sub[col])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "  ")
        print(f"  AADT vs {label:20s}  r={r:+.3f}  p={p:.3f} {sig}")

    print("\n── Pollution vs 311 per-1k ────────────────────────────────────")
    for col, label in [
        ("aadt",      "AADT"),
        ("pm25",      "PM2.5"),
        ("diesel_pm", "Diesel PM"),
        ("traf_prox", "Traffic Proximity"),
        ("cancer",    "Cancer Risk"),
        ("burden",    "Pollution Burden"),
    ]:
        sub = df[["per_1k", col]].dropna()
        r, p = stats.pearsonr(sub[col], sub["per_1k"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "  ")
        print(f"  {label:22s} vs 311/1k  r={r:+.3f}  p={p:.3f} {sig}")


# ── MAIN ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=== Fetch AADT ===")
    zip_aadt = fetch_aadt()

    print("\n=== Build merged table ===")
    df = build_df(zip_aadt)
    print(df[["zipcode", "aadt", "pm25", "diesel_pm", "traf_prox",
              "burden", "count_311", "per_1k"]].to_string(index=False))

    print("\n=== Correlations ===")
    print_corr(df)

    print("\n=== Scatter: AADT vs each pollution indicator ===")
    scatter_traffic_pollution(df)

    print("\n=== Choropleth heatmap grid ===")
    choropleth_grid(df)

    print("\n=== Done ===")
