"""
Boston 311 + Pollution Data Explorer
Sources:
  - 311 complaints:         https://311.boston.gov/open311/v2/requests.json
  - BERDO building emissions: resource 87521565-7f15-4b8d-a225-ac4df9e3f309
  - Community GHG emissions:  resource bd8dd4bb-867e-4ca2-b6c7-6c3bd9e6c290
  - Gov't operations GHG:     resource d2b83ed7-ce3c-448d-b20a-971ed911d9e7
  - Food establishment inspections: resource 4582bec6-2b4f-4f9e-bc55-cbaa73117f4c
  - Big Belly smart trash:    resource 15e7fa44-b9a8-42da-82e1-304e43460095

Outputs:
  311_environmental.csv       — filtered 311 env complaints
  berdo_emissions.csv         — BERDO building emissions
  boston_311_map.html         — interactive Folium map
  cross_correlation.png       — time-series + spatial correlation charts
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

CKAN_BASE  = "https://data.boston.gov/api/3/action/datastore_search"
OPEN311    = "https://311.boston.gov/open311/v2/requests.json"

ENVIRONMENTAL_311_TYPES = {
    "Rodent Activity", "Needle Pickup", "Needle Cleanup",
    "Litter", "Overflowing Trash Can", "Residential Trash out Illegally",
    "Illegal Dumping", "Graffiti Removal", "Encampment",
    "Dead Animal Collection", "Dead Animal Pickup",
}

# ── DATA FETCHING ─────────────────────────────────────────────────────────────

def fetch_311(days_back=730, max_records=2000):
    # Try pulling in 90-day chunks to avoid API timeouts on large ranges
    all_rows = []
    chunk = 90
    end   = datetime.now()
    start = end - timedelta(days=days_back)
    cursor = start
    while cursor < end and len(all_rows) < max_records:
        chunk_end   = min(cursor + timedelta(days=chunk), end)
        params = {
            "start_date": cursor.strftime("%Y-%m-%dT00:00:00Z"),
            "end_date":   chunk_end.strftime("%Y-%m-%dT00:00:00Z"),
            "page_size":  100,
            "status":     "open,closed",
            "page":       1,
        }
        while len(all_rows) < max_records:
            try:
                r = requests.get(OPEN311, params=params, timeout=15)
                r.raise_for_status()
                text = r.text.strip()
                if not text:
                    break
                batch = r.json()
            except Exception:
                break
            if not batch:
                break
            all_rows.extend(batch)
            if len(batch) < 100:
                break
            params["page"] += 1
        cursor = chunk_end

    # Fall back to saved CSV if API returned nothing
    if not all_rows:
        csv_path = "/home/user/ANEWFILES/311_environmental.csv"
        try:
            df = pd.read_csv(csv_path)
            print(f"  (loaded {len(df)} rows from {csv_path})")
        except FileNotFoundError:
            return pd.DataFrame()
    else:
        df = pd.DataFrame(all_rows)

    if not df.empty:
        for col in ["requested_datetime", "updated_datetime"]:
            if col in df.columns:
                df[col] = pd.to_datetime(df[col], errors="coerce", utc=True)
        for col in ["lat", "long"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def fetch_ckan(resource_id, limit=5000):
    rows, offset = [], 0
    while True:
        r = requests.get(CKAN_BASE, params={
            "resource_id": resource_id,
            "limit": min(limit - len(rows), 1000),
            "offset": offset,
        }, timeout=20)
        r.raise_for_status()
        result = r.json()["result"]
        batch  = result["records"]
        rows.extend(batch)
        if not batch or len(rows) >= result["total"] or len(rows) >= limit:
            break
        offset += len(batch)
    return pd.DataFrame(rows)


# ── ANALYSIS ──────────────────────────────────────────────────────────────────

def analyze_311(df):
    print(f"\n=== 311  (n={len(df)}) ===")
    env = df[df["service_name"].isin(ENVIRONMENTAL_311_TYPES)].copy()
    print(f"Env-related: {len(env)}")
    print("\nTop request types:\n", df["service_name"].value_counts().head(12).to_string())
    print("\nEnv breakdown:\n",     env["service_name"].value_counts().to_string())

    closed = env[env["status"] == "closed"].copy()
    if not closed.empty:
        closed["resolution_days"] = (
            closed["updated_datetime"] - closed["requested_datetime"]
        ).dt.total_seconds() / 86400
        print("\nAvg resolution (days):\n",
              closed.groupby("service_name")["resolution_days"].mean().sort_values().round(1).to_string())
    return env


def analyze_berdo(df):
    print(f"\n=== BERDO building emissions (n={len(df)}) ===")
    ghg_col  = next((c for c in df.columns if "natural gas emission" in c.lower()), None) \
             or next((c for c in df.columns if "ghg" in c.lower() or "emission" in c.lower()), None)
    type_col = next((c for c in df.columns if "largest property type" in c.lower()), None)
    zip_col  = next((c for c in df.columns if "zip" in c.lower()), None)
    if ghg_col:
        df[ghg_col] = pd.to_numeric(df[ghg_col], errors="coerce")
        print(f"\nTop emitting property types ({ghg_col}):\n",
              df.groupby(type_col)[ghg_col].mean().sort_values(ascending=False).head(8).round(0).to_string())
    return ghg_col, zip_col


def analyze_ghg_timeseries(df):
    print(f"\n=== Community GHG emissions (n={len(df)}) ===")
    print("Columns:", list(df.columns))
    year_col = next((c for c in df.columns if "year" in c.lower()), None)
    val_col  = next((c for c in df.columns
                     if any(k in c.lower() for k in ["emission", "ghg", "co2", "mtco2"])), None)
    if year_col and val_col:
        df[year_col] = pd.to_numeric(df[year_col], errors="coerce")
        df[val_col]  = pd.to_numeric(df[val_col],  errors="coerce")
        ts = df.groupby(year_col)[val_col].sum().sort_index()
        print("\nCommunity GHG by year (MTCO2e):\n", ts.round(0).to_string())
        return ts
    return None


def analyze_food_inspections(df):
    print(f"\n=== Food inspections (n={len(df)}) ===")
    print("Columns:", list(df.columns[:12]))
    date_col = next((c for c in df.columns if "date" in c.lower() or "dttm" in c.lower()), None)
    if date_col:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df["year_month"] = df[date_col].dt.to_period("M")
        monthly = df.groupby("year_month").size().sort_index()
        print(f"\nMonthly inspection counts (recent 6):\n{monthly.tail(6).to_string()}")
    return df, date_col


# ── CROSS-CORRELATION ─────────────────────────────────────────────────────────

def cross_correlate_zip(df_311_env, df_berdo, ghg_col, zip_col):
    """Spatial: BERDO emissions per zip vs 311 env complaints per zip."""
    print("\n=== SPATIAL CROSS-CORRELATION: emissions vs 311 complaints by ZIP ===")

    def extract_zip(addr):
        if pd.isna(addr):
            return None
        m = pd.Series(addr).str.extract(r"\b(0[12]\d{3})\b")
        return m.iloc[0, 0] if not m.empty and not pd.isna(m.iloc[0, 0]) else None

    df_311_env = df_311_env.copy()
    df_311_env["zip"] = df_311_env.get("address", pd.Series(dtype=str)).apply(extract_zip)

    complaints_by_zip = df_311_env.groupby("zip").size().rename("complaint_count")

    if zip_col and ghg_col:
        df_berdo[ghg_col] = pd.to_numeric(df_berdo[ghg_col], errors="coerce")
        # Normalise zip: strip spaces
        df_berdo["zip_clean"] = df_berdo[zip_col].astype(str).str.strip()
        df_311_env["zip_clean"] = df_311_env["zip"].astype(str).str.strip()

        emiss_by_zip = df_berdo.groupby("zip_clean")[ghg_col].sum().rename("total_emissions")
        comp_by_zip  = df_311_env.groupby("zip_clean").size().rename("complaint_count")
        merged = pd.concat([emiss_by_zip, comp_by_zip], axis=1).dropna()

        if len(merged) > 2:
            corr = merged["total_emissions"].corr(merged["complaint_count"])
            print(f"Pearson r (emissions vs complaints by zip): {corr:.3f}")
            print("\nZip-level table (sorted by emissions):\n",
                  merged.sort_values("total_emissions", ascending=False).head(10).round(0).to_string())
        return merged
    return pd.DataFrame()


def cross_correlate_time(df_311_env, df_food, food_date_col):
    """Temporal: monthly 311 env complaints vs food inspection violations."""
    print("\n=== TEMPORAL CROSS-CORRELATION: 311 env complaints vs food inspection volume ===")

    if df_311_env.empty:
        return

    df_311_env = df_311_env.copy()
    df_311_env["year_month"] = df_311_env["requested_datetime"].dt.to_period("M")
    monthly_311 = df_311_env.groupby("year_month").size().rename("env_complaints")

    if df_food is not None and not df_food.empty and food_date_col:
        df_food["year_month"] = pd.to_datetime(df_food[food_date_col], errors="coerce").dt.to_period("M")
        monthly_food = df_food.groupby("year_month").size().rename("inspections")
        combined = pd.concat([monthly_311, monthly_food], axis=1).dropna()
        if len(combined) > 2:
            corr = combined["env_complaints"].corr(combined["inspections"])
            print(f"Pearson r (monthly 311 env vs food inspections): {corr:.3f}")
            print("\nMonthly data (recent 12):\n",
                  combined.tail(12).to_string())
        return combined

    print("Monthly 311 env complaints:\n", monthly_311.tail(12).to_string())
    return monthly_311


# ── MAP ───────────────────────────────────────────────────────────────────────

def build_map(df_311, df_berdo, ghg_col, zip_col):
    try:
        import folium
        from folium.plugins import MarkerCluster, HeatMap
    except ImportError:
        print("Installing folium...")
        import subprocess, sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "folium", "-q"])
        import folium
        from folium.plugins import MarkerCluster, HeatMap

    print("\nBuilding interactive map...")
    m = folium.Map(location=[42.36, -71.06], zoom_start=12, tiles="CartoDB positron")

    # ── Layer 1: 311 complaints heatmap ──────────────────────────────────────
    map_cols = [c for c in ["lat", "long", "service_name"] if c in df_311.columns]
    heat_data = df_311[map_cols].dropna(subset=[c for c in ["lat", "long"] if c in map_cols]) if map_cols else pd.DataFrame()
    if not heat_data.empty and "lat" in heat_data and "long" in heat_data:
        HeatMap(
            heat_data[["lat", "long"]].values.tolist(),
            name="311 Complaints (heat)",
            radius=12, blur=15, min_opacity=0.3,
        ).add_to(m)

    # ── Layer 2: 311 env complaints clustered markers ────────────────────────
    if "service_name" not in df_311.columns or "lat" not in df_311.columns:
        env_df = pd.DataFrame()
    else:
        env_df = df_311[df_311["service_name"].isin(ENVIRONMENTAL_311_TYPES)].dropna(subset=["lat", "long"])
    color_map = {
        "Rodent Activity":                "#e74c3c",
        "Needle Cleanup":                 "#8e44ad",
        "Needle Pickup":                  "#8e44ad",
        "Litter":                         "#f39c12",
        "Overflowing Trash Can":          "#d35400",
        "Residential Trash out Illegally":"#c0392b",
        "Illegal Dumping":                "#922b21",
        "Graffiti Removal":               "#2980b9",
        "Encampment":                     "#27ae60",
        "Dead Animal Collection":         "#7f8c8d",
        "Dead Animal Pickup":             "#7f8c8d",
    }
    cluster = MarkerCluster(name="311 Env Complaints").add_to(m)
    for _, row in env_df.iterrows():
        color = color_map.get(row["service_name"], "#555")
        folium.CircleMarker(
            location=[row["lat"], row["long"]],
            radius=6, color=color, fill=True, fill_opacity=0.8,
            tooltip=f"{row['service_name']}<br>{row.get('address','')}<br>{str(row['requested_datetime'])[:10]}",
        ).add_to(cluster)

    # ── Layer 3: BERDO emissions bubbles by zip ───────────────────────────────
    ZIP_CENTROIDS = {
        "2101": (42.3601, -71.0549), "2102": (42.3601, -71.0549),
        "2108": (42.3588, -71.0707), "2109": (42.3601, -71.0518),
        "2110": (42.3572, -71.0504), "2111": (42.3630, -71.0596),
        "2112": (42.3658, -71.0614), "2113": (42.3658, -71.0560),
        "2114": (42.3614, -71.0706), "2115": (42.3447, -71.0875),
        "2116": (42.3508, -71.0757), "2118": (42.3361, -71.0714),
        "2119": (42.3199, -71.0875), "2120": (42.3308, -71.0931),
        "2121": (42.3077, -71.0864), "2122": (42.2993, -71.0691),
        "2124": (42.2935, -71.0797), "2125": (42.3077, -71.0704),
        "2126": (42.2801, -71.1003), "2127": (42.3341, -71.0376),
        "2128": (42.3700, -71.0155), "2129": (42.3839, -71.0576),
        "2130": (42.3096, -71.1148), "2131": (42.2906, -71.1249),
        "2132": (42.2786, -71.1519), "2134": (42.3508, -71.1102),
        "2135": (42.3536, -71.1502), "02215": (42.3468, -71.0997),
        "2210": (42.3361, -71.0360), "2467": (42.3219, -71.1826),
    }
    if ghg_col and zip_col and not df_berdo.empty:
        df_berdo[ghg_col] = pd.to_numeric(df_berdo[ghg_col], errors="coerce")
        emiss = df_berdo.groupby(zip_col)[ghg_col].sum().dropna()
        emiss_layer = folium.FeatureGroup(name="BERDO Emissions by ZIP (bubble)")
        max_val = emiss.max()
        for z, val in emiss.items():
            z_str = str(z).strip().lstrip("0")
            centroid = ZIP_CENTROIDS.get(z_str) or ZIP_CENTROIDS.get(str(z).strip())
            if centroid and max_val > 0:
                radius = 5 + 40 * (val / max_val)
                folium.CircleMarker(
                    location=centroid,
                    radius=radius,
                    color="#e74c3c", fill=True, fill_color="#e74c3c", fill_opacity=0.35,
                    tooltip=f"ZIP {z}<br>Natural gas emissions: {val:,.0f} kgCO2e",
                ).add_to(emiss_layer)
        emiss_layer.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)

    # Legend
    legend_html = """
    <div style="position:fixed;bottom:30px;left:30px;z-index:1000;background:white;
                padding:12px;border-radius:8px;border:1px solid #ccc;font-size:12px">
      <b>Boston 311 + Emissions</b><br>
      <span style="color:#e74c3c">&#9679;</span> Rodent Activity<br>
      <span style="color:#8e44ad">&#9679;</span> Needle Cleanup<br>
      <span style="color:#f39c12">&#9679;</span> Litter<br>
      <span style="color:#d35400">&#9679;</span> Overflowing Trash<br>
      <span style="color:#922b21">&#9679;</span> Illegal Dumping<br>
      <span style="color:#2980b9">&#9679;</span> Graffiti<br>
      <span style="color:#27ae60">&#9679;</span> Encampment<br>
      <span style="color:#e74c3c;opacity:0.4">&#11044;</span> BERDO emissions (bubble size = magnitude)
    </div>"""
    m.get_root().html.add_child(folium.Element(legend_html))

    out = "/home/user/ANEWFILES/boston_311_map.html"
    m.save(out)
    print(f"Map saved: {out}")
    return out


# ── CHARTS ────────────────────────────────────────────────────────────────────

def build_charts(df_311_env, zip_merged, time_merged):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        import subprocess, sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib", "-q"])
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Boston 311 Complaints × Pollution Data", fontsize=14, fontweight="bold")

    # 1. Bar: 311 env complaint types
    ax = axes[0, 0]
    if not df_311_env.empty:
        counts = df_311_env["service_name"].value_counts()
        counts.plot.barh(ax=ax, color="#e67e22")
        ax.set_title("311 Environmental Complaint Types")
        ax.set_xlabel("Count")

    # 2. Scatter: emissions vs complaints by ZIP
    ax = axes[0, 1]
    if not zip_merged.empty and len(zip_merged) > 1:
        ax.scatter(zip_merged["total_emissions"], zip_merged["complaint_count"],
                   alpha=0.7, color="#2980b9", edgecolors="white", s=80)
        for idx, row in zip_merged.iterrows():
            ax.annotate(str(idx), (row["total_emissions"], row["complaint_count"]),
                        fontsize=7, alpha=0.7)
        corr = zip_merged["total_emissions"].corr(zip_merged["complaint_count"])
        ax.set_title(f"Emissions vs 311 Complaints by ZIP  (r={corr:.2f})")
        ax.set_xlabel("Total Natural Gas Emissions (kgCO2e)")
        ax.set_ylabel("311 Env Complaint Count")

    # 3. Time series: monthly 311 complaints
    ax = axes[1, 0]
    if not df_311_env.empty and "requested_datetime" in df_311_env:
        monthly = df_311_env.groupby(df_311_env["requested_datetime"].dt.to_period("M")).size()
        monthly.index = monthly.index.to_timestamp()
        ax.plot(monthly.index, monthly.values, marker="o", color="#27ae60", linewidth=1.5)
        ax.fill_between(monthly.index, monthly.values, alpha=0.2, color="#27ae60")
        ax.set_title("Monthly 311 Env Complaints (2-year window)")
        ax.set_ylabel("Complaint Count")
        ax.tick_params(axis="x", rotation=30)

    # 4. Time series cross-correlation: 311 vs food inspections
    ax = axes[1, 1]
    if time_merged is not None and hasattr(time_merged, "columns") and "inspections" in time_merged.columns:
        t = time_merged.copy()
        t.index = t.index.to_timestamp()
        t["env_complaints_norm"] = (t["env_complaints"] - t["env_complaints"].mean()) / t["env_complaints"].std()
        t["inspections_norm"]    = (t["inspections"]    - t["inspections"].mean())    / t["inspections"].std()
        ax.plot(t.index, t["env_complaints_norm"], label="311 env (z-score)", color="#e74c3c")
        ax.plot(t.index, t["inspections_norm"],    label="Food inspections (z-score)", color="#3498db", linestyle="--")
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.legend(fontsize=8)
        ax.set_title("311 Env Complaints vs Food Inspections (normalised)")
        ax.tick_params(axis="x", rotation=30)

    plt.tight_layout()
    out = "/home/user/ANEWFILES/cross_correlation.png"
    plt.savefig(out, dpi=130, bbox_inches="tight")
    print(f"Chart saved: {out}")
    plt.close()


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    # 1. 311 (2-year window for better time-series)
    print("Fetching 311 (last 730 days)...")
    try:
        df_311     = fetch_311(days_back=730, max_records=2000)
        df_311_env = analyze_311(df_311)
    except Exception as e:
        print(f"311 failed: {e}")
        df_311 = df_311_env = pd.DataFrame()

    # 2. BERDO
    print("\nFetching BERDO building emissions...")
    try:
        df_berdo           = fetch_ckan("87521565-7f15-4b8d-a225-ac4df9e3f309", limit=3000)
        ghg_col, zip_col   = analyze_berdo(df_berdo)
    except Exception as e:
        print(f"BERDO failed: {e}")
        df_berdo = pd.DataFrame()
        ghg_col = zip_col = None

    # 3. Community GHG time-series
    print("\nFetching community GHG time-series...")
    try:
        df_ghg    = fetch_ckan("bd8dd4bb-867e-4ca2-b6c7-6c3bd9e6c290", limit=2000)
        ghg_ts    = analyze_ghg_timeseries(df_ghg)
    except Exception as e:
        print(f"Community GHG failed: {e}")
        ghg_ts = None

    # 4. Food inspections (sanitation proxy)
    print("\nFetching food establishment inspections...")
    try:
        df_food           = fetch_ckan("4582bec6-2b4f-4f9e-bc55-cbaa73117f4c", limit=5000)
        df_food, food_date = analyze_food_inspections(df_food)
    except Exception as e:
        print(f"Food inspections failed: {e}")
        df_food = None
        food_date = None

    # 5. Cross-correlations
    zip_merged  = cross_correlate_zip(df_311_env, df_berdo, ghg_col, zip_col)
    time_merged = cross_correlate_time(df_311_env, df_food, food_date)

    # 6. Map
    build_map(df_311, df_berdo, ghg_col, zip_col)

    # 7. Charts
    build_charts(df_311_env, zip_merged, time_merged)

    # 8. Save CSVs
    if not df_311_env.empty:
        df_311_env.to_csv("/home/user/ANEWFILES/311_environmental.csv", index=False)
    if not df_berdo.empty:
        df_berdo.to_csv("/home/user/ANEWFILES/berdo_emissions.csv", index=False)
    if not zip_merged.empty:
        zip_merged.to_csv("/home/user/ANEWFILES/zip_emissions_vs_complaints.csv")
    print("\nDone.")


if __name__ == "__main__":
    main()
