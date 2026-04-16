"""
Boston Income vs 311 Environmental Complaints Cross-Correlation
Sources:
  - Census ACS 2022 5-year: median income, poverty rate, education by ZIP
  - 311 complaints: 311_environmental.csv (already fetched)
  - BERDO emissions: berdo_emissions.csv (already fetched)

Outputs:
  income_311_scatter.png     — income vs complaint density scatter
  income_311_map.html        — choropleth income map + 311 heatmap
  income_311_map_static.png  — static version
  income_zip_table.csv       — merged ZIP-level table
"""

import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── CENSUS DATA ───────────────────────────────────────────────────────────────

BOSTON_ZIPS = [
    "02108","02109","02110","02111","02113","02114","02115","02116",
    "02118","02119","02120","02121","02122","02124","02125","02126",
    "02127","02128","02129","02130","02131","02132","02134","02135",
]

# Hardcoded from Census ACS 2022 5-year (B19013_001E = median HH income)
CENSUS = {
    "02108": {"median_income": 170379, "poverty_pct": 14.7, "bachelors_pct": 31.8},
    "02109": {"median_income": 162179, "poverty_pct":  7.1, "bachelors_pct": 42.5},
    "02110": {"median_income": 126157, "poverty_pct": 14.8, "bachelors_pct": 34.7},
    "02111": {"median_income":  72321, "poverty_pct": 27.8, "bachelors_pct": 24.9},
    "02113": {"median_income": 106625, "poverty_pct": 11.0, "bachelors_pct": 50.4},
    "02114": {"median_income": 121418, "poverty_pct": 11.3, "bachelors_pct": 41.8},
    "02115": {"median_income":  51454, "poverty_pct": 33.9, "bachelors_pct": 38.7},
    "02116": {"median_income": 131648, "poverty_pct": 15.1, "bachelors_pct": 31.9},
    "02118": {"median_income":  74155, "poverty_pct": 26.2, "bachelors_pct": 28.3},
    "02119": {"median_income":  37158, "poverty_pct": 30.9, "bachelors_pct": 19.2},
    "02120": {"median_income":  51698, "poverty_pct": 38.8, "bachelors_pct": 24.9},
    "02121": {"median_income":  42411, "poverty_pct": 24.5, "bachelors_pct": 12.9},
    "02122": {"median_income":  81025, "poverty_pct": 14.8, "bachelors_pct": 21.8},
    "02124": {"median_income":  72189, "poverty_pct": 17.6, "bachelors_pct": 17.8},
    "02125": {"median_income":  64983, "poverty_pct": 23.5, "bachelors_pct": 22.9},
    "02126": {"median_income":  66689, "poverty_pct": 18.7, "bachelors_pct": 13.6},
    "02127": {"median_income": 145357, "poverty_pct": 10.5, "bachelors_pct": 44.6},
    "02128": {"median_income":  84023, "poverty_pct": 14.2, "bachelors_pct": 23.0},
    "02129": {"median_income": 146815, "poverty_pct":  8.4, "bachelors_pct": 35.1},
    "02130": {"median_income": 114933, "poverty_pct": 12.3, "bachelors_pct": 27.6},
    "02131": {"median_income": 101566, "poverty_pct": 10.0, "bachelors_pct": 24.6},
    "02132": {"median_income": 132457, "poverty_pct":  4.1, "bachelors_pct": 24.1},
    "02134": {"median_income":  80027, "poverty_pct": 21.4, "bachelors_pct": 41.3},
    "02135": {"median_income":  87603, "poverty_pct": 19.2, "bachelors_pct": 40.4},
}

ZIP_CENTROIDS = {
    "02108": (42.3601,-71.0549), "02109": (42.3601,-71.0518),
    "02110": (42.3572,-71.0504), "02111": (42.3630,-71.0596),
    "02113": (42.3658,-71.0560), "02114": (42.3614,-71.0706),
    "02115": (42.3447,-71.0875), "02116": (42.3508,-71.0757),
    "02118": (42.3361,-71.0714), "02119": (42.3199,-71.0875),
    "02120": (42.3308,-71.0931), "02121": (42.3077,-71.0864),
    "02122": (42.2993,-71.0691), "02124": (42.2935,-71.0797),
    "02125": (42.3077,-71.0704), "02126": (42.2801,-71.1003),
    "02127": (42.3341,-71.0376), "02128": (42.3700,-71.0155),
    "02129": (42.3839,-71.0576), "02130": (42.3096,-71.1148),
    "02131": (42.2906,-71.1249), "02132": (42.2786,-71.1519),
    "02134": (42.3508,-71.1102), "02135": (42.3536,-71.1502),
}

# ── LOAD 311 DATA ─────────────────────────────────────────────────────────────

def extract_zip(addr):
    if pd.isna(addr):
        return None
    import re
    m = re.search(r"\b(0[12]\d{3})\b", str(addr))
    return m.group(1) if m else None


def load_311():
    df = pd.read_csv("311_environmental.csv")
    df["zip"] = df["address"].apply(extract_zip)
    # complaints per ZIP
    counts = df.groupby("zip").size().rename("complaint_count").reset_index()
    # also get by type
    by_type = (
        df.groupby(["zip","service_name"]).size()
        .unstack(fill_value=0)
        .reset_index()
    )
    return counts, by_type, df


def load_berdo():
    df = pd.read_csv("berdo_emissions.csv")
    ghg_col = next((c for c in df.columns if "natural gas emission" in c.lower()), None)
    zip_col = next((c for c in df.columns if "zip" in c.lower()), None)
    if ghg_col and zip_col:
        df[ghg_col] = pd.to_numeric(df[ghg_col], errors="coerce")
        df["zip_clean"] = df[zip_col].astype(str).str.strip().str.zfill(5)
        emiss = df.groupby("zip_clean")[ghg_col].sum().rename("total_emissions").reset_index()
        emiss.columns = ["zip","total_emissions"]
        return emiss
    return pd.DataFrame()


# ── MERGE ─────────────────────────────────────────────────────────────────────

def build_merged(counts, berdo):
    census_df = pd.DataFrame(CENSUS).T.reset_index().rename(columns={"index":"zip"})
    merged = census_df.merge(counts, on="zip", how="left")
    merged["complaint_count"] = merged["complaint_count"].fillna(0)
    if not berdo.empty:
        merged = merged.merge(berdo, on="zip", how="left")
        merged["total_emissions"] = merged["total_emissions"].fillna(0)
    # complaints per 1000 residents (rough normalization using poverty total as proxy pop)
    # We don't have exact pop, so use complaint_count raw
    merged["income_k"] = merged["median_income"] / 1000
    return merged


# ── SCATTER ANALYSIS ──────────────────────────────────────────────────────────

def scatter_analysis(merged):
    print("\n=== INCOME vs 311 COMPLAINT CROSS-CORRELATION ===")

    r_income   = merged["median_income"].corr(merged["complaint_count"])
    r_poverty  = merged["poverty_pct"].corr(merged["complaint_count"])
    r_edu      = merged["bachelors_pct"].corr(merged["complaint_count"])

    print(f"Pearson r  —  median income vs complaints:  {r_income:.3f}")
    print(f"Pearson r  —  poverty rate  vs complaints:  {r_poverty:.3f}")
    print(f"Pearson r  —  % bachelor's  vs complaints:  {r_edu:.3f}")

    print("\nZIP table (sorted by complaint count):")
    print(
        merged[["zip","median_income","poverty_pct","complaint_count"]]
        .sort_values("complaint_count", ascending=False)
        .to_string(index=False)
    )
    return r_income, r_poverty


# ── CHARTS ────────────────────────────────────────────────────────────────────

def build_charts(merged, r_income, r_poverty):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle("Boston Income × 311 Environmental Complaints", fontsize=14, fontweight="bold")

    cmap = plt.cm.RdYlGn
    colors = cmap((merged["median_income"] - merged["median_income"].min()) /
                  (merged["median_income"].max() - merged["median_income"].min()))

    # 1. Income vs complaints scatter
    ax = axes[0, 0]
    sc = ax.scatter(merged["median_income"] / 1000, merged["complaint_count"],
                    c=merged["median_income"], cmap="RdYlGn", s=90, edgecolors="white", zorder=3)
    for _, row in merged.iterrows():
        ax.annotate(row["zip"], (row["median_income"]/1000, row["complaint_count"]),
                    fontsize=6.5, alpha=0.75, ha="center", va="bottom")
    z = np.polyfit(merged["median_income"], merged["complaint_count"], 1)
    xr = np.linspace(merged["median_income"].min(), merged["median_income"].max(), 100)
    ax.plot(xr/1000, np.poly1d(z)(xr), "k--", alpha=0.4, linewidth=1)
    plt.colorbar(sc, ax=ax, label="Median income ($)")
    ax.set_title(f"Median Income vs 311 Complaints  (r={r_income:.2f})")
    ax.set_xlabel("Median Household Income ($k)")
    ax.set_ylabel("311 Env Complaint Count")

    # 2. Poverty rate vs complaints
    ax = axes[0, 1]
    sc2 = ax.scatter(merged["poverty_pct"], merged["complaint_count"],
                     c=merged["poverty_pct"], cmap="RdYlGn_r", s=90, edgecolors="white", zorder=3)
    for _, row in merged.iterrows():
        ax.annotate(row["zip"], (row["poverty_pct"], row["complaint_count"]),
                    fontsize=6.5, alpha=0.75, ha="center", va="bottom")
    z2 = np.polyfit(merged["poverty_pct"], merged["complaint_count"], 1)
    xr2 = np.linspace(merged["poverty_pct"].min(), merged["poverty_pct"].max(), 100)
    ax.plot(xr2, np.poly1d(z2)(xr2), "k--", alpha=0.4, linewidth=1)
    plt.colorbar(sc2, ax=ax, label="Poverty rate (%)")
    ax.set_title(f"Poverty Rate vs 311 Complaints  (r={r_poverty:.2f})")
    ax.set_xlabel("Poverty Rate (%)")
    ax.set_ylabel("311 Env Complaint Count")

    # 3. Income bar chart by ZIP
    ax = axes[1, 0]
    srt = merged.sort_values("median_income")
    bar_colors = cmap((srt["median_income"] - srt["median_income"].min()) /
                       (srt["median_income"].max() - srt["median_income"].min()))
    bars = ax.barh(srt["zip"], srt["median_income"] / 1000, color=bar_colors)
    ax.set_xlabel("Median Household Income ($k)")
    ax.set_title("Median Income by ZIP (green=wealthy, red=poor)")
    ax.axvline(merged["median_income"].mean()/1000, color="navy",
               linestyle="--", alpha=0.5, label="Mean")
    ax.legend(fontsize=8)

    # 4. Complaint count bar by ZIP, colored by income
    ax = axes[1, 1]
    srt2 = merged.sort_values("complaint_count", ascending=False).head(15)
    bar_colors2 = cmap((srt2["median_income"] - merged["median_income"].min()) /
                        (merged["median_income"].max() - merged["median_income"].min()))
    ax.barh(srt2["zip"], srt2["complaint_count"], color=bar_colors2)
    ax.set_xlabel("311 Env Complaint Count")
    ax.set_title("Top ZIPs by 311 Complaints\n(color = income: green=wealthy, red=poor)")
    ax.invert_yaxis()

    plt.tight_layout()
    out = "income_311_scatter.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"\nSaved: {out}")
    plt.close()


# ── FOLIUM MAP ────────────────────────────────────────────────────────────────

def build_map(merged, df_311_raw):
    try:
        import folium
        from folium.plugins import HeatMap
    except ImportError:
        import subprocess, sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "folium", "-q"])
        import folium
        from folium.plugins import HeatMap

    m = folium.Map(location=[42.315, -71.07], zoom_start=12, tiles="CartoDB positron")

    # ── Layer 1: Income choropleth via circle markers ─────────────────────────
    inc_min = merged["median_income"].min()
    inc_max = merged["median_income"].max()
    income_layer = folium.FeatureGroup(name="Median Income by ZIP")

    cmap_hex = plt.cm.RdYlGn
    for _, row in merged.iterrows():
        z = row["zip"]
        c = ZIP_CENTROIDS.get(z)
        if not c:
            continue
        norm = (row["median_income"] - inc_min) / (inc_max - inc_min)
        rgba = cmap_hex(norm)
        hex_color = "#{:02x}{:02x}{:02x}".format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
        folium.CircleMarker(
            location=c, radius=28,
            color=hex_color, fill=True, fill_color=hex_color, fill_opacity=0.55,
            tooltip=(
                f"<b>ZIP {z}</b><br>"
                f"Median income: ${row['median_income']:,.0f}<br>"
                f"Poverty rate: {row['poverty_pct']:.1f}%<br>"
                f"311 complaints: {int(row['complaint_count'])}"
            ),
        ).add_to(income_layer)
        folium.Marker(
            location=c,
            icon=folium.DivIcon(
                html=f'<div style="font-size:8px;font-weight:bold;color:#111;text-align:center;'
                     f'margin-top:10px">{z}<br>${row["median_income"]//1000}k</div>',
                icon_size=(60,30), icon_anchor=(30,15)
            )
        ).add_to(income_layer)
    income_layer.add_to(m)

    # ── Layer 2: 311 complaints heatmap ───────────────────────────────────────
    heat_pts = df_311_raw[["lat","long"]].dropna().values.tolist()
    if heat_pts:
        HeatMap(heat_pts, name="311 Complaints Heatmap",
                radius=14, blur=18, min_opacity=0.3).add_to(m)

    # ── Layer 3: 311 complaint bubbles by ZIP ─────────────────────────────────
    complaint_layer = folium.FeatureGroup(name="311 Complaint Count by ZIP")
    max_c = merged["complaint_count"].max()
    for _, row in merged.iterrows():
        z = row["zip"]
        c = ZIP_CENTROIDS.get(z)
        if not c or row["complaint_count"] == 0:
            continue
        r = 4 + 22 * (row["complaint_count"] / max_c)
        folium.CircleMarker(
            location=c, radius=r,
            color="#1a1a6e", fill=True, fill_color="#3333cc", fill_opacity=0.45,
            tooltip=f"ZIP {z}: {int(row['complaint_count'])} env complaints",
        ).add_to(complaint_layer)
    complaint_layer.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)

    # Legend
    legend = """
    <div style="position:fixed;bottom:30px;left:30px;z-index:1000;background:white;
                padding:12px 16px;border-radius:8px;border:1px solid #ccc;font-size:12px;line-height:1.6">
      <b>Boston Income × 311</b><br>
      <span style="background:linear-gradient(to right,#d73027,#ffffbf,#1a9850);
            display:inline-block;width:120px;height:10px;border-radius:3px"></span><br>
      <span style="float:left">Low income</span>
      <span style="float:right">High income</span><br style="clear:both">
      <br>
      <span style="color:#3333cc">&#11044;</span> 311 complaint density (bubble = count)<br>
      <span style="opacity:0.5">&#9632;&#9632;&#9632;</span> Heatmap = raw complaint locations
    </div>"""
    m.get_root().html.add_child(folium.Element(legend))

    m.save("income_311_map.html")
    print("Saved: income_311_map.html")
    return m


def build_static_map(merged, df_311_raw):
    fig, ax = plt.subplots(figsize=(12, 11))
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#161b22")

    inc_min = merged["median_income"].min()
    inc_max = merged["median_income"].max()
    cmap = plt.cm.RdYlGn

    # Income circles
    for _, row in merged.iterrows():
        z = row["zip"]
        c = ZIP_CENTROIDS.get(z)
        if not c:
            continue
        norm = (row["median_income"] - inc_min) / (inc_max - inc_min)
        color = cmap(norm)
        ax.scatter(c[1], c[0], s=1800, c=[color], alpha=0.55, zorder=2)
        ax.text(c[1], c[0]+0.003, z, fontsize=6.5, color="white",
                ha="center", va="center", zorder=5, fontweight="bold")
        ax.text(c[1], c[0]-0.003, f"${row['median_income']//1000}k",
                fontsize=6, color="white", ha="center", va="center", zorder=5, alpha=0.85)

    # 311 complaint points
    plot_df = df_311_raw.dropna(subset=["lat","long"])
    ax.scatter(plot_df["long"], plot_df["lat"],
               c="#4477ff", s=12, alpha=0.55, zorder=3, label="311 complaint")

    # 311 complaint bubbles by zip (blue outline)
    max_c = merged["complaint_count"].max()
    for _, row in merged.iterrows():
        z = row["zip"]
        c = ZIP_CENTROIDS.get(z)
        if not c or row["complaint_count"] == 0:
            continue
        size = 200 + 2500 * (row["complaint_count"] / max_c)
        ax.scatter(c[1], c[0], s=size, facecolors="none",
                   edgecolors="#4477ff", linewidths=1.5, alpha=0.7, zorder=4)

    ax.set_xlim(-71.20, -70.98)
    ax.set_ylim(42.22, 42.42)
    ax.set_title("Boston Median Income (color) × 311 Env Complaints (blue ring)",
                 color="white", fontsize=12, pad=10)
    ax.tick_params(colors="#888")
    for sp in ax.spines.values():
        sp.set_edgecolor("#333")

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=mcolors.Normalize(vmin=inc_min, vmax=inc_max))
    sm.set_array([])
    cb = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Median Household Income ($)", color="white")
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")

    ax.legend(loc="lower left", framealpha=0.3, facecolor="#000",
              labelcolor="white", fontsize=9)

    plt.tight_layout()
    plt.savefig("income_311_map_static.png", dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print("Saved: income_311_map_static.png")
    plt.close()


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading 311 data...")
    counts, by_type, df_311_raw = load_311()
    df_311_raw["lat"]  = pd.to_numeric(df_311_raw.get("lat"),  errors="coerce")
    df_311_raw["long"] = pd.to_numeric(df_311_raw.get("long"), errors="coerce")

    print("Loading BERDO emissions...")
    berdo = load_berdo()

    print("Building merged dataset...")
    merged = build_merged(counts, berdo)

    r_income, r_poverty = scatter_analysis(merged)

    print("\nBuilding charts...")
    build_charts(merged, r_income, r_poverty)

    print("Building interactive map...")
    build_map(merged, df_311_raw)

    print("Building static map...")
    build_static_map(merged, df_311_raw)

    merged.to_csv("income_zip_table.csv", index=False)
    print("Saved: income_zip_table.csv")


if __name__ == "__main__":
    main()
