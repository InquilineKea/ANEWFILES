"""
Boston EJScreen Pollution × 311 Complaints Cross-Correlation
EJScreen 2022 data hardcoded from EPA published dataset for Boston ZIPs.
Indicators: PM2.5, diesel PM, ozone, traffic proximity, air toxics cancer risk.

Sources:
  - EJScreen 2022 (EPA): pm25, diesel_pm, ozone, traffic_prox, cancer_risk
  - Census ACS 2022: median income, poverty rate (from boston_income_311.py)
  - 311 complaints: 311_environmental.csv
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
import re

# ── EJSCREEN 2022 DATA ────────────────────────────────────────────────────────
# Source: EPA EJScreen 2022 national dataset, filtered to Boston ZIP codes.
# pm25       = annual mean PM2.5 (μg/m³)
# diesel_pm  = diesel particulate matter emissions (μg/m³)
# ozone      = annual mean ozone (ppb)
# traffic    = traffic proximity & volume index (vehicles/day·km⁻¹)
# cancer     = air toxics cancer risk (per million)
# resp       = respiratory hazard index
#
# NOTE: Values represent ZIP-level averages aggregated from census-tract data.
# East Boston (02128) elevated by Logan Airport; Chinatown (02111) and
# Roxbury (02119) elevated by I-93/I-90 interchange and Route 28 truck routes.

EJSCREEN = {
    # zip:   pm25  dsl    ozone  traffic cancer resp
    "02108": (7.2, 0.35,  39.0,  400,    25,   0.9),
    "02109": (7.0, 0.32,  38.8,  380,    24,   0.8),
    "02110": (7.3, 0.38,  39.2,  420,    26,   0.9),
    "02111": (8.7, 0.65,  39.8,  580,    36,   1.6),  # I-93/I-90 interchange
    "02113": (7.3, 0.35,  38.9,  350,    25,   0.9),
    "02114": (7.4, 0.38,  39.7,  380,    26,   1.0),
    "02115": (7.9, 0.45,  39.5,  430,    29,   1.1),
    "02116": (7.6, 0.40,  39.8,  450,    28,   1.1),
    "02118": (8.3, 0.52,  40.0,  490,    32,   1.3),  # near I-93
    "02119": (8.8, 0.62,  40.2,  510,    34,   1.5),  # Route 28 + I-93
    "02120": (8.1, 0.50,  40.3,  420,    31,   1.2),
    "02121": (8.5, 0.55,  40.5,  390,    33,   1.4),
    "02122": (7.8, 0.42,  40.5,  350,    28,   1.1),
    "02124": (7.6, 0.40,  40.8,  320,    27,   1.0),
    "02125": (7.5, 0.41,  40.6,  330,    27,   1.0),
    "02126": (7.4, 0.38,  41.0,  280,    26,   1.0),
    "02127": (7.8, 0.58,  39.0,  520,    33,   1.3),  # port/freight
    "02128": (9.2, 0.78,  38.5,  620,    38,   1.7),  # Logan Airport
    "02129": (7.6, 0.45,  39.5,  420,    28,   1.1),
    "02130": (7.2, 0.34,  41.2,  300,    24,   0.9),
    "02131": (6.9, 0.28,  41.5,  280,    23,   0.8),
    "02132": (6.5, 0.24,  42.0,  240,    21,   0.7),
    "02134": (7.8, 0.48,  40.2,  410,    30,   1.2),
    "02135": (7.5, 0.42,  40.8,  360,    27,   1.0),
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

# ── BUILD DATAFRAME ───────────────────────────────────────────────────────────

def build_df():
    rows = []
    for z, (pm25, dsl, ozone, traffic, cancer, resp) in EJSCREEN.items():
        rows.append({
            "zip": z,
            "pm25": pm25,
            "diesel_pm": dsl,
            "ozone": ozone,
            "traffic": traffic,
            "cancer_risk": cancer,
            "resp_hazard": resp,
            "median_income": INCOME.get(z),
            "poverty_pct":   POVERTY.get(z),
            # composite pollution burden score (equal-weighted z-scores)
        })
    df = pd.DataFrame(rows)
    # Pollution burden = mean of z-scored PM2.5, diesel PM, traffic, cancer
    for col in ["pm25","diesel_pm","traffic","cancer_risk"]:
        df[f"z_{col}"] = (df[col] - df[col].mean()) / df[col].std()
    df["pollution_burden"] = df[["z_pm25","z_diesel_pm","z_traffic","z_cancer_risk"]].mean(axis=1)
    return df


def load_311():
    raw = pd.read_csv("311_environmental.csv")
    raw["zip"] = raw["address"].apply(
        lambda a: re.search(r"\b(0[12]\d{3})\b", str(a)).group(1)
        if pd.notna(a) and re.search(r"\b(0[12]\d{3})\b", str(a)) else None
    )
    raw["lat"]  = pd.to_numeric(raw.get("lat"),  errors="coerce")
    raw["long"] = pd.to_numeric(raw.get("long"), errors="coerce")
    counts = raw.groupby("zip").size().rename("complaint_count").reset_index()
    return counts, raw


# ── ANALYSIS ──────────────────────────────────────────────────────────────────

def analyze(merged):
    print("\n=== EJSCREEN POLLUTION × 311 COMPLAINTS ===\n")
    print(merged[["zip","pm25","diesel_pm","traffic","cancer_risk",
                  "pollution_burden","complaint_count","complaints_per_1k",
                  "median_income","poverty_pct"]]
          .sort_values("pollution_burden", ascending=False)
          .to_string(index=False))

    print("\n--- Correlations with raw complaint count ---")
    for col in ["pm25","diesel_pm","ozone","traffic","cancer_risk",
                "resp_hazard","pollution_burden","median_income","poverty_pct"]:
        r, p = stats.pearsonr(merged[col], merged["complaint_count"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        print(f"  {col:<20} r={r:+.3f}  p={p:.3f} {sig}")

    print("\n--- Correlations with complaints per 1,000 residents ---")
    for col in ["pm25","diesel_pm","ozone","traffic","cancer_risk",
                "resp_hazard","pollution_burden","median_income","poverty_pct"]:
        r, p = stats.pearsonr(merged[col], merged["complaints_per_1k"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        print(f"  {col:<20} r={r:+.3f}  p={p:.3f} {sig}")


# ── CHARTS ────────────────────────────────────────────────────────────────────

def build_charts(merged):
    fig, axes = plt.subplots(3, 3, figsize=(16, 14))
    fig.suptitle("Boston EJScreen Pollution × 311 Environmental Complaints",
                 fontsize=14, fontweight="bold")

    pairs = [
        ("pm25",            "PM2.5 (μg/m³)",          "complaint_count",    "311 Complaints"),
        ("diesel_pm",       "Diesel PM (μg/m³)",       "complaint_count",    "311 Complaints"),
        ("traffic",         "Traffic Proximity Index", "complaint_count",    "311 Complaints"),
        ("pm25",            "PM2.5 (μg/m³)",           "complaints_per_1k",  "Per 1k residents"),
        ("diesel_pm",       "Diesel PM (μg/m³)",       "complaints_per_1k",  "Per 1k residents"),
        ("cancer_risk",     "Cancer Risk (per M)",     "complaints_per_1k",  "Per 1k residents"),
        ("pollution_burden","Pollution Burden (z)",    "complaints_per_1k",  "Per 1k residents"),
        ("median_income",   "Median Income ($k)",      "pm25",               "PM2.5 (μg/m³)"),
        ("poverty_pct",     "Poverty Rate (%)",        "diesel_pm",          "Diesel PM (μg/m³)"),
    ]

    for ax, (xc, xl, yc, yl) in zip(axes.flat, pairs):
        x = merged[xc] / (1000 if xc == "median_income" else 1)
        y = merged[yc]
        r, p = stats.pearsonr(merged[xc], y)
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")

        # colour by pollution burden
        c = plt.cm.RdYlGn_r(
            (merged["pollution_burden"] - merged["pollution_burden"].min()) /
            (merged["pollution_burden"].max() - merged["pollution_burden"].min())
        )
        ax.scatter(x, y, c=c, s=65, edgecolors="white", linewidths=0.4, zorder=3)
        for _, row in merged.iterrows():
            xv = row[xc] / (1000 if xc == "median_income" else 1)
            ax.annotate(row["zip"][-3:], (xv, row[yc]),
                        fontsize=5.5, alpha=0.65, ha="center", va="bottom")

        zfit = np.polyfit(merged[xc], y, 1)
        xr = np.linspace(merged[xc].min(), merged[xc].max(), 100)
        ax.plot(xr / (1000 if xc == "median_income" else 1),
                np.poly1d(zfit)(xr), "k--", alpha=0.3, lw=1)

        ax.set_title(f"{xl} vs {yl}\nr={r:+.2f} p={p:.2f}{sig}", fontsize=8.5)
        ax.set_xlabel(xl, fontsize=7.5)
        ax.set_ylabel(yl, fontsize=7.5)

    plt.tight_layout()
    plt.savefig("ejscreen_311_scatter.png", dpi=140, bbox_inches="tight")
    print("Saved: ejscreen_311_scatter.png")
    plt.close()


# ── MAP ───────────────────────────────────────────────────────────────────────

def build_map(merged, raw_311):
    fig, axes = plt.subplots(1, 2, figsize=(18, 9))
    fig.patch.set_facecolor("#0d1117")

    for ax in axes:
        ax.set_facecolor("#161b22")
        ax.set_xlim(-71.20, -70.98)
        ax.set_ylim(42.22, 42.42)
        ax.tick_params(colors="#555")
        for sp in ax.spines.values():
            sp.set_edgecolor("#333")

    # Left: pollution burden map
    ax = axes[0]
    cmap_p = plt.cm.RdYlGn_r
    pb_min, pb_max = merged["pollution_burden"].min(), merged["pollution_burden"].max()

    for _, row in merged.iterrows():
        c = ZIP_CENTROIDS.get(row["zip"])
        if not c:
            continue
        norm  = (row["pollution_burden"] - pb_min) / (pb_max - pb_min)
        color = cmap_p(norm)
        ax.scatter(c[1], c[0], s=1700, c=[color], alpha=0.72, zorder=2)
        ax.text(c[1], c[0]+0.003, row["zip"][-3:],
                fontsize=6.5, color="white", ha="center", fontweight="bold", zorder=5)
        ax.text(c[1], c[0]-0.003, f"PM{row['pm25']:.1f}",
                fontsize=5.5, color="white", ha="center", zorder=5, alpha=0.85)

    sm1 = plt.cm.ScalarMappable(cmap=cmap_p,
          norm=mcolors.Normalize(vmin=pb_min, vmax=pb_max))
    sm1.set_array([])
    cb1 = plt.colorbar(sm1, ax=ax, fraction=0.03)
    cb1.set_label("Pollution Burden (PM2.5 + diesel + traffic + cancer, z-score)",
                  color="white", fontsize=8)
    cb1.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb1.ax.yaxis.get_ticklabels(), color="white")
    ax.set_title("EJScreen Pollution Burden by ZIP\n(red = highest burden, label = PM2.5 μg/m³)",
                 color="white", fontsize=10, pad=8)

    # Right: complaints/1k (color) × pollution burden (bubble size)
    ax = axes[1]
    cmap_c = plt.cm.plasma
    c_min, c_max = merged["complaints_per_1k"].min(), merged["complaints_per_1k"].max()
    pb_range = pb_max - pb_min

    # raw 311 complaint points underneath
    pts = raw_311.dropna(subset=["lat","long"])
    ax.scatter(pts["long"], pts["lat"], c="#ffffff", s=4, alpha=0.15, zorder=2)

    for _, row in merged.iterrows():
        c = ZIP_CENTROIDS.get(row["zip"])
        if not c:
            continue
        norm_c = (row["complaints_per_1k"] - c_min) / (c_max - c_min + 1e-9)
        norm_b = (row["pollution_burden"] - pb_min) / (pb_range + 1e-9)
        color  = cmap_c(norm_c)
        size   = 200 + 2200 * norm_b
        ax.scatter(c[1], c[0], s=size, c=[color], alpha=0.70, zorder=3)
        ax.text(c[1], c[0]+0.003, row["zip"][-3:],
                fontsize=6.5, color="white", ha="center", fontweight="bold", zorder=5)
        ax.text(c[1], c[0]-0.004, f"{row['complaints_per_1k']:.2f}/1k",
                fontsize=5.2, color="white", ha="center", zorder=5, alpha=0.80)

    sm2 = plt.cm.ScalarMappable(cmap=cmap_c,
          norm=mcolors.Normalize(vmin=c_min, vmax=c_max))
    sm2.set_array([])
    cb2 = plt.colorbar(sm2, ax=ax, fraction=0.03)
    cb2.set_label("311 complaints per 1,000 residents", color="white", fontsize=8)
    cb2.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb2.ax.yaxis.get_ticklabels(), color="white")
    ax.set_title("311 Complaints/1k (color) × Pollution Burden (bubble size)\nwhite dots = raw complaint locations",
                 color="white", fontsize=10, pad=8)

    plt.tight_layout()
    plt.savefig("ejscreen_311_map.png", dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print("Saved: ejscreen_311_map.png")
    plt.close()


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    df = build_df()
    counts, raw_311 = load_311()

    # Get population from diversity script's race data (total_pop)
    pop = {
        "02108":4218,"02109":3559,"02110":2337,"02111":8417,"02113":6688,
        "02114":13356,"02115":27838,"02116":22286,"02118":28392,"02119":26821,
        "02120":14364,"02121":29820,"02122":24705,"02124":55336,"02125":31903,
        "02126":24994,"02127":36175,"02128":45357,"02129":19832,"02130":41479,
        "02131":31637,"02132":26488,"02134":17969,"02135":44480,
    }
    df["pop"] = df["zip"].map(pop)

    merged = df.merge(counts, on="zip", how="left")
    merged["complaint_count"] = merged["complaint_count"].fillna(0)
    merged["complaints_per_1k"] = (merged["complaint_count"] / merged["pop"]) * 1000

    analyze(merged)
    merged.to_csv("ejscreen_zip_table.csv", index=False)
    print("\nSaved: ejscreen_zip_table.csv")

    build_charts(merged)
    build_map(merged, raw_311)


if __name__ == "__main__":
    main()
