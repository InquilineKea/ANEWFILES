"""
Boston Racial Diversity × 311 Complaints
Computes alpha diversity (Shannon, Simpson, Richness, Evenness),
beta diversity (Bray-Curtis dissimilarity matrix), and cross-correlates
with 311 environmental complaint counts by ZIP code.

Census ACS 2022 5-year race data (B02001):
  _002E White, _003E Black, _004E Native, _005E Asian,
  _006E Pacific Islander, _007E Other, _008E Two or more
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
from scipy.spatial.distance import braycurtis
import re

# ── RACE DATA (Census ACS 2022 5-year) ───────────────────────────────────────

RACE = {
    #  zip:  [total, white, black, native, asian, pac_isl, other, two_plus]
    "02108": [4218,  3598,  176,   0,    107,   0,    186,   151],
    "02109": [3559,  2893,   42,   0,    381,   0,     11,   232],
    "02110": [2337,  2000,  129,   0,    103,  13,     20,    72],
    "02111": [8417,  3347,  202,  27,   3958,   0,    153,   730],
    "02113": [6688,  5977,  238,   8,    213,   0,    179,    73],
    "02114": [13356,10145,  413,  13,   1158,   0,    555,  1072],
    "02115": [27838,16024, 2595,  28,   5302,  88,   1178,  2623],
    "02116": [22286,15072, 1150,  45,   3442, 126,    695,  1756],
    "02118": [28392,13979, 4843, 229,   4452,   0,   2506,  2383],
    "02119": [26821, 4972,13328, 148,    721,   0,   3948,  3704],
    "02120": [14364, 5850, 2892, 137,   2525,   4,   1880,  1076],
    "02121": [29820, 3386,16227,  27,    203,   0,   4439,  5538],
    "02122": [24705, 8681, 7074, 132,   3756,   0,   1988,  3074],
    "02124": [55336,11458,27518, 113,   3672,   7,   3740,  8828],
    "02125": [31903,11150, 9656, 292,   3508,   0,   3507,  3790],
    "02126": [24994, 1838,18306,   9,    289,   0,   1132,  3420],
    "02127": [36175,28507, 1679,  32,   2527,   0,   1119,  2311],
    "02128": [45357,21623, 1437, 284,   1844,  39,   3249, 16881],
    "02129": [19832,14998, 1394,   2,   1602,   0,    397,  1439],
    "02130": [41479,26372, 5282, 219,   2640,  19,   3010,  3937],
    "02131": [31637,17556, 7669,  62,    684,   0,   2542,  3124],
    "02132": [26488,19931, 1309,  78,   2302,   0,    769,  2099],
    "02134": [17969,10829,  985, 101,   4129,   0,    737,  1188],
    "02135": [44480,28348, 2461, 149,   7546,  80,   2873,  3023],
}
RACE_COLS = ["white","black","native","asian","pac_isl","other","two_plus"]

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

# ── BUILD RACE DATAFRAME ──────────────────────────────────────────────────────

def build_race_df():
    rows = []
    for z, vals in RACE.items():
        total = vals[0]
        counts = np.array(vals[1:], dtype=float)
        props  = counts / total if total > 0 else counts

        # Alpha diversity metrics
        # Shannon: H = -sum(p * ln(p)), ignore zeros
        p_nz = props[props > 0]
        shannon = -np.sum(p_nz * np.log(p_nz))

        # Simpson: D = 1 - sum(p^2)
        simpson = 1 - np.sum(props ** 2)

        # Richness: number of groups with >1% representation
        richness = int(np.sum(props > 0.01))

        # Evenness (Pielou's J): H / ln(S), S = richness
        evenness = shannon / np.log(richness) if richness > 1 else 0.0

        # Dominance: proportion of largest group
        dominance = props.max()

        rows.append({
            "zip": z,
            "total_pop": total,
            "shannon":   round(shannon, 4),
            "simpson":   round(simpson, 4),
            "richness":  richness,
            "evenness":  round(evenness, 4),
            "dominance": round(dominance, 4),
            "pct_white": round(props[0]*100, 1),
            "pct_black": round(props[1]*100, 1),
            "pct_asian": round(props[3]*100, 1),
            "pct_hispanic_other": round((props[5]+props[4])*100, 1),
            **{f"prop_{c}": round(p, 4) for c, p in zip(RACE_COLS, props)},
        })
    return pd.DataFrame(rows)


# ── BETA DIVERSITY (Bray-Curtis) ──────────────────────────────────────────────

def beta_diversity(df):
    """Compute pairwise Bray-Curtis dissimilarity between all ZIP pairs."""
    prop_cols = [f"prop_{c}" for c in RACE_COLS]
    zips  = df["zip"].tolist()
    mat   = df[prop_cols].values
    n     = len(zips)
    bc    = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = braycurtis(mat[i], mat[j])
            bc[i,j] = bc[j,i] = d
    return pd.DataFrame(bc, index=zips, columns=zips)


# ── LOAD 311 DATA ─────────────────────────────────────────────────────────────

def load_311():
    df = pd.read_csv("311_environmental.csv")
    df["zip"] = df["address"].apply(
        lambda a: re.search(r"\b(0[12]\d{3})\b", str(a)).group(1)
        if pd.notna(a) and re.search(r"\b(0[12]\d{3})\b", str(a)) else None
    )
    df["lat"]  = pd.to_numeric(df.get("lat"),  errors="coerce")
    df["long"] = pd.to_numeric(df.get("long"), errors="coerce")
    counts = df.groupby("zip").size().rename("complaint_count").reset_index()
    return counts, df


# ── MERGED ANALYSIS ───────────────────────────────────────────────────────────

def analyze(df_div, counts):
    merged = df_div.merge(counts, on="zip", how="left")
    merged["complaint_count"] = merged["complaint_count"].fillna(0)

    # Complaints per 1000 residents
    merged["complaints_per_1k"] = (merged["complaint_count"] / merged["total_pop"]) * 1000

    print("\n=== RACIAL DIVERSITY × 311 COMPLAINTS ===")
    print("\nAlpha diversity metrics by ZIP (sorted by Shannon):")
    print(
        merged[["zip","shannon","simpson","richness","evenness","dominance",
                "pct_white","pct_black","pct_asian","complaint_count","complaints_per_1k"]]
        .sort_values("shannon", ascending=False)
        .to_string(index=False)
    )

    print("\n--- Pearson correlations with raw complaint count ---")
    for col in ["shannon","simpson","richness","evenness","dominance",
                "pct_white","pct_black","pct_asian"]:
        r, p = stats.pearsonr(merged[col], merged["complaint_count"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        print(f"  {col:<22} r={r:+.3f}  p={p:.3f} {sig}")

    print("\n--- Pearson correlations with complaints per 1,000 residents ---")
    for col in ["shannon","simpson","richness","evenness","dominance",
                "pct_white","pct_black","pct_asian"]:
        r, p = stats.pearsonr(merged[col], merged["complaints_per_1k"])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        print(f"  {col:<22} r={r:+.3f}  p={p:.3f} {sig}")

    return merged


# ── CHARTS ────────────────────────────────────────────────────────────────────

def build_charts(merged, bc_matrix):
    fig = plt.figure(figsize=(16, 14))
    fig.suptitle("Boston Racial Diversity × 311 Environmental Complaints",
                 fontsize=14, fontweight="bold")

    gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.38)

    cmap_div = plt.cm.plasma

    def scatter_vs_complaints(ax, x_col, x_label, y_col="complaint_count", y_label=None):
        y_label = y_label or ("311 Complaints" if y_col == "complaint_count"
                              else "Complaints per 1k residents")
        c = cmap_div((merged[x_col] - merged[x_col].min()) /
                     (merged[x_col].max() - merged[x_col].min() + 1e-9))
        ax.scatter(merged[x_col], merged[y_col], c=c, s=70, edgecolors="white",
                   linewidths=0.5, zorder=3)
        for _, row in merged.iterrows():
            ax.annotate(row["zip"][-3:], (row[x_col], row[y_col]),
                        fontsize=5.5, alpha=0.65, ha="center", va="bottom")
        z = np.polyfit(merged[x_col], merged[y_col], 1)
        xr = np.linspace(merged[x_col].min(), merged[x_col].max(), 100)
        ax.plot(xr, np.poly1d(z)(xr), "k--", alpha=0.35, linewidth=1)
        r, p = stats.pearsonr(merged[x_col], merged[y_col])
        sig = "**" if p < 0.05 else ("*" if p < 0.10 else "")
        ax.set_title(f"{x_label} vs {y_label}\nr={r:+.2f} p={p:.2f}{sig}", fontsize=9)
        ax.set_xlabel(x_label, fontsize=8)
        ax.set_ylabel(y_label, fontsize=8)

    # Row 1: Shannon, Simpson, Richness vs raw complaints
    scatter_vs_complaints(fig.add_subplot(gs[0,0]), "shannon",  "Shannon Diversity (H)")
    scatter_vs_complaints(fig.add_subplot(gs[0,1]), "simpson",  "Simpson Diversity (1-D)")
    scatter_vs_complaints(fig.add_subplot(gs[0,2]), "evenness", "Pielou Evenness (J)")

    # Row 2: same metrics vs complaints per 1k
    scatter_vs_complaints(fig.add_subplot(gs[1,0]), "shannon",  "Shannon H",
                          y_col="complaints_per_1k", y_label="Per 1k residents")
    scatter_vs_complaints(fig.add_subplot(gs[1,1]), "dominance","Dominance (max group %)",
                          y_col="complaints_per_1k", y_label="Per 1k residents")
    scatter_vs_complaints(fig.add_subplot(gs[1,2]), "pct_white","% White",
                          y_col="complaints_per_1k", y_label="Per 1k residents")

    # Row 3: Beta diversity heatmap + Shannon bar chart
    ax_bc = fig.add_subplot(gs[2, :2])
    zips  = bc_matrix.index.tolist()
    # sort by Shannon for readability
    order = merged.set_index("zip")["shannon"].reindex(zips).sort_values(ascending=False).index.tolist()
    bc_sorted = bc_matrix.loc[order, order]
    im = ax_bc.imshow(bc_sorted.values, cmap="YlOrRd", aspect="auto")
    ax_bc.set_xticks(range(len(order)))
    ax_bc.set_yticks(range(len(order)))
    ax_bc.set_xticklabels([z[-3:] for z in order], fontsize=6, rotation=90)
    ax_bc.set_yticklabels([z[-3:] for z in order], fontsize=6)
    plt.colorbar(im, ax=ax_bc, label="Bray-Curtis dissimilarity")
    ax_bc.set_title("Beta Diversity: Bray-Curtis Dissimilarity\n(ZIPs sorted by Shannon, high→low)", fontsize=9)

    ax_bar = fig.add_subplot(gs[2, 2])
    srt = merged.sort_values("shannon")
    bar_c = cmap_div((srt["shannon"] - srt["shannon"].min()) /
                      (srt["shannon"].max() - srt["shannon"].min()))
    ax_bar.barh(srt["zip"], srt["shannon"], color=bar_c)
    ax_bar.set_xlabel("Shannon H", fontsize=8)
    ax_bar.set_title("Shannon Diversity by ZIP\n(high = more diverse)", fontsize=9)

    plt.savefig("diversity_311_analysis.png", dpi=140, bbox_inches="tight")
    print("Saved: diversity_311_analysis.png")
    plt.close()


# ── STATIC MAP ────────────────────────────────────────────────────────────────

def build_map(merged):
    fig, axes = plt.subplots(1, 2, figsize=(18, 9))
    fig.patch.set_facecolor("#0d1117")

    for ax in axes:
        ax.set_facecolor("#161b22")
        ax.set_xlim(-71.20, -70.98)
        ax.set_ylim(42.22, 42.42)
        ax.tick_params(colors="#888")
        for sp in ax.spines.values():
            sp.set_edgecolor("#333")

    # Left: Shannon diversity map
    ax = axes[0]
    cmap_sh = plt.cm.plasma
    sh_min, sh_max = merged["shannon"].min(), merged["shannon"].max()
    for _, row in merged.iterrows():
        c = ZIP_CENTROIDS.get(row["zip"])
        if not c:
            continue
        norm = (row["shannon"] - sh_min) / (sh_max - sh_min)
        color = cmap_sh(norm)
        ax.scatter(c[1], c[0], s=1600, c=[color], alpha=0.7, zorder=2)
        ax.text(c[1], c[0]+0.003, row["zip"][-3:],
                fontsize=6.5, color="white", ha="center", fontweight="bold", zorder=5)
        ax.text(c[1], c[0]-0.003, f"H={row['shannon']:.2f}",
                fontsize=5.5, color="white", ha="center", zorder=5, alpha=0.85)
    sm1 = plt.cm.ScalarMappable(cmap=cmap_sh,
                                 norm=mcolors.Normalize(vmin=sh_min, vmax=sh_max))
    sm1.set_array([])
    cb1 = plt.colorbar(sm1, ax=ax, fraction=0.03)
    cb1.set_label("Shannon H (higher = more diverse)", color="white")
    cb1.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb1.ax.yaxis.get_ticklabels(), color="white")
    ax.set_title("Racial Diversity (Shannon H) by ZIP", color="white", fontsize=11, pad=8)

    # Right: Complaints per 1k, sized by Shannon
    ax = axes[1]
    cmap_c = plt.cm.RdYlGn_r
    c_min  = merged["complaints_per_1k"].min()
    c_max  = merged["complaints_per_1k"].max()
    max_sh = merged["shannon"].max()
    for _, row in merged.iterrows():
        c = ZIP_CENTROIDS.get(row["zip"])
        if not c:
            continue
        norm  = (row["complaints_per_1k"] - c_min) / (c_max - c_min + 1e-9)
        color = cmap_c(norm)
        size  = 300 + 1800 * (row["shannon"] / max_sh)
        ax.scatter(c[1], c[0], s=size, c=[color], alpha=0.65, zorder=2)
        ax.text(c[1], c[0]+0.003, row["zip"][-3:],
                fontsize=6.5, color="white", ha="center", fontweight="bold", zorder=5)
        ax.text(c[1], c[0]-0.004, f"{row['complaints_per_1k']:.1f}/1k",
                fontsize=5.5, color="white", ha="center", zorder=5, alpha=0.85)
    sm2 = plt.cm.ScalarMappable(cmap=cmap_c,
                                 norm=mcolors.Normalize(vmin=c_min, vmax=c_max))
    sm2.set_array([])
    cb2 = plt.colorbar(sm2, ax=ax, fraction=0.03)
    cb2.set_label("311 complaints per 1,000 residents", color="white")
    cb2.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb2.ax.yaxis.get_ticklabels(), color="white")
    ax.set_title("311 Complaints/1k (color) × Diversity (bubble size)",
                 color="white", fontsize=11, pad=8)
    ax.text(0.02, 0.02, "Larger bubble = more diverse ZIP",
            transform=ax.transAxes, color="#aaa", fontsize=8)

    plt.tight_layout()
    plt.savefig("diversity_311_map.png", dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print("Saved: diversity_311_map.png")
    plt.close()


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("Building race/diversity dataframe...")
    df_div = build_race_df()

    print("Computing beta diversity (Bray-Curtis)...")
    bc_matrix = beta_diversity(df_div)

    print("Loading 311 data...")
    counts, df_311_raw = load_311()

    merged = analyze(df_div, counts)
    merged.to_csv("diversity_zip_table.csv", index=False)
    print("Saved: diversity_zip_table.csv")

    print("\nBuilding charts...")
    build_charts(merged, bc_matrix)

    print("Building map...")
    build_map(merged)


if __name__ == "__main__":
    main()
