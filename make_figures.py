#!/usr/bin/env python3
"""
Figures for the aging brain composition analysis (Zeng-Aging-Mouse-10Xv3).

    python3 aging_composition.py aging_cells.csv     # writes the result CSVs
    python3 make_figures.py

Design notes
------------
* Colour does one job per figure. Figures 1, 2 and 4 use EMPHASIS (one accent
  hue plus a de-emphasis grey) rather than a categorical palette, because the
  story is "these went up and this one did not" -- a categorical scheme would
  spread attention evenly across series and bury exactly that contrast.
* Adult -> aged is a before/after per item, so figures 2 and 3 are dumbbells.
* Figure 3 uses a log x-axis because the immune subtypes span 0.25% to 96%,
  a ~400x range that a linear axis renders as one bar and three invisible ones.
* Palette validated with the dataviz validator (light surface, categorical):
  all checks pass; the aqua slot carries a contrast warning and is therefore
  used only with visible direct labels.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aging_composition import FOCUS, clr, donor_region_props, load  # noqa: E402

# --- palette (validated) ---
BLUE = "#2a78d6"      # accent / aged
ORANGE = "#eb6834"    # the highlighted null result
GREY = "#9a9a94"      # de-emphasis / adult
GREY_L = "#d5d5cf"    # connectors, grid
INK = "#0b0b0b"
INK_2 = "#52514e"
SURFACE = "#fcfcfb"

# Deterministic SVG output: without a fixed hashsalt matplotlib mints random
# element IDs each run, so re-running produces a large but meaningless diff.
matplotlib.rcParams["svg.hashsalt"] = "abc-aging-figures"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 10,
    "axes.edgecolor": "#c9c9c3", "axes.labelcolor": INK_2,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans",
})

REGION_ORDER = ["Isocortex", "HPF", "STR", "PAL", "HY", "MB", "P"]
REGION_FULL = {"Isocortex": "Isocortex", "HPF": "Hippocampal form.", "STR": "Striatum",
               "PAL": "Pallidum", "HY": "Hypothalamus", "MB": "Midbrain", "P": "Pons"}


def _tidy(ax, xlabel: str = "", title: str = "", sub: str = "") -> None:
    ax.grid(axis="x", color=GREY_L, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9.5)
    if title:
        ax.set_title(title, fontsize=12.5, color=INK, fontweight="bold",
                     loc="left", pad=16 if sub else 8)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=9.2,
                color=INK_2, va="bottom")


def bootstrap_ci(a: np.ndarray, b: np.ndarray, n: int = 4000,
                 seed: int = 0) -> tuple[float, float]:
    """Percentile CI for the median difference (aged - adult)."""
    rng = np.random.default_rng(seed)
    d = [np.median(rng.choice(b, len(b))) - np.median(rng.choice(a, len(a)))
         for _ in range(n)]
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


# ---------------------------------------------------------------- figure 1
def fig1_effects(props, meta, outdir):
    """Which cell types shift with age: region-stratified CLR effect + 95% CI."""
    clr_all = pd.DataFrame(clr(props.to_numpy()), columns=props.columns)
    pooled = pd.read_csv(os.path.join(outdir, "results_aging_pooled.csv"))

    rows = []
    for _, r in pooled.iterrows():
        t = r.subclass
        adj_a, adj_b = [], []
        for region in meta.region.unique():
            m = (meta.region == region).to_numpy()
            ad, ag = m & (meta.age == "adult").to_numpy(), m & (meta.age == "aged").to_numpy()
            if ad.sum() < 3 or ag.sum() < 3:
                continue
            base = np.median(clr_all.loc[ad, t])
            adj_a.append(clr_all.loc[ad, t].to_numpy() - base)
            adj_b.append(clr_all.loc[ag, t].to_numpy() - base)
        a, b = np.concatenate(adj_a), np.concatenate(adj_b)
        lo, hi = bootstrap_ci(a, b)
        rows.append({"label": r.label, "effect": r.clr_shift, "lo": lo, "hi": hi,
                     "q": r.q, "sig": r.q < 0.05})
    d = pd.DataFrame(rows).sort_values("effect")

    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    y = np.arange(len(d))
    ax.axvline(0, color=INK_2, lw=1.1, zorder=1)
    for i, (_, r) in enumerate(d.iterrows()):
        is_endo = r.label == "Endothelial"
        c = ORANGE if is_endo else (BLUE if r.sig else GREY)
        ax.plot([r.lo, r.hi], [i, i], color=c, lw=2.0, alpha=0.55,
                solid_capstyle="round", zorder=2)
        ax.scatter([r["effect"]], [i], s=110, color=c, zorder=3,
                   edgecolor=SURFACE, linewidth=2)
        lab = "q=%.3f" % r.q if r.q >= 0.001 else "q<0.001"
        ax.text(r.hi + 0.045, i, lab, va="center", fontsize=8.6,
                color=INK if (r.sig or is_endo) else INK_2)
    ax.set_yticks(y)
    ax.set_yticklabels(d.label, fontsize=10)
    for tick, (_, r) in zip(ax.get_yticklabels(), d.iterrows()):
        if r.label == "Endothelial":
            tick.set_color(ORANGE)
            tick.set_fontweight("bold")
        elif r.sig:
            tick.set_color(INK)
    ax.set_xlim(-0.85, 1.35)
    _tidy(ax, "CLR shift, aged − adult  (→ higher share in aged brain)",
          "Vascular mural and immune cells rise with age; endothelium does not",
          "58 donors, region-stratified · dot = median, bar = 95% bootstrap CI")
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=BLUE, ms=9, label="Significant (q<0.05)"),
        Line2D([], [], marker="o", ls="", color=GREY, ms=9, label="Not significant"),
        Line2D([], [], marker="o", ls="", color=ORANGE, ms=9, label="Endothelial (the null)"),
    ], loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig1_age_effects.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig1_age_effects.png")


# ---------------------------------------------------------------- figure 2
def fig2_regions(props, meta, outdir):
    """Microglia vs endothelium per region: the contrast, as dumbbells."""
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), sharey=True)
    # Colour encodes AGE here, so it must mean the same thing in both panels.
    # The verdict is carried by the subtitle, never by recolouring the variable.
    panels = [("334 Microglia NN", "Microglia", "Rises in 6 of 7 regions"),
              ("333 Endo NN", "Endothelial", "Flat — no consistent direction")]
    regions = [r for r in REGION_ORDER if r in set(meta.region)]

    for ax, (col, name, verdict) in zip(axes, panels):
        for i, reg in enumerate(regions):
            m = (meta.region == reg).to_numpy()
            a = 100 * props.loc[m & (meta.age == "adult").to_numpy(), col].median()
            b = 100 * props.loc[m & (meta.age == "aged").to_numpy(), col].median()
            ax.plot([a, b], [i, i], color=GREY_L, lw=3.5, solid_capstyle="round", zorder=1)
            ax.scatter([a], [i], s=95, color=GREY, zorder=3,
                       edgecolor=SURFACE, linewidth=1.8)
            ax.scatter([b], [i], s=95, color=BLUE, zorder=4,
                       edgecolor=SURFACE, linewidth=1.8)
            ax.text(max(a, b) + 0.7, i, "%.1f→%.1f" % (a, b), va="center",
                    fontsize=8.6, color=INK_2)
        ax.set_yticks(range(len(regions)))
        ax.set_yticklabels([REGION_FULL[r] for r in regions], fontsize=10)
        ax.set_xlim(0, 26)
        ax.set_ylim(len(regions) - 0.4, -0.6)   # explicit: top-to-bottom order
        _tidy(ax, "% of all cells in that region", name, verdict)

    axes[1].legend(handles=[
        Line2D([], [], marker="o", ls="", color=GREY, ms=9, label="Adult (~P56)"),
        Line2D([], [], marker="o", ls="", color=BLUE, ms=9, label="Aged (~P550)"),
    ], loc="lower right", frameon=True, framealpha=0.95, edgecolor="#dedbd2",
        fontsize=9)
    fig.suptitle("Age effect is cell-type specific, not a global vascular or glial shift",
                 fontsize=13, fontweight="bold", color=INK, x=0.008, ha="left", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig2_region_dumbbell.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig2_region_dumbbell.png")


# ---------------------------------------------------------------- figure 3
def fig3_immune(outdir):
    """Immune compartment make-up, adult vs aged. Log axis: 0.25%..96%."""
    ip = pd.read_csv(os.path.join(outdir, "aging_immune_props.csv"))
    cols = [c for c in ip.columns if c.startswith(("33", "334"))]
    order = ["334 Microglia NN", "335 BAM NN", "338 Lymphoid NN", "337 DC NN"]
    order = [c for c in order if c in ip.columns]
    sparse = {"338 Lymphoid NN", "337 DC NN"}

    fig, ax = plt.subplots(figsize=(10.4, 3.9))
    for i, c in enumerate(order):
        a = 100 * ip.loc[ip.donor_age_category == "adult", c].median()
        b = 100 * ip.loc[ip.donor_age_category == "aged", c].median()
        ax.plot([a, b], [i, i], color=GREY_L, lw=3.5, solid_capstyle="round", zorder=1)
        ax.scatter([a], [i], s=95, color=GREY, zorder=3, edgecolor=SURFACE, linewidth=1.8)
        ax.scatter([b], [i], s=95, color=BLUE, zorder=4, edgecolor=SURFACE, linewidth=1.8)
        # Labels sit left of the marks for the widest row so nothing overflows.
        txt = "%.2f%% → %.2f%%" % (a, b)
        if c in sparse:
            txt += "  ⚠ sparse"
        if max(a, b) > 20:
            ax.text(min(a, b) * 0.72, i, txt, va="center", ha="right",
                    fontsize=9, color=INK_2)
        else:
            ax.text(max(a, b) * 1.35, i, txt, va="center", fontsize=9, color=INK_2)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([FOCUS.get(c, c) for c in order], fontsize=10)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 260)
    ax.set_xticks([0.25, 1, 5, 25, 100])
    ax.set_xticklabels(["0.25%", "1%", "5%", "25%", "100%"])
    ax.set_ylim(len(order) - 0.4, -0.6)
    _tidy(ax, "share of the immune compartment (log scale)",
          "Within immune cells, microglia crowd out border macrophages",
          "Medians across 58 donors · ⚠ = absent from >20% of samples, provisional")
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=GREY, ms=9, label="Adult"),
        Line2D([], [], marker="o", ls="", color=BLUE, ms=9, label="Aged"),
    ], loc="center right", frameon=True, framealpha=0.95, edgecolor="#dedbd2",
        fontsize=9)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig3_immune_composition.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig3_immune_composition.png")


# ---------------------------------------------------------------- figure 4
def fig4_donors(props, meta, outdir):
    """Per-donor spread: the actual evidence behind figure 1."""
    clr_all = pd.DataFrame(clr(props.to_numpy()), columns=props.columns)
    targets = [("334 Microglia NN", "Microglia", "q=0.004"),
               ("331 Peri NN", "Pericyte", "q=0.001"),
               ("333 Endo NN", "Endothelial", "q=0.83, n.s.")]
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.6), sharey=True)
    rng = np.random.default_rng(0)

    # Shared, robust y-limits. A couple of near-zero observations sit around
    # CLR -4 and would otherwise squash every distribution into a third of the
    # panel; they are clipped to the axis edge and counted in a footnote
    # rather than silently dropped.
    series = {}
    for col, _, _ in targets:
        vals = {}
        for region in meta.region.unique():
            m = (meta.region == region).to_numpy()
            ad, ag = m & (meta.age == "adult").to_numpy(), m & (meta.age == "aged").to_numpy()
            if ad.sum() < 3 or ag.sum() < 3:
                continue
            base = np.median(clr_all.loc[ad, col])
            vals.setdefault("adult", []).append(clr_all.loc[ad, col].to_numpy() - base)
            vals.setdefault("aged", []).append(clr_all.loc[ag, col].to_numpy() - base)
        series[col] = (np.concatenate(vals["adult"]), np.concatenate(vals["aged"]))
    allv = np.concatenate([np.concatenate(v) for v in series.values()])
    lo, hi = np.percentile(allv, 2.0), np.percentile(allv, 98.0)
    pad = 0.15 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    n_off = int((allv < lo).sum() + (allv > hi).sum())

    for ax, (col, name, qtxt) in zip(axes, targets):
        a, b = series[col]
        for j, (grp, c) in enumerate([(a, GREY), (b, BLUE)]):
            x = j + rng.uniform(-0.14, 0.14, len(grp))
            inside = (grp >= lo) & (grp <= hi)
            ax.scatter(x[inside], grp[inside], s=34, color=c, alpha=0.75, zorder=3,
                       edgecolor=SURFACE, linewidth=0.8)
            # Off-scale points get a caret at the edge, so they are visibly
            # "beyond the axis" rather than masquerading as boundary values.
            for sign, edge in ((-1, lo), (1, hi)):
                off = ~inside & (np.sign(grp - hi if sign > 0 else grp - lo) == sign)
                if off.any():
                    ax.scatter(x[off], np.full(off.sum(), edge), s=44, color=c,
                               marker="v" if sign < 0 else "^", zorder=3,
                               edgecolor=SURFACE, linewidth=0.8)
            ax.plot([j - 0.27, j + 0.27], [np.median(grp)] * 2, color=INK,
                    lw=2.6, zorder=5, solid_capstyle="round")
        ax.axhline(0, color=GREY_L, lw=1.2, zorder=1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Adult\n(n=%d)" % len(a), "Aged\n(n=%d)" % len(b)], fontsize=9.5)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(lo, hi)
        ax.grid(axis="y", color=GREY_L, lw=0.6, alpha=0.7)
        ax.set_axisbelow(True)
        ax.set_title(name, fontsize=11.5, color=INK, fontweight="bold", loc="left")
        ax.text(1.0, 1.015, qtxt, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=9, color=INK if "n.s." not in qtxt else ORANGE)
        if ax is axes[0]:
            ax.set_ylabel("CLR, centred on each region's adult median", fontsize=9.2)
    fig.suptitle("Every donor×region observation — thick bar = median",
                 fontsize=12.5, fontweight="bold", color=INK, x=0.008, ha="left", y=0.995)
    fig.text(0.008, 0.005, "Grey = adult, blue = aged. %d value(s) fall outside the axis and are shown "
             "as carets at the edge." % n_off,
             fontsize=8.4, color=INK_2)
    fig.tight_layout(rect=[0.0, 0.035, 1, 0.94])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig4_donor_spread.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig4_donor_spread.png")



# ---------------------------------------------------------------- figure 5
def fig5_confound(outdir):
    """What survives the neuronal-denominator confound.

    Two analyses of the same donors: share of ALL cells (confounded by
    differential neuronal recovery) versus share of NON-NEURONAL cells only
    (which conditions that denominator away). Two series -> categorical
    palette slots 1 and 2, both validated.
    """
    m = pd.read_csv(os.path.join(outdir, "results_aging_confound.csv"))
    m = m.sort_values("clr_nn", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.6, 5.8))
    ax.axvline(0, color=INK_2, lw=1.1, zorder=2)
    for i, r in m.iterrows():
        ax.plot([r.clr_whole, r.clr_nn], [i, i], color=GREY_L, lw=3.2,
                solid_capstyle="round", zorder=1)
        ax.scatter([r.clr_whole], [i], s=100, color=BLUE, zorder=3,
                   edgecolor=SURFACE, linewidth=1.8)
        ax.scatter([r.clr_nn], [i], s=100, color=ORANGE, zorder=4,
                   edgecolor=SURFACE, linewidth=1.8)
        note = ""
        if r.q_whole < 0.05 and r.q_nn >= 0.05:
            note = "was sig. — did not survive"
        elif r.q_whole >= 0.05 and r.q_nn < 0.05:
            note = "emerges once confound removed"
        if note:
            ax.text(max(r.clr_whole, r.clr_nn) + 0.05, i, note, va="center",
                    fontsize=8.4,
                    color=ORANGE if r.q_nn < 0.05 else INK_2)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels(m.label_nn, fontsize=10)
    for tick, (_, r) in zip(ax.get_yticklabels(), m.iterrows()):
        if r.q_nn < 0.05:
            tick.set_color(ORANGE)
            tick.set_fontweight("bold")
    ax.set_xlim(-0.95, 1.5)
    ax.set_ylim(len(m) - 0.4, -0.6)
    _tidy(ax, "CLR shift, aged − adult",
          "Most age effects were the neuronal denominator moving",
          "Same 58 donors, two denominators · only endothelium and border macrophages survive")
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=BLUE, ms=9,
               label="Share of ALL cells (confounded)"),
        Line2D([], [], marker="o", ls="", color=ORANGE, ms=9,
               label="Share of NON-NEURONAL cells only"),
    ], loc="lower right", frameon=True, framealpha=0.95, edgecolor="#dedbd2",
        fontsize=9)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig5_confound.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig5_confound.png")


def main() -> int:
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    df = load("aging_cells.csv")
    props, meta = donor_region_props(df, 150)
    print("generating figures ...")
    fig1_effects(props, meta, outdir)
    fig2_regions(props, meta, outdir)
    fig3_immune(outdir)
    fig4_donors(props, meta, outdir)
    fig5_confound(outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
