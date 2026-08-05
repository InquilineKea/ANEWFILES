#!/usr/bin/env python3
"""
Figures for the spatial / compositional entropy analysis.

    python3 entropy_analysis.py
    python3 make_entropy_figures.py

fig6: the spatial-entropy instrument, validated on real MERFISH coordinates.
fig7: compositional entropy vs age -- which runs OPPOSITE to the hypothesis.

Colour jobs, per the form heuristic:
* fig6 left is a magnitude map -> SEQUENTIAL, one hue light->dark. Never a
  rainbow: a rainbow ramp invents category boundaries where the data is
  continuous, and reverses under common colour-vision deficiencies.
* fig6 right and fig7 encode adult vs aged / rank -> emphasis, grey + one accent.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "abc-entropy-figures"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entropy_analysis import shannon  # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"
GREY, GREY_L = "#9a9a94", "#d5d5cf"
INK, INK_2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.size": 10,
    "axes.edgecolor": "#c9c9c3", "axes.labelcolor": INK_2,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans",
})

# Single-hue sequential ramp anchored on the validated blue: light = ordered,
# dark = mixed. Monotonic in lightness, which is what makes it readable in
# greyscale and under CVD.
SEQ = LinearSegmentedColormap.from_list(
    "entropy_seq", ["#eef3fa", "#c3d8f0", "#8fb7e4", "#5892d8", "#2a78d6", "#17457d"])


def section_entropy(csv: str, section: str, k: int = 20) -> pd.DataFrame:
    cols = ["brain_section_label", "subclass", "parcellation_division",
            "x_reconstructed", "y_reconstructed"]
    it = pd.read_csv(csv, usecols=cols, chunksize=500_000)
    parts = [c[c.brain_section_label == section] for c in it]
    g = pd.concat(parts, ignore_index=True)
    xy = g[["x_reconstructed", "y_reconstructed"]].to_numpy(float)
    codes = g.subclass.astype(str).to_numpy()
    _, idx = cKDTree(xy).query(xy, k=k + 1)
    neigh = codes[idx[:, 1:]]
    g["entropy"] = [shannon(np.unique(r, return_counts=True)[1]) for r in neigh]
    return g


def fig6(outdir: str, merfish_csv: str) -> None:
    S = pd.read_csv(os.path.join(outdir, "results_entropy_spatial_merfish.csv.gz"))
    sec = S.section.value_counts().idxmax()
    g = section_entropy(merfish_csv, sec)

    fig = plt.figure(figsize=(13.2, 5.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.32, 1.0], wspace=0.42)

    # -- left: the map --
    ax = fig.add_subplot(gs[0, 0])
    order = np.argsort(g.entropy.to_numpy())  # draw ordered cells last
    sc = ax.scatter(g.x_reconstructed.to_numpy()[order],
                    g.y_reconstructed.to_numpy()[order],
                    c=g.entropy.to_numpy()[order], s=1.1, cmap=SEQ,
                    vmin=np.nanpercentile(g.entropy, 2),
                    vmax=np.nanpercentile(g.entropy, 98), linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    # Horizontal colorbar beneath the map: a vertical one lands between the
    # panels and collides with the right panel's tick labels.
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.042,
                      pad=0.03, aspect=38)
    cb.set_label("neighbourhood entropy (nats) — pale = locally uniform, dark = locally mixed",
                 fontsize=9)
    cb.outline.set_visible(False)
    ax.set_title("Local cell-type entropy, one coronal section",
                 fontsize=12.5, color=INK, fontweight="bold", loc="left", pad=16)
    ax.text(0, 1.015, "%s · k=20 neighbours" % sec, transform=ax.transAxes,
            fontsize=9.2, color=INK_2, va="bottom")

    # -- right: HPF ranking --
    ax2 = fig.add_subplot(gs[0, 1])
    hpf = S[S.division == "HPF"]
    r = (hpf.groupby("structure").entropy.agg(["median", "count"])
         .query("count >= 500").sort_values("median"))
    y = np.arange(len(r))
    colors = [BLUE if s in ("DG", "CA3", "CA1", "CA2") else GREY for s in r.index]
    ax2.barh(y, r["median"], color=colors, height=0.66, zorder=3)
    for i, (name, row) in enumerate(r.iterrows()):
        ax2.text(row["median"] + 0.02, i, "%.2f" % row["median"], va="center",
                 fontsize=8.4, color=INK_2)
    ax2.set_yticks(y)
    ax2.set_yticklabels(r.index, fontsize=9.5)
    # Direct-label by colouring the tick text: the region names already
    # identify each bar, so a legend box would be redundant and it collides
    # with the shortest bar's value label.
    for tick, name in zip(ax2.get_yticklabels(), r.index):
        if name in ("DG", "CA3", "CA1", "CA2"):
            tick.set_color(BLUE)
            tick.set_fontweight("bold")
    ax2.set_xlim(0, 1.85)
    ax2.grid(axis="x", color=GREY_L, lw=0.6, alpha=0.7)
    ax2.set_axisbelow(True)
    ax2.set_xlabel("median neighbourhood entropy (nats)", fontsize=9.5)
    ax2.set_title("Hippocampal subregions, ordered", fontsize=12.5, color=INK,
                  fontweight="bold", loc="left", pad=16)
    ax2.text(0, 1.015, "Dentate gyrus is the most ordered structure — the measure "
             "recovers known anatomy. Blue = dentate / CA fields.", transform=ax2.transAxes, fontsize=9.2,
             color=INK_2, va="bottom")

    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig6_spatial_entropy.%s" % ext), dpi=200,
                    bbox_inches="tight",
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig6_spatial_entropy.png")


def fig7(outdir: str) -> None:
    E = pd.read_csv(os.path.join(outdir, "results_entropy_aging.csv"))
    R = pd.read_csv(os.path.join(outdir, "results_entropy_aging_test.csv"))
    scopes = ["all cells", "non-neuronal only"]
    titles = {"all cells": "All cells",
              "non-neuronal only": "Non-neuronal only"}
    subs = {"all cells": "confounded by neuronal dropout",
            "non-neuronal only": "confound removed"}

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8), sharey=True)
    rng = np.random.default_rng(0)
    for ax, scope in zip(axes, scopes):
        sub = E[E.scope == scope]
        row = R[R.scope == scope].iloc[0]
        for j, (age, c) in enumerate([("adult", GREY), ("aged", BLUE)]):
            v = sub[sub.age == age].entropy.to_numpy()
            x = j + rng.uniform(-0.15, 0.15, len(v))
            ax.scatter(x, v, s=36, color=c, alpha=0.78, zorder=3,
                       edgecolor=SURFACE, linewidth=0.8)
            ax.plot([j - 0.27, j + 0.27], [np.median(v)] * 2, color=INK, lw=2.6,
                    zorder=5, solid_capstyle="round")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Adult", "Aged"], fontsize=10)
        ax.set_xlim(-0.5, 1.5)
        ax.grid(axis="y", color=GREY_L, lw=0.6, alpha=0.7)
        ax.set_axisbelow(True)
        ax.set_title(titles[scope], fontsize=12, color=INK, fontweight="bold",
                     loc="left", pad=18)
        ax.text(0, 1.012, "%s · shift %+.3f nats, q=%.3f" % (subs[scope], row.shift_nats, row.q),
                transform=ax.transAxes, fontsize=9, color=INK_2, va="bottom")
        ax.annotate("", xy=(1.34, np.median(sub[sub.age == "aged"].entropy)),
                    xytext=(1.34, np.median(sub[sub.age == "adult"].entropy)),
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.0))
    axes[0].set_ylabel("cell-type entropy, rarefied to 300 cells (nats)", fontsize=9.4)
    fig.suptitle("Cell-type entropy DECREASES with age — opposite to the hypothesis",
                 fontsize=13, fontweight="bold", color=INK, x=0.008, ha="left", y=0.995)
    fig.text(0.008, 0.005, "Compositional entropy, NOT spatial entropy — this dataset has no "
             "coordinates. Each point is one donor×region sample.",
             fontsize=8.4, color=INK_2)
    fig.tight_layout(rect=[0, 0.04, 1, 0.93])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig7_entropy_age.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig7_entropy_age.png")


def main() -> int:
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    merfish = sys.argv[2] if len(sys.argv) > 2 else "merfish_cells.csv"
    print("generating entropy figures ...")
    fig7(outdir)
    if os.path.exists(merfish):
        fig6(outdir, merfish)
    else:
        print("  fig6 skipped: %s not found" % merfish)
    return 0


if __name__ == "__main__":
    sys.exit(main())
