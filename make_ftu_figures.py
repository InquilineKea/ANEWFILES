#!/usr/bin/env python3
"""
Figure 8: FTU innervation — what soma-based spatial data can and cannot show.

    python3 ftu_innervation.py merfish_cells.csv
    python3 make_ftu_figures.py

Left panel: intrinsic excitatory/inhibitory balance per hippocampal FTU.
  Job is magnitude with one item as the story (DG) -> EMPHASIS, accent + grey.
Right panel: distance from each FTU to the nearest modulatory soma.
  Job is magnitude on a grid -> HEATMAP with a SEQUENTIAL single-hue ramp.
  Five modulators as five categorical colours would breach the all-pairs cap
  and encode magnitude as identity, which is the wrong job.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "abc-ftu-figures"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

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

SEQ = LinearSegmentedColormap.from_list(
    "dist_seq", ["#eef3fa", "#c3d8f0", "#8fb7e4", "#5892d8", "#2a78d6", "#17457d"])

ORDER = ["DG", "CA3", "CA2", "CA1", "SUB", "ProS", "POST", "PRE", "PAR", "ENTl", "ENTm"]
FULL = {"DG": "DG (dentate)", "CA3": "CA3", "CA2": "CA2", "CA1": "CA1",
        "SUB": "Subiculum", "ProS": "Prosubiculum", "POST": "Postsubiculum",
        "PRE": "Presubiculum", "PAR": "Parasubiculum",
        "ENTl": "Entorhinal (lat.)", "ENTm": "Entorhinal (med.)"}


def main() -> int:
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    I = pd.read_csv(os.path.join(outdir, "results_ftu_intrinsic.csv")).set_index("FTU")
    M = pd.read_csv(os.path.join(outdir, "results_ftu_extrinsic.csv")).set_index("FTU")
    order = [f for f in ORDER if f in I.index]

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6),
                             gridspec_kw={"width_ratios": [1.0, 1.12], "wspace": 0.30})

    # ---- left: E/I balance ----
    ax = axes[0]
    y = np.arange(len(order))
    vals = I.loc[order, "EI_ratio"].to_numpy()
    colors = [ORANGE if f == "DG" else GREY for f in order]
    ax.barh(y, vals, color=colors, height=0.66, zorder=3)
    for i, (f, v) in enumerate(zip(order, vals)):
        ax.text(v + 0.25, i, "%.1f" % v, va="center", fontsize=8.8,
                color=ORANGE if f == "DG" else INK_2,
                fontweight="bold" if f == "DG" else "normal")
    ax.set_yticks(y)
    ax.set_yticklabels([FULL[f] for f in order], fontsize=9.5)
    ax.get_yticklabels()[order.index("DG")].set_color(ORANGE)
    ax.get_yticklabels()[order.index("DG")].set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 20)
    ax.grid(axis="x", color=GREY_L, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel("excitatory : inhibitory neuron ratio", fontsize=9.5)
    ax.set_title("Intrinsic circuit composition", fontsize=12.5, color=INK,
                 fontweight="bold", loc="left", pad=17)
    ax.text(0, 1.018, "Measurable from cell bodies · dentate gyrus is far the most "
            "excitation-dominated", transform=ax.transAxes, fontsize=9.2,
            color=INK_2, va="bottom")

    # ---- right: distance to nearest modulatory soma ----
    ax2 = axes[1]
    mods = [c.split("_")[1] for c in M.columns if c.startswith("d_") and c.endswith("median_mm")]
    D = np.array([[M.loc[f, "d_%s_median_mm" % m] for m in mods] for f in order])
    im = ax2.imshow(D, cmap=SEQ, aspect="auto", vmin=0, vmax=np.nanmax(D))
    for i in range(D.shape[0]):
        for j in range(D.shape[1]):
            v = D[i, j]
            ax2.text(j, i, "%.1f" % v, ha="center", va="center", fontsize=8.6,
                     color="white" if v > 0.62 * np.nanmax(D) else INK)
    ax2.set_xticks(range(len(mods)))
    ax2.set_xticklabels(["Cholin.", "Dopa.", "Seroton.", "Noradr.", "Histam."][:len(mods)],
                        fontsize=9.5)
    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels([FULL[f] for f in order], fontsize=9.5)
    for s in ax2.spines.values():
        s.set_visible(False)
    ax2.set_xticks(np.arange(-0.5, len(mods), 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax2.grid(which="minor", color=SURFACE, linewidth=2)
    ax2.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax2, fraction=0.036, pad=0.02)
    cb.set_label("distance (mm)", fontsize=9)
    cb.outline.set_visible(False)
    ax2.set_title("Distance to nearest modulatory cell body", fontsize=12.5,
                  color=INK, fontweight="bold", loc="left", pad=17)
    ax2.text(0, 1.018, "Zero modulatory somata inside ANY hippocampal FTU — "
             "all modulation is extrinsic", transform=ax2.transAxes,
             fontsize=9.2, color=INK_2, va="bottom")

    fig.text(0.008, 0.005,
             "Soma-based data cannot measure innervation: axons are not segmented. "
             "These distances prove the input is extrinsic; they do not quantify how densely it arrives.",
             fontsize=8.4, color=INK_2)
    fig.tight_layout(rect=[0, 0.035, 1, 1])
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, "fig8_ftu_innervation.%s" % ext), dpi=200,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)
    print("  fig8_ftu_innervation.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
