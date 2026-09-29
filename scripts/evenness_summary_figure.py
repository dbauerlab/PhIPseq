"""One figure summarizing library evenness per the plain-language stats
already discussed: for each sample, the P10-P90 "typical range" of
per-peptide read counts, the median, and any near-dropout outliers (<=10
reads) called out explicitly rather than buried in a max/min ratio.
Coloured by treatment and drawn once for all samples (results/figures/)
and once per IP condition (results/figures/<ip>/, see sample_groups.py).

Usage: evenness_summary_figure.py [samples_tsv] [samples_meta_csv]
  (default to reference/samples.tsv and reference/samples_meta.csv)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from sample_groups import figure_groups, legend_outside, load_samples, sample_label, treatment_legend_handles

PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())
RUN_DIR = Path(_paths["run_dir"])
RESULTS_DIR = RUN_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

# Figure thresholds live in reference/figure_params.json (shared with
# overview_figure.py) rather than being duplicated as a literal here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
OUTLIER_THRESHOLD = _fig_params["outlier_threshold"]  # reads; peptides at or below this are called out individually
XTICK_ROTATION = _fig_params["xtick_rotation"]

# Near-dropout markers are black so they stay distinct from the treatment palette.
OUTLIER_COLOR = "black"

# Room below the axes for the rotated "<ip> <tc_id> (N ND, N miss)" tick labels.
TICK_LABEL_HEIGHT_IN = 3.2


def plot_evenness(mat, lib_of_ref, group, colors):
    samples = group["samples"]
    fig, ax = plt.subplots(figsize=(1.4 * len(samples) + 2, 5.5 + TICK_LABEL_HEIGHT_IN))
    xtick_labels = []

    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]
        detected = counts[counts > 0]
        n_missing = int((counts == 0).sum())
        outliers = detected[detected <= OUTLIER_THRESHOLD]
        if len(detected):
            p10, p50, p90 = np.percentile(detected, [10, 50, 90])

            # the tight "typical" band most peptides fall in
            ax.plot([i, i], [p10, p90], color=colors[srow["treatment"]], linewidth=6, solid_capstyle="butt",
                    alpha=0.85)
            ax.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=8, zorder=5)

            # near-dropout outliers, called out individually rather than folded into a min/max ratio
            ax.plot(np.full(len(outliers), i), outliers.values, "v", color=OUTLIER_COLOR, markersize=9,
                    linestyle="", zorder=6)
        label_bits = []
        if len(outliers):
            label_bits.append(f"{len(outliers)} ND")
        if n_missing:
            label_bits.append(f"{n_missing} miss")
        extra = " (" + ", ".join(label_bits) + ")" if label_bits else " (clean)"
        # Kept as a single line: rotated multi-line tick labels lay out very
        # differently from rotated single-line ones (each line becomes its
        # own column) and need much more room than is practical here.
        xtick_labels.append(f"{sample_label(srow, group)}{extra}")

    ax.set_yscale("log")
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(xtick_labels, rotation=XTICK_ROTATION, ha="center", fontsize=8)
    ax.set_ylabel("Reads per peptide (log scale)")
    ax.set_title("Library evenness: how tightly clustered is representation?" + group["title_suffix"]
                 + f"\n(ND = near-dropout, <={OUTLIER_THRESHOLD} reads)")
    handles = [Line2D([], [], color="gray", linewidth=6, label="typical range (10th-90th percentile)"),
               Line2D([], [], marker="o", linestyle="", color="white", markeredgecolor="black", markersize=8,
                      label="median"),
               Line2D([], [], marker="v", linestyle="", color=OUTLIER_COLOR, markersize=9,
                      label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)")]
    legend_outside(ax, handles + treatment_legend_handles(samples, colors), fontsize=8)
    # Fixed widths/heights (in inches, converted to figure fractions) are
    # reserved for the outside legend and the rotated tick labels. Not using
    # bbox_inches="tight" here: with this matplotlib version it can blow up
    # the canvas to a pathological size when combined with a log-scaled axis
    # and many plotted markers.
    fig.subplots_adjust(right=1 - 3.2 / fig.get_figwidth(), bottom=TICK_LABEL_HEIGHT_IN / fig.get_figheight())
    fig.savefig(group["fig_dir"] / "evenness_summary.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {group['fig_dir'] / 'evenness_summary.png'}", file=sys.stderr)


def main():
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    for group in figure_groups(samples, FIG_DIR):
        plot_evenness(mat, lib_of_ref, group, colors)


if __name__ == "__main__":
    main()
