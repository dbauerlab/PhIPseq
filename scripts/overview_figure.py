"""One combined figure comparing samples together: coverage (%% of
designed peptides detected at >=1 read and at >= coverage_min_count reads,
see count_filters.py) and evenness (typical range + near-dropout
outliers) side by side, coloured by treatment. Drawn once for all samples
(results/figures/) and once per IP condition (results/figures/<ip>/, see
sample_groups.py).

Usage: overview_figure.py [samples_tsv] [samples_meta_csv]
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
from matplotlib.patches import Patch

from count_filters import COVERAGE_MIN_COUNT
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

# Figure thresholds/styling live in reference/figure_params.json (shared
# with evenness_summary_figure.py) rather than being duplicated as literals
# here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
OUTLIER_THRESHOLD = _fig_params["outlier_threshold"]  # reads
XTICK_ROTATION = _fig_params["xtick_rotation"]

# Near-dropout markers are black so they stay distinct from the treatment palette.
OUTLIER_COLOR = "black"

# Room below the evenness panel for the rotated "<ip> <tc_id> (N ND, N miss)" tick labels.
TICK_LABEL_HEIGHT_IN = 3.2


def plot_overview(mat, lib_of_ref, group, colors):
    samples = group["samples"]
    fig, (ax_cov, ax_even) = plt.subplots(2, 1, figsize=(0.7 * len(samples) + 5, 8.5 + TICK_LABEL_HEIGHT_IN))

    cov_labels, even_labels = [], []
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        color = colors[srow["treatment"]]
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]
        total = len(own_refs)
        detected = counts[counts > 0]
        n_min = int((counts >= COVERAGE_MIN_COUNT).sum())
        pct = 100 * len(detected) / total if total else 0
        pct_min = 100 * n_min / total if total else 0

        # --- coverage panel: light bar >=1 read, solid bar >= COVERAGE_MIN_COUNT reads ---
        ax_cov.bar(i, pct, color=color, alpha=0.35, edgecolor="black", linewidth=0.5)
        ax_cov.bar(i, pct_min, color=color, edgecolor="black", linewidth=0.5)
        ax_cov.text(i, pct + 1.5, f"{n_min}", ha="center", va="bottom", rotation=90, fontsize=7)

        # --- evenness panel ---
        outliers = detected[detected <= OUTLIER_THRESHOLD]
        if len(detected):
            p10, p50, p90 = np.percentile(detected, [10, 50, 90])
            ax_even.plot([i, i], [p10, p90], color=color, linewidth=6, solid_capstyle="butt", alpha=0.85)
            ax_even.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=7, zorder=5)
            ax_even.plot(np.full(len(outliers), i), outliers.values, "v", color=OUTLIER_COLOR, markersize=8,
                         linestyle="", zorder=6)

        n_missing = int((counts == 0).sum())
        # Kept as a single line: rotated multi-line tick labels lay out very
        # differently from rotated single-line ones (each line becomes its
        # own column) and need much more room than is practical here.
        extra_bits = [b for b in [f"{len(outliers)} ND" if len(outliers) else "",
                                  f"{n_missing} miss" if n_missing else ""] if b]
        label = sample_label(srow, group)
        cov_labels.append(label)
        even_labels.append(label + (" (" + ", ".join(extra_bits) + ")" if extra_bits else ""))

    ax_cov.set_ylim(0, 125)
    ax_cov.set_xticks(range(len(samples)))
    ax_cov.set_xticklabels(cov_labels, rotation=XTICK_ROTATION, ha="center")
    ax_cov.set_ylabel("% of designed peptides detected")
    ax_cov.set_title(f"Coverage (bar labels: peptides with >= {COVERAGE_MIN_COUNT} reads, of {total} designed)")
    legend_outside(ax_cov, treatment_legend_handles(samples, colors) + [
        Line2D([], [], linestyle="", label="— detected at —"),
        Patch(facecolor="gray", alpha=0.35, edgecolor="black", label=">= 1 read"),
        Patch(facecolor="gray", edgecolor="black", label=f">= {COVERAGE_MIN_COUNT} reads"),
    ], fontsize=8)

    ax_even.set_yscale("log")
    ax_even.set_xlim(-0.6, len(samples) - 0.4)
    ax_even.set_xticks(range(len(samples)))
    ax_even.set_xticklabels(even_labels, fontsize=8, rotation=XTICK_ROTATION, ha="center")
    ax_even.set_ylabel("Reads per peptide (log scale)")
    ax_even.set_title(f"Evenness (10th-90th percentile + outliers)\n(ND = near-dropout, <={OUTLIER_THRESHOLD} reads)")
    legend_outside(ax_even, [
        Line2D([], [], color="gray", linewidth=6, label="typical range (10th-90th percentile)"),
        Line2D([], [], marker="o", linestyle="", color="white", markeredgecolor="black", markersize=7,
               label="median"),
        Line2D([], [], marker="v", linestyle="", color=OUTLIER_COLOR, markersize=8,
               label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)"),
    ], fontsize=8)

    fig.suptitle("Library completeness and evenness overview" + group["title_suffix"], fontsize=13)
    # Fixed widths/heights (in inches, converted to figure fractions) are
    # reserved for the outside legends, the rotated tick labels, and a gap
    # between the two stacked panels (their rotated labels/title would
    # otherwise collide). Not using
    # bbox_inches="tight" here: with this matplotlib version it can blow up
    # the canvas to a pathological size when combined with a log-scaled axis
    # and many plotted markers.
    fig.subplots_adjust(right=1 - 3.2 / fig.get_figwidth(), bottom=TICK_LABEL_HEIGHT_IN / fig.get_figheight(),
                        top=1 - 0.9 / fig.get_figheight(), hspace=0.9)
    fig.savefig(group["fig_dir"] / "overview_all_samples.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {group['fig_dir'] / 'overview_all_samples.png'}", file=sys.stderr)


def main():
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    for group in figure_groups(samples, FIG_DIR):
        plot_overview(mat, lib_of_ref, group, colors)


if __name__ == "__main__":
    main()
