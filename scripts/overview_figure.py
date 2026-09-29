"""One combined figure comparing all samples together: coverage (%% of
designed peptides detected) and evenness (typical range + near-dropout
outliers) side by side, colored by library so the two T7CoV vs T7Vir3.2
batches are visually grouped.

Usage: overview_figure.py [samples_tsv]  (defaults to reference/samples.tsv,
i.e. all 6 samples once the full run has completed)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_DIR = Path(__file__).resolve().parent.parent
SAMPLES_TSV = PROJECT_DIR / "reference" / "samples.tsv"

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
LIBRARY_COLORS = _fig_params["library_colors"]
LEGEND_BBOX_ANCHOR = tuple(_fig_params["legend_bbox_anchor"])
XTICK_ROTATION = _fig_params["xtick_rotation"]


def main():
    samples_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    samples = pd.read_csv(samples_path, sep="\t")
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    fig, (ax_cov, ax_even) = plt.subplots(2, 1, figsize=(0.7 * len(samples) + 5, 11))

    xtick_labels = []
    seen_libs = []
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        color = LIBRARY_COLORS.get(library, "gray")
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]
        total = len(own_refs)
        detected = counts[counts > 0]
        pct = 100 * len(detected) / total if total else 0

        # --- coverage panel ---
        ax_cov.bar(i, pct, color=color, label=library if library not in seen_libs else None)
        ax_cov.text(i, pct + 1.5, f"{len(detected)}/{total}", ha="center", va="bottom", fontsize=8)

        # --- evenness panel ---
        if len(detected):
            p10, p50, p90 = np.percentile(detected, [10, 50, 90])
            ax_even.plot([i, i], [p10, p90], color=color, linewidth=6, solid_capstyle="butt", alpha=0.85,
                         label=library if library not in seen_libs else None)
            ax_even.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=7, zorder=5)
            outliers = detected[detected <= OUTLIER_THRESHOLD]
            for val in outliers.values:
                ax_even.plot(i, val, "v", color="firebrick", markersize=8, zorder=6,
                             label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)" if "outlier" not in seen_libs else None)
                seen_libs.append("outlier")

        seen_libs.append(library)
        n_missing = int((counts == 0).sum())
        n_outlier = int((detected <= OUTLIER_THRESHOLD).sum()) if len(detected) else 0
        # Kept as a single line: rotated multi-line tick labels lay out very
        # differently from rotated single-line ones (each line becomes its
        # own column) and need much more room than is practical here.
        extra_bits = [b for b in [f"{n_outlier} ND" if n_outlier else "", f"{n_missing} miss" if n_missing else ""] if b]
        xtick_labels.append(tc_id + (" (" + ", ".join(extra_bits) + ")" if extra_bits else ""))

    ax_cov.set_ylim(0, 108)
    ax_cov.set_xticks(range(len(samples)))
    ax_cov.set_xticklabels(samples["tc_id"], rotation=XTICK_ROTATION, ha="center")
    ax_cov.set_ylabel("% of designed peptides detected")
    ax_cov.set_title("Coverage")
    ax_cov.legend(loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR, borderaxespad=0, fontsize=8)

    ax_even.set_yscale("log")
    ax_even.set_xlim(-0.6, len(samples) - 0.4)
    ax_even.set_xticks(range(len(samples)))
    ax_even.set_xticklabels(xtick_labels, fontsize=8, rotation=XTICK_ROTATION, ha="center")
    ax_even.set_ylabel("Reads per peptide (log scale)")
    ax_even.set_title("Evenness (10th-90th percentile + outliers)\n(ND = near-dropout, <=10 reads)")
    ax_even.legend(loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR, borderaxespad=0, fontsize=8)

    fig.suptitle("All samples: library completeness and evenness overview", fontsize=13)
    # Fixed fractions of the figure size are reserved for the outside
    # legends, the rotated tick labels, and a gap between the two stacked
    # panels (their rotated labels/title would otherwise collide). Not using
    # bbox_inches="tight" here: with this matplotlib version it can blow up
    # the canvas to a pathological size when combined with a log-scaled axis
    # and many plotted markers.
    fig.subplots_adjust(right=1 - 2.2 / fig.get_figwidth(), bottom=0.22, top=0.92, hspace=0.9)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "overview_all_samples.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIG_DIR / 'overview_all_samples.png'}", file=sys.stderr)


if __name__ == "__main__":
    main()
