"""One figure summarizing library evenness per the plain-language stats
already discussed: for each sample, the P10-P90 "typical range" of
per-peptide read counts, the median, and any near-dropout outliers (<=10
reads) called out explicitly rather than buried in a max/min ratio.

Usage: evenness_summary_figure.py [samples_tsv]  (defaults to scripts/samples.tsv)
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
SAMPLES_TSV = PROJECT_DIR / "scripts" / "samples.tsv"

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


def main():
    samples_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    samples = pd.read_csv(samples_path, sep="\t")
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    fig, ax = plt.subplots(figsize=(1.4 * len(samples) + 2, 5.5))
    palette = plt.get_cmap("tab10")
    xtick_labels = []

    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]
        detected = counts[counts > 0]
        if detected.empty:
            continue
        p10, p50, p90 = np.percentile(detected, [10, 50, 90])
        color = palette(i)

        # the tight "typical" band most peptides fall in
        ax.plot([i, i], [p10, p90], color=color, linewidth=6, solid_capstyle="butt", alpha=0.85,
                label="typical range (10th-90th percentile)" if i == 0 else None)
        ax.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=8, zorder=5,
                label="median" if i == 0 else None)

        # near-dropout outliers, called out individually rather than folded into a min/max ratio
        outliers = detected[detected <= OUTLIER_THRESHOLD]
        n_missing = int((counts == 0).sum())
        for val in outliers.values:
            ax.plot(i, val, "v", color="firebrick", markersize=9, zorder=6,
                     label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)" if i == 0 and val == outliers.values[0] else None)
        label_bits = []
        if len(outliers):
            label_bits.append(f"{len(outliers)} near-dropout")
        if n_missing:
            label_bits.append(f"{n_missing} missing")
        extra = "\n" + ", ".join(label_bits) if label_bits else "\nall clean"
        xtick_labels.append(f"{tc_id}{extra}")

    ax.set_yscale("log")
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(xtick_labels)
    ax.set_ylabel("Reads per peptide (log scale)")
    ax.set_title("Library evenness: how tightly clustered is representation?")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper right", fontsize=8)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "evenness_summary.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIG_DIR / 'evenness_summary.png'}", file=sys.stderr)


if __name__ == "__main__":
    main()
