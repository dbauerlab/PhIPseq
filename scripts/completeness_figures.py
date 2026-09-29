"""Plainer, completeness-focused figures to complement summarize.py's
rank-abundance/Lorenz plots:
  1. coverage_bar.png       -- simple % of designed library detected, per sample
  2. rarefaction_curves.png -- did we sequence deep enough to have found
                                everything there is to find? (exact expected
                                unique-peptides-detected vs. reads subsampled,
                                via the standard ecology rarefaction formula)
  3. detection_histogram.png -- read-count distribution per peptide, with the
                                 "never detected" peptides called out explicitly

Usage: completeness_figures.py [samples_tsv]  (defaults to reference/samples.tsv)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import gammaln

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
# with the other figure-generating scripts) rather than being duplicated as
# literals here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
LEGEND_BBOX_ANCHOR = tuple(_fig_params["legend_bbox_anchor"])
XTICK_ROTATION = _fig_params["xtick_rotation"]
DETECTION_HISTOGRAM_NCOLS = _fig_params["detection_histogram_ncols"]


def expected_unique_at_depths(counts, depths):
    """Exact expected number of species with >=1 read when subsampling
    `depth` reads without replacement from the observed pool (Hurlbert 1971
    rarefaction). Computed in log-space via gammaln for numerical stability
    at hundreds-of-millions-of-reads scale."""
    counts = np.asarray(counts, dtype=np.float64)
    present = counts > 0
    Ni = counts[present]
    N = counts.sum()
    out = []
    for n in depths:
        n = min(n, N)
        can_be_absent = (N - Ni) >= n
        log_p_absent = np.zeros_like(Ni)
        log_p_absent[~can_be_absent] = -np.inf  # too few "other" reads -> guaranteed present
        idx = can_be_absent
        log_p_absent[idx] = (
            (gammaln(N - Ni[idx] + 1) - gammaln(N - Ni[idx] - n + 1))
            - (gammaln(N + 1) - gammaln(N - n + 1))
        )
        p_absent = np.exp(log_p_absent)
        out.append(float(np.sum(1 - p_absent)))
    return np.array(out)


def main():
    samples_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    samples = pd.read_csv(samples_path, sep="\t")
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]

    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # --- 1. Coverage bar chart ---
    fig, ax = plt.subplots(figsize=(1.2 * len(samples) + 2, 4))
    labels, pcts, colors = [], [], []
    palette = plt.get_cmap("tab10")
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index
        own_refs = own_refs.intersection(mat.index)
        total = len(own_refs)
        detected = int((mat.loc[own_refs, tc_id] > 0).sum())
        labels.append(tc_id)
        pcts.append(100 * detected / total if total else 0)
        colors.append(palette(i))
    ax.bar(labels, pcts, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=XTICK_ROTATION, ha="center")
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of designed peptides detected (>=1 read)")
    ax.set_title("Library coverage: how complete is each batch?")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "coverage_bar.png", dpi=150)
    plt.close(fig)

    # --- 2. Rarefaction / saturation curves ---
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index
        own_refs = own_refs.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id].values
        total_lib_size = len(own_refs)
        total_reads = counts.sum()
        if total_reads == 0:
            continue
        depths = np.unique(np.round(np.logspace(np.log10(max(total_reads / 1e5, 100)),
                                                  np.log10(total_reads), 12)).astype(np.int64))
        unique_at_depth = expected_unique_at_depths(counts, depths)
        ax.plot(depths, unique_at_depth, marker="o", markersize=3, color=palette(i), label=tc_id)
        ax.axhline(total_lib_size, color=palette(i), linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("Reads subsampled")
    ax.set_ylabel("Expected unique peptides detected")
    ax.set_title("Have we sequenced deep enough? (rarefaction)\ncurve flattening = fully saturated at this depth")
    ax.legend(loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR, borderaxespad=0, fontsize=8)
    fig.savefig(FIG_DIR / "rarefaction_curves.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- 3. Detection histogram with explicit "missing" bin ---
    n_samples = len(samples)
    ncols = min(DETECTION_HISTOGRAM_NCOLS, n_samples)
    nrows = int(np.ceil(n_samples / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), sharey=True, squeeze=False)
    axes_flat = axes.flatten()
    legend_handles = None
    for ax, (_, srow) in zip(axes_flat, samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index
        own_refs = own_refs.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id].values
        n_missing = int((counts == 0).sum())
        nonzero = counts[counts > 0]
        bins = np.logspace(0, np.log10(max(nonzero.max(), 10)), 30) if len(nonzero) else [1]
        ax.bar([0.5], [n_missing], width=0.9, color="firebrick", label="never detected (0 reads)")
        if len(nonzero):
            ax.hist(nonzero, bins=bins, color="steelblue", label="detected peptides")
        ax.set_xscale("symlog", linthresh=1)
        ax.set_title(f"{tc_id}\n{n_missing} missing / {len(own_refs)} designed", fontsize=9)
        ax.set_xlabel("Reads per peptide")
        if legend_handles is None:
            legend_handles, _ = ax.get_legend_handles_labels()
    for row in range(nrows):
        axes_flat[row * ncols].set_ylabel("Number of peptides")
    for ax in axes_flat[n_samples:]:
        ax.axis("off")
    fig.suptitle("Per-peptide detection: how many are missing vs. how well are the rest covered?")
    if legend_handles:
        fig.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.0, 0.5), borderaxespad=0, fontsize=9)
    # rect leaves headroom above the top row of subplot titles for the
    # suptitle -- tight_layout() doesn't know about the suptitle otherwise
    # and the two overlap.
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG_DIR / "detection_histogram.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote coverage_bar.png, rarefaction_curves.png, detection_histogram.png to {FIG_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
