"""Plainer, completeness-focused figures to complement summarize.py's
rank-abundance/Lorenz plots, coloured by treatment and drawn once for all
samples (results/figures/) and once per IP condition
(results/figures/<ip>/, see sample_groups.py):
  1. coverage_bar.png       -- % of designed library detected, per sample, at
                                >=1 read (light bar) and at >= coverage_min_count
                                reads (solid bar; see count_filters.py)
  2. rarefaction_curves.png -- did we sequence deep enough to have found
                                everything there is to find? (exact expected
                                unique-peptides-detected vs. reads subsampled,
                                via the standard ecology rarefaction formula)
     rarefaction_curves_min<N>.png -- the same for peptides reaching
                                >= coverage_min_count reads
  3. detection_histogram.png -- read-count distribution per peptide, with the
                                 "never detected" peptides called out explicitly
                                 and the below-threshold range shaded

Usage: completeness_figures.py [samples_tsv] [samples_meta_csv]
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
from scipy.special import gammaln
from scipy.stats import binom

from count_filters import COVERAGE_MIN_COUNT
from sample_groups import (LINESTYLE_CYCLE, figure_groups, legend_outside, load_samples, sample_label,
                           style_legend_handles, style_map, treatment_legend_handles)

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
# with the other figure-generating scripts) rather than being duplicated as
# literals here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
XTICK_ROTATION = _fig_params["xtick_rotation"]
DETECTION_HISTOGRAM_NCOLS = _fig_params["detection_histogram_ncols"]

# Neutral colour for the "never detected" bin, kept distinct from the
# treatment palette.
MISSING_COLOR = "dimgray"
# Shading for the 1..(coverage_min_count - 1) read range: detected, but
# below the coverage threshold.
BELOW_THRESHOLD_COLOR = "gold"


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


def expected_at_least_at_depths(counts, depths, min_count):
    """Expected number of peptides with >= min_count reads when keeping each
    read independently with probability depth / total (binomial thinning --
    the standard approximation to subsampling without replacement, exact at
    full depth)."""
    counts = np.asarray(counts, dtype=np.int64)
    Ni = counts[counts >= min_count]
    N = counts.sum()
    return np.array([float(binom.sf(min_count - 1, Ni, min(n / N, 1.0)).sum()) for n in depths])


def own_library_counts(mat, lib_of_ref, srow):
    """One sample's counts over its own library's peptides."""
    own_refs = lib_of_ref[lib_of_ref == srow["library"]].index.intersection(mat.index)
    return mat.loc[own_refs, srow["tc_id"]]


def plot_coverage(mat, lib_of_ref, group, colors):
    samples = group["samples"]
    fig, ax = plt.subplots(figsize=(1.2 * len(samples) + 2, 4))
    labels, pcts_any, pcts_min = [], [], []
    for _, srow in samples.iterrows():
        counts = own_library_counts(mat, lib_of_ref, srow)
        labels.append(sample_label(srow, group))
        pcts_any.append(100 * (counts > 0).mean() if len(counts) else 0)
        pcts_min.append(100 * (counts >= COVERAGE_MIN_COUNT).mean() if len(counts) else 0)
    bar_colors = [colors[t] for t in samples["treatment"]]
    ax.bar(range(len(labels)), pcts_any, color=bar_colors, alpha=0.35, edgecolor="black", linewidth=0.5)
    ax.bar(range(len(labels)), pcts_min, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=XTICK_ROTATION, ha="center")
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of designed peptides detected")
    ax.set_title("Library coverage: how complete is each batch?" + group["title_suffix"])
    threshold_handles = [Line2D([], [], linestyle="", label="— detected at —"),
                         Patch(facecolor="gray", alpha=0.35, edgecolor="black", label=">= 1 read"),
                         Patch(facecolor="gray", edgecolor="black", label=f">= {COVERAGE_MIN_COUNT} reads")]
    legend_outside(ax, treatment_legend_handles(samples, colors) + threshold_handles, fontsize=8)
    fig.savefig(group["fig_dir"] / "coverage_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_rarefaction(mat, lib_of_ref, group, colors, min_count=1):
    """Rarefaction curves for peptides reaching >= min_count reads: exact
    (Hurlbert) for min_count=1, binomial thinning otherwise."""
    linestyles = style_map(group, LINESTYLE_CYCLE)
    fig, ax = plt.subplots(figsize=(7, 5))
    lib_sizes = set()
    for _, srow in group["samples"].iterrows():
        counts = own_library_counts(mat, lib_of_ref, srow).values
        total_reads = counts.sum()
        if total_reads == 0:
            continue
        lib_sizes.add(len(counts))
        depths = np.unique(np.round(np.logspace(np.log10(max(total_reads / 1e5, 100)),
                                                  np.log10(total_reads), 12)).astype(np.int64))
        if min_count == 1:
            unique_at_depth = expected_unique_at_depths(counts, depths)
        else:
            unique_at_depth = expected_at_least_at_depths(counts, depths, min_count)
        ax.plot(depths, unique_at_depth, marker="o", markersize=3, color=colors[srow["treatment"]],
                linestyle=linestyles[srow[group["style_col"]]], linewidth=1.2)
    for size in lib_sizes:
        ax.axhline(size, color="black", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("Reads subsampled")
    ax.set_ylabel("Expected unique peptides detected" + (f" (>= {min_count} reads)" if min_count > 1 else ""))
    ax.set_title(f"Have we sequenced deep enough? (rarefaction, >= {min_count} read{'s' if min_count > 1 else ''})"
                 + group["title_suffix"]
                 + "\ncurve flattening = fully saturated at this depth; dotted = designed library size")
    legend_outside(ax, treatment_legend_handles(group["samples"], colors, kind="line")
                   + style_legend_handles(group, linestyles, kind="line"), fontsize=8)
    name = "rarefaction_curves.png" if min_count == 1 else f"rarefaction_curves_min{min_count}.png"
    fig.savefig(group["fig_dir"] / name, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_detection_histogram(mat, lib_of_ref, group, colors):
    samples = group["samples"]
    n_samples = len(samples)
    ncols = min(DETECTION_HISTOGRAM_NCOLS, n_samples)
    nrows = int(np.ceil(n_samples / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), sharey=True, squeeze=False)
    axes_flat = axes.flatten()
    for ax, (_, srow) in zip(axes_flat, samples.iterrows()):
        counts = own_library_counts(mat, lib_of_ref, srow).values
        n_missing = int((counts == 0).sum())
        n_below = int(((counts > 0) & (counts < COVERAGE_MIN_COUNT)).sum())
        nonzero = counts[counts > 0]
        bins = np.logspace(0, np.log10(max(nonzero.max(), 10)), 30) if len(nonzero) else [1]
        ax.bar([0.5], [n_missing], width=0.9, color=MISSING_COLOR)
        if len(nonzero):
            ax.hist(nonzero, bins=bins, color=colors[srow["treatment"]], edgecolor="black", linewidth=0.3)
        ax.axvspan(1, COVERAGE_MIN_COUNT, color=BELOW_THRESHOLD_COLOR, alpha=0.3, zorder=0)
        ax.set_xscale("symlog", linthresh=1)
        ax.set_title(f"{sample_label(srow, group)} - {srow['treatment']}\n"
                     f"{n_missing} missing, {n_below} below {COVERAGE_MIN_COUNT} / {len(counts)} designed",
                     fontsize=9)
        ax.set_xlabel("Reads per peptide")
    for row in range(nrows):
        axes_flat[row * ncols].set_ylabel("Number of peptides")
    for ax in axes_flat[n_samples:]:
        ax.axis("off")
    fig.suptitle("Per-peptide detection: how many are missing vs. how well are the rest covered?"
                 + group["title_suffix"])
    handles = [Patch(color=MISSING_COLOR, label="never detected (0 reads)"),
               Patch(color=BELOW_THRESHOLD_COLOR, alpha=0.3,
                     label=f"below coverage threshold (1-{COVERAGE_MIN_COUNT - 1} reads)"),
               Line2D([], [], linestyle="", label="detected peptides, by")] + treatment_legend_handles(samples, colors)
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), borderaxespad=0, fontsize=9)
    # rect leaves headroom above the top row of subplot titles for the
    # suptitle -- tight_layout() doesn't know about the suptitle otherwise
    # and the two overlap.
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(group["fig_dir"] / "detection_histogram.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    for group in figure_groups(samples, FIG_DIR):
        plot_coverage(mat, lib_of_ref, group, colors)
        plot_rarefaction(mat, lib_of_ref, group, colors)
        plot_rarefaction(mat, lib_of_ref, group, colors, min_count=COVERAGE_MIN_COUNT)
        plot_detection_histogram(mat, lib_of_ref, group, colors)
        print(f"Wrote coverage_bar.png, rarefaction_curves{{,_min{COVERAGE_MIN_COUNT}}}.png, detection_histogram.png "
              f"to {group['fig_dir']}", file=sys.stderr)


if __name__ == "__main__":
    main()
