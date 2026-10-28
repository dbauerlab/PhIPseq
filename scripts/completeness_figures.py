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
# --- Imports -----------------------------------------------------------------
import json  # to read the reference/*.json config files
import sys  # sys.stderr for progress messages (CLI args are read inside sample_groups)
from pathlib import Path  # object-oriented file paths ("/" joins path parts)

import numpy as np  # vectorised maths on whole arrays of counts at once
import pandas as pd  # DataFrames: the count matrix and metadata tables
import matplotlib

# Select the non-interactive "Agg" backend BEFORE importing pyplot: it renders
# straight to PNG files and needs no display, so the script runs on a cluster
# node. (This is why the pyplot import sits below a statement.)
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # the plotting interface (figures, axes)
from matplotlib.lines import Line2D  # used to build custom legend entries
from matplotlib.patches import Patch  # filled-rectangle legend entries
from scipy.special import gammaln  # log(Gamma(x)) = log((x-1)!), safe for huge x
from scipy.stats import binom  # binomial distribution (used for its survival function)

# Local sibling modules in scripts/ (importable because the script is run as
# `python scripts/completeness_figures.py`, which puts scripts/ on sys.path).
# COVERAGE_MIN_COUNT = reads needed for a peptide to count as "detected" in the
# stricter coverage metric (filter_params.json "coverage_min_count", 10 at time of writing).
from count_filters import COVERAGE_MIN_COUNT
# Shared sample-sheet loading, figure grouping (all samples / per IP) and
# legend/colour/line-style helpers, so every figure script looks the same.
from sample_groups import (LINESTYLE_CYCLE, figure_groups, legend_outside, load_samples, sample_label,
                           style_legend_handles, style_map, treatment_legend_handles)

# Repo root: this file is scripts/completeness_figures.py, so go up two levels
# (.resolve() first turns a relative path into an absolute one).
PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())  # dict of paths
RUN_DIR = Path(_paths["run_dir"])  # the data directory (may live anywhere)
RESULTS_DIR = RUN_DIR / "results"  # count_matrix.csv is read from here
FIG_DIR = RESULTS_DIR / "figures"  # all-samples figures here; per-IP ones in FIG_DIR/<ip>/
# Per-peptide annotation (ref_id -> library, organism, ...) from build_reference.py.
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

# Figure thresholds/styling live in reference/figure_params.json (shared
# with the other figure-generating scripts) rather than being duplicated as
# literals here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
XTICK_ROTATION = _fig_params["xtick_rotation"]  # angle (degrees) of x-axis sample labels
# Max number of subplot columns in the detection-histogram grid (more samples -> more rows).
DETECTION_HISTOGRAM_NCOLS = _fig_params["detection_histogram_ncols"]

# Neutral colour for the "never detected" bin, kept distinct from the
# treatment palette.
MISSING_COLOR = "dimgray"
# Shading for the 1..(coverage_min_count - 1) read range: detected, but
# below the coverage threshold.
BELOW_THRESHOLD_COLOR = "gold"


# ---------------------------------------------------------------------------
# expected_unique_at_depths(counts, depths) -> 1-D numpy array (one value per depth)
#   In:  counts = read counts for every peptide of one sample's library
#        (array-like); depths = read depths to "subsample" to.
#   Out: expected number of distinct peptides with >= 1 read at each depth.
#   Called by plot_rarefaction() when min_count == 1.
#
# Hurlbert (1971) rarefaction: if we drew n of the N reads at random WITHOUT
# replacement, peptide i (with Ni reads) is missed only if all n reads come
# from the other N - Ni reads. That probability is
#     P(absent_i) = C(N - Ni, n) / C(N, n)
# and the expected richness is E[S_n] = sum_i (1 - P(absent_i)).
# With N in the hundreds of millions the binomial coefficients overflow any
# float, so we work with logs: log C(a, b) = gammaln(a+1) - gammaln(b+1) -
# gammaln(a-b+1). In the ratio the gammaln(n+1) terms cancel, leaving
#     log P = [gammaln(N-Ni+1) - gammaln(N-Ni-n+1)] - [gammaln(N+1) - gammaln(N-n+1)]
# which is what the code below computes, then exponentiates once at the end.
# ---------------------------------------------------------------------------
def expected_unique_at_depths(counts, depths):
    """Exact expected number of species with >=1 read when subsampling
    `depth` reads without replacement from the observed pool (Hurlbert 1971
    rarefaction). Computed in log-space via gammaln for numerical stability
    at hundreds-of-millions-of-reads scale."""
    # float64 because gammaln and -inf need floating point.
    counts = np.asarray(counts, dtype=np.float64)
    present = counts > 0  # boolean mask: peptides seen at least once
    # Ni = counts of observed peptides only; unseen ones can never be drawn
    # (expected contribution 0), so dropping them changes nothing.
    Ni = counts[present]
    N = counts.sum()  # total reads in this sample (the pool we subsample from)
    out = []  # one expected-richness value per depth
    for n in depths:
        n = min(n, N)  # can't subsample more reads than exist
        # Peptide i can be missed only if there are at least n "other" reads to
        # draw from; otherwise C(N - Ni, n) = 0 and it is certainly present.
        can_be_absent = (N - Ni) >= n  # boolean mask, one entry per peptide
        log_p_absent = np.zeros_like(Ni)  # same shape/dtype as Ni, filled below
        log_p_absent[~can_be_absent] = -np.inf  # too few "other" reads -> guaranteed present
        idx = can_be_absent  # mask of peptides where the formula applies
        # Vectorised log P(absent) for all applicable peptides at once
        # (the Hurlbert formula in log space, see block comment above).
        log_p_absent[idx] = (
            (gammaln(N - Ni[idx] + 1) - gammaln(N - Ni[idx] - n + 1))
            - (gammaln(N + 1) - gammaln(N - n + 1))
        )
        p_absent = np.exp(log_p_absent)  # back to probabilities; exp(-inf) = 0
        # Expected richness = sum of each peptide's probability of being seen.
        out.append(float(np.sum(1 - p_absent)))
    return np.array(out)  # array aligned with `depths`


# ---------------------------------------------------------------------------
# expected_at_least_at_depths(counts, depths, min_count) -> 1-D numpy array
#   In:  one sample's per-peptide counts, the depths to evaluate, and the
#        read threshold (COVERAGE_MIN_COUNT in practice).
#   Out: expected number of peptides with >= min_count reads at each depth.
#   Called by plot_rarefaction() when min_count > 1.
#
# There is no tidy closed form like Hurlbert's for ">= k reads", so we use
# binomial thinning: subsampling to depth n is approximated by keeping each
# read independently with probability p = n / N. A peptide with Ni reads then
# keeps Binomial(Ni, p) reads, and P(>= min_count) is the binomial survival
# function sf(min_count - 1) = P(X > min_count - 1). Summing over peptides
# gives the expected count. At p = 1 every peptide keeps all its reads, so the
# value is exact at full depth.
# ---------------------------------------------------------------------------
def expected_at_least_at_depths(counts, depths, min_count):
    """Expected number of peptides with >= min_count reads when keeping each
    read independently with probability depth / total (binomial thinning --
    the standard approximation to subsampling without replacement, exact at
    full depth)."""
    # Integer counts: binom needs whole-number trial counts.
    counts = np.asarray(counts, dtype=np.int64)
    # Thinning only ever removes reads, so peptides already below min_count
    # can never reach it -- keep just the ones that can.
    Ni = counts[counts >= min_count]
    N = counts.sum()  # total reads in this sample
    # For each depth n: p = n/N (capped at 1); binom.sf broadcasts over the
    # whole Ni array, giving P(>= min_count) per peptide; .sum() = expectation.
    return np.array([float(binom.sf(min_count - 1, Ni, min(n / N, 1.0)).sum()) for n in depths])


# ---------------------------------------------------------------------------
# own_library_counts(mat, lib_of_ref, srow) -> pandas Series (ref_id -> count)
#   In:  mat = full peptides x samples count matrix; lib_of_ref = Series
#        mapping ref_id -> library name; srow = one row of the sample sheet.
#   Out: that sample's counts restricted to peptides of its own library, so
#        coverage is measured against the library it was actually screened
#        with. Called by every plot_* function below.
# ---------------------------------------------------------------------------
def own_library_counts(mat, lib_of_ref, srow):
    """One sample's counts over its own library's peptides."""
    # Boolean mask on lib_of_ref keeps this sample's library; .index gives the
    # ref_ids; .intersection keeps only ids that are actually rows of mat.
    own_refs = lib_of_ref[lib_of_ref == srow["library"]].index.intersection(mat.index)
    return mat.loc[own_refs, srow["tc_id"]]  # rows = own peptides, column = this sample


# ---------------------------------------------------------------------------
# plot_coverage(mat, lib_of_ref, group, colors) -> None (writes a PNG)
#   In:  count matrix, ref_id -> library map, one figure group (dict from
#        sample_groups.figure_groups) and the treatment -> colour dict.
#   Out: <group fig_dir>/coverage_bar.png -- one bar per sample showing the %
#        of its library detected at >= 1 read (pale bar) with the % at
#        >= COVERAGE_MIN_COUNT reads (solid bar) drawn over it.
#   Called by main() once per figure group.
# ---------------------------------------------------------------------------
def plot_coverage(mat, lib_of_ref, group, colors):
    samples = group["samples"]  # sample-sheet rows belonging to this group
    # Width grows with the number of samples (1.2 inch per bar + margin); 4 in tall.
    fig, ax = plt.subplots(figsize=(1.2 * len(samples) + 2, 4))
    labels, pcts_any, pcts_min = [], [], []  # x labels, % at >=1 read, % at >=threshold
    for _, srow in samples.iterrows():  # iterrows yields (index, row Series) pairs
        counts = own_library_counts(mat, lib_of_ref, srow)  # this sample's own-library counts
        labels.append(sample_label(srow, group))  # tick label (ip-prefixed in the "all" group)
        # Mean of a boolean mask = fraction of True values; x100 = percent.
        # Guard against an empty library (mean of empty -> NaN) by using 0.
        pcts_any.append(100 * (counts > 0).mean() if len(counts) else 0)
        pcts_min.append(100 * (counts >= COVERAGE_MIN_COUNT).mean() if len(counts) else 0)
    bar_colors = [colors[t] for t in samples["treatment"]]  # one colour per bar, by treatment
    # Two bar layers at the same x positions: a translucent (alpha=0.35) ">=1
    # read" bar, then the opaque ">= threshold" bar on top. The >=threshold
    # value is always <= the >=1 value, so the pale part shows the difference.
    ax.bar(range(len(labels)), pcts_any, color=bar_colors, alpha=0.35, edgecolor="black", linewidth=0.5)
    ax.bar(range(len(labels)), pcts_min, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))  # one tick per bar...
    ax.set_xticklabels(labels, rotation=XTICK_ROTATION, ha="center")  # ...labelled by sample
    ax.set_ylim(0, 105)  # a little headroom above 100%
    ax.set_ylabel("% of designed peptides detected")
    # title_suffix is e.g. " (all samples)" or " (IP4)".
    ax.set_title("Library coverage: how complete is each batch?" + group["title_suffix"])
    # Extra legend entries explaining pale vs solid bars. The first is an
    # invisible Line2D (no line, no marker) used purely as a text header.
    threshold_handles = [Line2D([], [], linestyle="", label="— detected at —"),
                         Patch(facecolor="gray", alpha=0.35, edgecolor="black", label=">= 1 read"),
                         Patch(facecolor="gray", edgecolor="black", label=f">= {COVERAGE_MIN_COUNT} reads")]
    # Treatment colour key + the threshold key, placed outside the axes on the right.
    legend_outside(ax, treatment_legend_handles(samples, colors) + threshold_handles, fontsize=8)
    # bbox_inches="tight" grows the saved image so the outside legend isn't cropped.
    fig.savefig(group["fig_dir"] / "coverage_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)  # free the figure's memory (many figures are made in one run)


# ---------------------------------------------------------------------------
# plot_rarefaction(mat, lib_of_ref, group, colors, min_count=1) -> None
#   In:  as plot_coverage, plus min_count (1 = any read; >1 = stricter threshold).
#   Out: <group fig_dir>/rarefaction_curves.png (min_count == 1) or
#        rarefaction_curves_min<min_count>.png -- one curve per sample of
#        expected peptides detected vs. reads subsampled (log x axis). A curve
#        that flattens out means more sequencing would find few new peptides.
#   Called by main() twice per group (min_count=1 and min_count=COVERAGE_MIN_COUNT).
# ---------------------------------------------------------------------------
def plot_rarefaction(mat, lib_of_ref, group, colors, min_count=1):
    """Rarefaction curves for peptides reaching >= min_count reads: exact
    (Hurlbert) for min_count=1, binomial thinning otherwise."""
    # Map each value of the group's style column (ip in "all", replicate per
    # IP) to a line style ("-", "--", ...), so colour = treatment, dash = ip/replicate.
    linestyles = style_map(group, LINESTYLE_CYCLE)
    fig, ax = plt.subplots(figsize=(7, 5))  # a single 7x5-inch panel
    lib_sizes = set()  # distinct library sizes seen -> dotted reference lines
    for _, srow in group["samples"].iterrows():
        # .values -> plain numpy array (the maths functions don't need the index).
        counts = own_library_counts(mat, lib_of_ref, srow).values
        total_reads = counts.sum()  # this sample's own-library read total
        if total_reads == 0:
            continue  # nothing to rarefy (and would divide by zero); skip sample
        lib_sizes.add(len(counts))  # number of own-library peptides = max possible richness
        # 12 depths evenly spaced on a log scale, from max(total/100,000, 100)
        # reads up to the full total. np.logspace takes log10 endpoints. Depths
        # are rounded to whole reads, and np.unique drops duplicates (and sorts).
        depths = np.unique(np.round(np.logspace(np.log10(max(total_reads / 1e5, 100)),
                                                  np.log10(total_reads), 12)).astype(np.int64))
        # Exact Hurlbert formula for ">=1 read"; binomial-thinning approximation otherwise.
        if min_count == 1:
            unique_at_depth = expected_unique_at_depths(counts, depths)
        else:
            unique_at_depth = expected_at_least_at_depths(counts, depths, min_count)
        # One line per sample: colour by treatment, line style by ip/replicate.
        ax.plot(depths, unique_at_depth, marker="o", markersize=3, color=colors[srow["treatment"]],
                linestyle=linestyles[srow[group["style_col"]]], linewidth=1.2)
    # Dotted horizontal line at each library size: the ceiling a curve could reach.
    for size in lib_sizes:
        ax.axhline(size, color="black", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_xscale("log")  # depths span several orders of magnitude
    ax.set_xlabel("Reads subsampled")
    # Add the threshold to the y label only for the min_count > 1 version.
    ax.set_ylabel("Expected unique peptides detected" + (f" (>= {min_count} reads)" if min_count > 1 else ""))
    # Title: "read" vs "reads" pluralised by min_count, then the group suffix,
    # then a second line explaining how to read the plot.
    ax.set_title(f"Have we sequenced deep enough? (rarefaction, >= {min_count} read{'s' if min_count > 1 else ''})"
                 + group["title_suffix"]
                 + "\ncurve flattening = fully saturated at this depth; dotted = designed library size")
    # Legend: treatment colours (as lines) followed by the line-style key.
    legend_outside(ax, treatment_legend_handles(group["samples"], colors, kind="line")
                   + style_legend_handles(group, linestyles, kind="line"), fontsize=8)
    # File name depends on min_count, e.g. rarefaction_curves_min10.png.
    name = "rarefaction_curves.png" if min_count == 1 else f"rarefaction_curves_min{min_count}.png"
    fig.savefig(group["fig_dir"] / name, dpi=150, bbox_inches="tight")
    plt.close(fig)  # release memory


# ---------------------------------------------------------------------------
# plot_detection_histogram(mat, lib_of_ref, group, colors) -> None
#   In:  as plot_coverage.
#   Out: <group fig_dir>/detection_histogram.png -- a grid with one subplot per
#        sample: a grey bar for peptides with 0 reads, a log-binned histogram
#        of reads per detected peptide, and a gold band over the 1 to
#        (COVERAGE_MIN_COUNT - 1) read range (detected but below threshold).
#   Called by main() once per figure group.
# ---------------------------------------------------------------------------
def plot_detection_histogram(mat, lib_of_ref, group, colors):
    samples = group["samples"]  # this group's sample-sheet rows
    n_samples = len(samples)  # one subplot per sample
    ncols = min(DETECTION_HISTOGRAM_NCOLS, n_samples)  # don't make empty columns for small groups
    nrows = int(np.ceil(n_samples / ncols))  # enough rows to fit every sample (round up)
    # A nrows x ncols grid of 4x4-inch panels. sharey=True puts every panel on
    # the same y scale so heights are comparable; squeeze=False always returns
    # a 2-D array of axes, even for a 1x1 or 1xN grid.
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), sharey=True, squeeze=False)
    axes_flat = axes.flatten()  # 1-D list of axes, filled row by row
    # zip pairs each axes with one sample; it stops at the shorter sequence,
    # leaving any spare axes at the end untouched.
    for ax, (_, srow) in zip(axes_flat, samples.iterrows()):
        counts = own_library_counts(mat, lib_of_ref, srow).values  # numpy array of counts
        n_missing = int((counts == 0).sum())  # summing a boolean mask counts the Trues
        # Detected (>0) AND below threshold; "&" combines masks element-wise.
        n_below = int(((counts > 0) & (counts < COVERAGE_MIN_COUNT)).sum())
        nonzero = counts[counts > 0]  # counts of detected peptides only
        # 30 log-spaced bin edges from 10^0 = 1 read up to the max count (at
        # least 10, so the range is never degenerate). If nothing was detected
        # a dummy [1] is used (the histogram is skipped below anyway).
        bins = np.logspace(0, np.log10(max(nonzero.max(), 10)), 30) if len(nonzero) else [1]
        # The "never detected" peptides as a single bar centred at x = 0.5,
        # i.e. left of the first histogram bin (which starts at 1 read).
        ax.bar([0.5], [n_missing], width=0.9, color=MISSING_COLOR)
        if len(nonzero):
            # Histogram of reads per detected peptide, coloured by treatment.
            ax.hist(nonzero, bins=bins, color=colors[srow["treatment"]], edgecolor="black", linewidth=0.3)
        # Shade x from 1 to COVERAGE_MIN_COUNT (detected, but below threshold);
        # zorder=0 draws it behind the bars.
        ax.axvspan(1, COVERAGE_MIN_COUNT, color=BELOW_THRESHOLD_COLOR, alpha=0.3, zorder=0)
        # "symlog" = linear between -1 and 1 (linthresh=1) and log beyond, so
        # the 0-read bar at x = 0.5 can sit on the same axis as the log-scale
        # histogram (a pure log axis cannot show values <= 0).
        ax.set_xscale("symlog", linthresh=1)
        # Two-line subplot title: sample + treatment, then the counts summary.
        ax.set_title(f"{sample_label(srow, group)} - {srow['treatment']}\n"
                     f"{n_missing} missing, {n_below} below {COVERAGE_MIN_COUNT} / {len(counts)} designed",
                     fontsize=9)
        ax.set_xlabel("Reads per peptide")
    # y label only on the first (left-most) panel of each row; row * ncols is
    # that panel's position in the flattened list.
    for row in range(nrows):
        axes_flat[row * ncols].set_ylabel("Number of peptides")
    # Hide unused panels in the last row (grid slots beyond the sample count).
    for ax in axes_flat[n_samples:]:
        ax.axis("off")
    # Figure-level title above all panels.
    fig.suptitle("Per-peptide detection: how many are missing vs. how well are the rest covered?"
                 + group["title_suffix"])
    # One shared legend for the whole figure: the 0-read bar, the shaded band,
    # an invisible-line text header, then one entry per treatment colour.
    handles = [Patch(color=MISSING_COLOR, label="never detected (0 reads)"),
               Patch(color=BELOW_THRESHOLD_COLOR, alpha=0.3,
                     label=f"below coverage threshold (1-{COVERAGE_MIN_COUNT - 1} reads)"),
               Line2D([], [], linestyle="", label="detected peptides, by")] + treatment_legend_handles(samples, colors)
    # Figure legend anchored just right of the figure (x = 1.0), centred vertically.
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), borderaxespad=0, fontsize=9)
    # rect leaves headroom above the top row of subplot titles for the
    # suptitle -- tight_layout() doesn't know about the suptitle otherwise
    # and the two overlap.
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    # bbox_inches="tight" also captures the legend that sits outside the figure.
    fig.savefig(group["fig_dir"] / "detection_histogram.png", dpi=150, bbox_inches="tight")
    plt.close(fig)  # release memory


# ---------------------------------------------------------------------------
# main() -- entry point: load inputs, then draw every figure for every group.
#   In:  optional CLI args (sample sheet, sample metadata; read by
#        sample_groups.load_samples), combined_metadata.csv and
#        results/count_matrix.csv (written by summarize.py).
#   Out: the four PNGs per figure group (all samples + one per IP).
# ---------------------------------------------------------------------------
def main():
    # samples: merged, plot-ordered sample sheet; colors: treatment -> colour.
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)  # dtype=str: keep every column as text
    # Series indexed by ref_id whose values are library names (ref_id -> library).
    lib_of_ref = meta.set_index("ref_id")["library"]
    # Raw (unfiltered) peptides x samples counts; first column (ref_id) becomes
    # the index. Coverage must use the raw matrix: the filtered one has
    # already dropped the undetected peptides.
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    # figure_groups yields the "all samples" group then one group per IP, each
    # with its own output folder (created if missing).
    for group in figure_groups(samples, FIG_DIR):
        plot_coverage(mat, lib_of_ref, group, colors)
        plot_rarefaction(mat, lib_of_ref, group, colors)  # >= 1 read (exact)
        plot_rarefaction(mat, lib_of_ref, group, colors, min_count=COVERAGE_MIN_COUNT)  # thinning
        plot_detection_histogram(mat, lib_of_ref, group, colors)
        # Progress message on stderr. "{{" / "}}" are literal braces in an
        # f-string, so this prints shell-style "rarefaction_curves{,_min10}.png".
        print(f"Wrote coverage_bar.png, rarefaction_curves{{,_min{COVERAGE_MIN_COUNT}}}.png, detection_histogram.png "
              f"to {group['fig_dir']}", file=sys.stderr)


# Run main() only when executed as a script, not when imported.
if __name__ == "__main__":
    main()
