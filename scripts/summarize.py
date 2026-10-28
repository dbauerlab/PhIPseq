"""Aggregate per-sample counts into a peptide x sample matrix and a QC report:
coverage, dropouts, evenness, cross-library mapping, batch reproducibility.

Reads <run_dir>/results/<TC_ID>/counts.tsv (one paired-end fragment count
per peptide, produced by run_sample.sh's merged-lane, paired-mode
cutadapt/bowtie2 pipeline) + its cutadapt.log/bowtie2.log,
<run_dir>/reference/combined_metadata.csv, and the per-sample metadata
sheet (reference/samples_meta.csv by default: sample, ip, treatment).
Figures are coloured by treatment and drawn once for all samples and
once per IP condition (see sample_groups.py). Low-count filtering (see
count_filters.py): QC coverage columns are reported at >=1 read and at
coverage_min_count reads; the PCA and correlation heatmaps apply the
matrix filter within each figure group, on that group's samples only.
Writes, all under <run_dir>:
  results/count_matrix.csv          (raw, unfiltered)
  results/count_matrix_filtered.csv (peptides passing the matrix filter)
  results/filter_report.csv         (per-sample reads removed by it)
  results/figures/{,<ip>/}pca_scores.csv (per-group PCA scores)
  results/qc_summary.csv
  results/figures/*.png        (all samples)
  results/figures/<ip>/*.png   (one folder per IP condition)
  results/summary.md
<run_dir> is configured in reference/paths.json, separate from this repo
checkout.

Usage: summarize.py [samples_tsv] [samples_meta_csv]
  samples_tsv       defaults to reference/samples.tsv
  samples_meta_csv  defaults to reference/samples_meta.csv
"""
# Standard library: JSON config files, regular expressions (log parsing),
# sys (stderr printing / argv, the latter read inside sample_groups.load_samples),
# and pathlib for OS-independent path handling.
import json
import re
import sys
from pathlib import Path

# Third-party: numpy for array maths, pandas for tables, matplotlib for figures.
import numpy as np
import pandas as pd
import matplotlib

# Select the non-interactive "Agg" backend (renders straight to PNG files, no
# display needed). Must happen before pyplot is imported, hence the odd ordering.
matplotlib.use("Agg")
import matplotlib.colors as mcolors  # colour-conversion helpers (e.g. name -> RGBA tuple)
import matplotlib.pyplot as plt  # the usual figure/axes API
# make_axes_locatable lets us glue extra thin axes (colour strips, a colorbar)
# onto the sides of an existing axes so they line up exactly with it.
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Sibling modules in scripts/ (importable because the script is run as
# `python scripts/summarize.py`, which puts scripts/ on sys.path).
# count_filters: the low-count thresholds from reference/filter_params.json
# (coverage_min_count=10, matrix_filter_min_count=10, matrix_filter_min_samples=2),
# human-readable descriptions of them, and the filter functions themselves.
from count_filters import (COVERAGE_MIN_COUNT, FILTER_DESCRIPTION, GROUP_FILTER_DESCRIPTION,
                           MATRIX_FILTER_MIN_COUNT, MATRIX_FILTER_MIN_SAMPLES, filter_matrix, passes_filter)
# sample_groups: loads/merges the sample sheets, fixes the treatment -> colour map,
# defines the figure groups (all samples + one per IP) and shared legend helpers.
from sample_groups import (LINESTYLE_CYCLE, MARKER_CYCLE, figure_groups, legend_outside, load_samples,
                           sample_label, style_legend_handles, style_map, treatment_legend_handles)

# Repo checkout root: this file is scripts/summarize.py, so go up two levels.
# .resolve() turns it into an absolute path (following symlinks) first.
PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())  # whole file -> dict
RUN_DIR = Path(_paths["run_dir"])  # e.g. .../phipseq_lib2
RESULTS_DIR = RUN_DIR / "results"  # per-sample folders + all aggregated outputs
FIG_DIR = RESULTS_DIR / "figures"  # all-samples figures; per-IP ones go in FIG_DIR/<ip>/
# Peptide metadata written by build_reference.py: one row per reference peptide,
# with columns incl. ref_id (e.g. "Vir3|12345"), library ("CoV"/"Vir3"), organism...
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

# Figure styling shared across all figure-generating scripts.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
XTICK_ROTATION = _fig_params["xtick_rotation"]  # degrees to rotate bar-chart sample labels (90)
HEATMAP_TICK_FONTSIZE = _fig_params["heatmap_tick_fontsize"]  # font size of heatmap labels (12)
# How many of the most variable peptides feed the PCA (5000).
PCA_N_TOP_VARIABLE_PEPTIDES = _fig_params["pca_n_top_variable_peptides"]


# load_counts: read one sample's counts.tsv.
#   In:  tc_id -- sample ID, i.e. the name of its folder under results/.
#   Out: dict {ref_id: read-pair count}; only peptides seen at least once are present.
#        An empty dict if the file does not exist (sample not run / failed).
#   Used by build_matrix() to assemble the peptide x sample matrix.
# counts.tsv is written by run_sample.sh as "<ref_id>\t<count>" lines, e.g.
#   "CoV|BCA25643.133.BC1\t61"
def load_counts(tc_id):
    path = RESULTS_DIR / tc_id / "counts.tsv"
    # A missing file is silently treated as "no reads at all" (all-zero column).
    if not path.exists():
        return {}
    counts = {}
    with open(path) as f:  # context manager: file is closed automatically afterwards
        for line in f:  # iterating a file object yields one line at a time
            # Strip the newline, split on the tab, and unpack the two fields at once.
            ref_id, n = line.rstrip("\n").split("\t")
            counts[ref_id] = int(n)  # counts are text in the file -> convert to int
    return counts


# build_matrix: assemble every sample's counts into one table.
#   In:  samples -- the merged sample sheet (one row per sample, column tc_id);
#        meta    -- combined_metadata.csv as a DataFrame (column ref_id).
#   Out: DataFrame, rows = every designed peptide (ref_id), columns = tc_ids,
#        values = read-pair counts (0 where a peptide was never seen).
#   Used by main(); written out as count_matrix.csv and used by everything else.
def build_matrix(samples, meta):
    """One peptide x sample matrix of paired-end fragment counts (one
    count per aligned read pair -- see run_sample.sh's counting logic)."""
    # Load every sample's {ref_id: count} dict, keyed by tc_id (dicts keep insertion
    # order, so columns end up in sample-sheet order).
    cols = {}
    for tc_id in samples["tc_id"]:
        cols[tc_id] = load_counts(tc_id)
    # Empty frame with the right row index / column labels. Note: with no data,
    # pandas ignores dtype="int64" and fills with NaN (float64); each column is
    # overwritten with ints just below, so the end result is int64 anyway.
    df = pd.DataFrame(index=meta["ref_id"], columns=list(cols.keys()), dtype="int64")
    # .items() yields (key, value) pairs, unpacked straight into two names.
    for tc_id, counts in cols.items():
        # For each peptide ID in the row index, look it up in this sample's dict,
        # defaulting to 0 if absent. Index.map applies the function element-wise.
        # Any ref_id in counts.tsv but not in the metadata is silently dropped here.
        df[tc_id] = df.index.map(lambda r: counts.get(r, 0)).astype("int64")
    return df


# parse_cutadapt_paired_log: pull read-pair totals out of a cutadapt report.
#   In:  path -- Path to results/<tc_id>/cutadapt.log.
#   Out: (total read pairs in, read pairs kept after trimming) as ints.
#   Used by sample_readstats() to compute raw_read_pairs and trim_rate.
def parse_cutadapt_paired_log(path):
    """Parse a paired-mode cutadapt summary (pair-level fields, not the
    single-end 'Total reads processed'/'Reads written' field names)."""
    text = path.read_text()  # whole log as one string
    # Matches e.g. "Total read pairs processed:         57,535,032".
    # \s* = any amount of whitespace; ([\d,]+) captures the digits-and-commas
    # number as group 1; the commas are then removed so int() can parse it.
    # (If the line is absent, re.search returns None and .group raises AttributeError.)
    total = int(re.search(r"Total read pairs processed:\s*([\d,]+)", text).group(1).replace(",", ""))
    # Matches e.g. "Pairs written (passing filters):    54,486,460 (94.7%)".
    # The literal parentheses must be escaped as \( \) in the regex; the
    # trailing "(94.7%)" is simply not captured.
    written = int(re.search(r"Pairs written \(passing filters\):\s*([\d,]+)", text).group(1).replace(",", ""))
    return total, written


# parse_bowtie2_paired_log: pull alignment counts out of bowtie2's stderr summary.
#   In:  path -- Path to results/<tc_id>/bowtie2.log.
#   Out: (total pairs, pairs aligned concordantly 0 / exactly 1 / >1 times,
#         overall alignment rate in %).
#   Used by sample_readstats() for the mapping-rate QC columns.
# Background: a pair aligns "concordantly" when both mates map to the same
# reference in the expected orientation and at a plausible distance. "Exactly
# 1 time" = unique hit; ">1 times" = multimapping (several equally good hits;
# bowtie2 reports one of them). The overall rate is per *mate* (individual
# read), so it also counts discordant pairs and single mates that aligned.
# Example log (first lines and last line):
#   54486460 reads; of these:
#     54486460 (100.00%) were paired; of these:
#       2178132 (4.00%) aligned concordantly 0 times
#       52282532 (95.96%) aligned concordantly exactly 1 time
#       25796 (0.05%) aligned concordantly >1 times
#   ...
#   97.14% overall alignment rate
def parse_bowtie2_paired_log(path):
    """Parse a paired-mode bowtie2 summary: concordant-pair alignment
    counts plus the overall (mate-level) alignment rate it reports on its
    last line."""
    text = path.read_text()
    # "54486460 reads; of these:" -- in paired mode bowtie2 says "reads" but means
    # pairs. (\d+) captures the integer; re.search returns the first match.
    total_pairs = int(re.search(r"(\d+) reads; of these:", text).group(1))
    # Each of the next three: capture the count before a "(NN.NN%)" percentage
    # ([\d.]+ = digits and dots), followed by the literal category text.
    concordant_0 = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly 0 times", text).group(1))
    concordant_1 = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", text).group(1))
    concordant_multi = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly >1 times", text).group(1))
    # "97.14% overall alignment rate" -> 97.14 (still a percentage, not a fraction).
    overall_pct = float(re.search(r"([\d.]+)% overall alignment rate", text).group(1))
    return total_pairs, concordant_0, concordant_1, concordant_multi, overall_pct


# sample_readstats: per-sample read/alignment QC numbers.
#   In:  tc_id -- sample ID.
#   Out: dict with raw_read_pairs, trim_rate, concordant_map_rate,
#        concordant_multimap_rate, overall_map_rate (all fractions 0-1 except the
#        first); NaN for anything whose log file is missing.
#   Used by main() to fill these columns of qc_summary.csv.
def sample_readstats(tc_id):
    """cutadapt/bowtie2 pair-level QC stats for one sample. One log file
    per tool per sample now (lanes already merged and R1/R2 already
    trimmed/aligned together upstream in run_sample.sh), so there's nothing
    left to sum across lanes/mates the way this used to."""
    # Start with NaN ("not a number" = missing) for everything, so a missing
    # log simply leaves its columns empty rather than crashing.
    stats = {
        "raw_read_pairs": float("nan"),
        "trim_rate": float("nan"),
        "concordant_map_rate": float("nan"),
        "concordant_multimap_rate": float("nan"),
        "overall_map_rate": float("nan"),
    }
    cpath = RESULTS_DIR / tc_id / "cutadapt.log"
    if cpath.exists():
        raw, written = parse_cutadapt_paired_log(cpath)  # unpack the returned 2-tuple
        stats["raw_read_pairs"] = raw  # pairs straight off the sequencer
        # Fraction of pairs that kept an adapter match in both mates (run_sample.sh
        # uses --discard-untrimmed). Conditional expression guards against /0.
        stats["trim_rate"] = written / raw if raw else float("nan")
    bpath = RESULTS_DIR / tc_id / "bowtie2.log"
    if bpath.exists():
        # Unpack all five values; the leading underscore marks concordant_0 as unused.
        total_pairs, _concordant_0, concordant_1, concordant_multi, overall_pct = parse_bowtie2_paired_log(bpath)
        # Fraction of (trimmed) pairs aligning concordantly at all = unique + multimapped.
        # Denominator is pairs given to bowtie2, i.e. after trimming, not raw pairs.
        stats["concordant_map_rate"] = (concordant_1 + concordant_multi) / total_pairs if total_pairs else float("nan")
        # Fraction of pairs that multimapped (ambiguous placement between peptides).
        stats["concordant_multimap_rate"] = concordant_multi / total_pairs if total_pairs else float("nan")
        stats["overall_map_rate"] = overall_pct / 100.0  # percent -> fraction
    return stats


# gini: Gini coefficient of a sample's per-peptide counts.
#   In:  counts -- 1-D array-like of counts (one sample, own-library peptides).
#   Out: float in [0, 1): 0 = every peptide has the same count (perfectly even
#        library), towards 1 = all reads on a single peptide. NaN if no reads.
#   Used by main() for the gini_own_library QC column.
# It is 1 - 2 x (area under the Lorenz curve), computed here with the closed
# form G = (n + 1 - 2 * sum_i(cum_i) / cum_n) / n over ascending-sorted values.
def gini(counts):
    x = np.sort(np.asarray(counts, dtype=float))  # to a float numpy array, sorted ascending
    x = x[x >= 0]  # boolean mask: keep non-negative values (counts always are; just a guard)
    n = len(x)
    if n == 0 or x.sum() == 0:  # undefined for an empty or all-zero sample
        return float("nan")
    cum = np.cumsum(x)  # running total: cum[i] = x[0] + ... + x[i]; cum[-1] = grand total
    return (n + 1 - 2 * np.sum(cum) / cum[-1]) / n


# lorenz_curve: points of the Lorenz curve for one sample.
#   In:  counts -- 1-D array-like of counts.
#   Out: (frac_pop, cum) -- x = cumulative fraction of peptides (0..1),
#        y = cumulative fraction of reads they hold (0..1), poorest peptides first.
#   Used by plot_lorenz(). A perfectly even library gives the diagonal y = x;
#   the further the curve sags below it, the more uneven the library.
def lorenz_curve(counts):
    x = np.sort(np.asarray(counts, dtype=float))  # ascending: least-read peptides first
    cum = np.cumsum(x)  # running total of reads
    # Prepend a 0 so the curve starts at the origin, then divide by the total
    # (the last element) to turn read totals into fractions. (All-zero sample -> 0/0 = NaN.)
    cum = np.insert(cum, 0, 0) / cum[-1]
    frac_pop = np.linspace(0, 1, len(cum))  # n+1 evenly spaced x-values from 0 to 1
    return frac_pop, cum


# shannon_entropy: Shannon entropy H = -sum(p * log(p)) of the read distribution.
#   In:  counts -- 1-D array-like of counts; base -- log base (2 -> bits).
#   Out: float; 0 = all reads on one peptide, maximum log2(n) = all n peptides equal.
#        NaN if the sample has no reads.
#   Used by main(), which also divides by log2(designed library size) to get a
#   normalised 0-1 score comparable between libraries of different sizes.
def shannon_entropy(counts, base=2):
    """Shannon diversity (in bits, base 2) of a sample's read counts across
    its own peptide library: how evenly are reads spread across peptides,
    vs. concentrated on a few?"""
    x = np.asarray(counts, dtype=float)
    total = x.sum()
    if total <= 0:
        return float("nan")
    # Proportions of reads per peptide; zero-count peptides are dropped because
    # p*log(p) -> 0 as p -> 0 (and log(0) would give -inf).
    p = x[x > 0] / total
    # Natural-log entropy, converted to the requested base via log_b(y) = ln(y)/ln(b).
    return float(-np.sum(p * np.log(p)) / np.log(base))


# compute_pca: principal component analysis of one figure group's samples.
#   In:  mat -- full raw peptide x sample count matrix; own_refs -- ref_ids of the
#        (single) library's peptides; tc_ids -- this group's samples (in plot order).
#   Out: (scores [samples x PCs array], explained_var [fraction per PC],
#        n_peptides_eligible), or None if too few samples (<3) / peptides (<2).
#   Used by main() for each figure group; output goes to plot_pca/write_pca_scores.
def compute_pca(mat, own_refs, tc_ids):
    """PCA of just these samples (one figure group). Peptides are the
    group's own library peptides passing the matrix filter *within the
    group* (count_filters.passes_filter on these samples; zeros allowed),
    in log10(CPM + 1) -- CPM scaled by
    each sample's raw own-library read total, since depth varies several-
    fold -- then the top N most variable, mean-centered per peptide.
    Returns (scores, explained_var, n_peptides_eligible), or None if there
    isn't enough signal to run PCA on."""
    # .loc[row labels, column labels]: this library's peptides x this group's samples.
    sub = mat.loc[own_refs, tc_ids]
    # CPM (counts per million): divide each column by its column total
    # (sum(axis=0) = sum down the rows, one value per sample; pandas broadcasts the
    # division across columns by label) and scale to 1e6, removing depth differences.
    # Totals are taken BEFORE filtering, so CPM reflects each sample's real depth.
    cpm = sub / sub.sum(axis=0) * 1e6
    # Boolean Series per peptide: >= 10 reads in >= 2 of these samples.
    keep = passes_filter(sub)
    # PCA needs at least 2 features and 3 samples to give anything meaningful.
    # keep.sum() counts True values.
    if keep.sum() < 2 or len(tc_ids) < 3:
        return None
    # log10(x + 1) compresses the huge dynamic range of counts (and maps 0 -> 0),
    # so a few very abundant peptides don't dominate the variance.
    log_cpm = np.log10(cpm.loc[keep] + 1)
    # Method chain: per-peptide variance across samples (axis=1 = along each row),
    # sort largest first, take the top N, keep just their ref_id labels (.index).
    top_ids = log_cpm.var(axis=1).sort_values(ascending=False).head(PCA_N_TOP_VARIABLE_PEPTIDES).index

    # Transpose (.T) so rows = samples (the observations) and columns = peptides
    # (the features), then convert to a plain numpy array.
    X = log_cpm.loc[top_ids].T.to_numpy()  # samples x peptides
    # Mean-centre each peptide column (axis=0 = mean down the rows); keepdims keeps
    # shape (1, n_peptides) so it broadcasts against X row by row.
    X = X - X.mean(axis=0, keepdims=True)
    # PCA via singular value decomposition X = U * diag(S) * Vt. full_matrices=False
    # gives the compact form (only min(n_samples, n_peptides) components).
    # Vt (peptide loadings) is not needed, hence the _ prefix.
    U, S, _Vt = np.linalg.svd(X, full_matrices=False)
    # Sample coordinates on each PC = U scaled by the singular values
    # (broadcasting S across U's columns), equivalent to X @ Vt.T.
    scores = U * S
    # Variance along each PC is proportional to S^2; divide by the total to get
    # each PC's fraction of the overall variance.
    explained_var = (S ** 2) / np.sum(S ** 2)
    # int(...) converts numpy's integer to a plain Python int.
    return scores, explained_var, int(keep.sum())

# own_library_counts: one sample's counts restricted to its own library.
#   In:  mat -- count matrix; lib_of_ref -- Series ref_id -> library name;
#        srow -- one row of the sample sheet (needs "library" and "tc_id").
#   Out: Series of counts indexed by ref_id, own-library peptides only.
#   Used by the rank-abundance and Lorenz plots (each sample is only judged on the
#   library it was made from; cross-library hits would be contamination).
def own_library_counts(mat, lib_of_ref, srow):
    """One sample's counts over its own library's peptides."""
    # Boolean mask on the Series selects peptides whose library matches the
    # sample's; .index gives their ref_ids.
    own_refs = lib_of_ref[lib_of_ref == srow["library"]].index
    # Pick the sample's column, then its own-library rows by label.
    return mat[srow["tc_id"]].loc[own_refs]


# line_legend: legend for the line plots (rank-abundance, Lorenz), placed outside
#   the axes on the right: one block of treatment colours, then one block
#   explaining the line styles (ip in the all-samples figure, replicate per IP).
#   In: ax, figure group dict, treatment -> colour dict, style value -> linestyle dict.
def line_legend(ax, group, colors, linestyles):
    # The two handle lists are concatenated with + into one legend.
    legend_outside(ax, treatment_legend_handles(group["samples"], colors, kind="line")
                   + style_legend_handles(group, linestyles, kind="line"), fontsize=8)


# plot_rank_abundance: log-log rank-abundance curve per sample.
#   In:  mat, lib_of_ref, group (dict from figure_groups), colors (treatment -> colour).
#   Out: writes <group fig_dir>/rank_abundance.png. Called once per figure group.
# Peptides are ranked by count (rank 1 = most abundant). A flat curve = even
# library; a steep drop = a few peptides dominate. Zero counts are left out
# (log scale cannot show 0), so where each line ends = number of peptides detected.
def plot_rank_abundance(mat, lib_of_ref, group, colors):
    """Rank-abundance (log-log) per sample, restricted to each sample's own library."""
    # Map each value of the group's style column (ip or replicate) to a line style.
    linestyles = style_map(group, LINESTYLE_CYCLE)
    fig, ax = plt.subplots(figsize=(7, 5))  # one figure with one axes; size in inches
    # iterrows() yields (index, row Series) pairs; the index isn't needed (_).
    for _, srow in group["samples"].iterrows():
        # .values -> numpy array; np.sort is ascending, [::-1] reverses to descending.
        vals = np.sort(own_library_counts(mat, lib_of_ref, srow).values)[::-1]
        vals = vals[vals > 0]  # drop undetected peptides (log axis)
        # x = ranks 1..n, y = counts; colour by treatment, line style by ip/replicate.
        ax.plot(np.arange(1, len(vals) + 1), vals, color=colors[srow["treatment"]],
                linestyle=linestyles[srow[group["style_col"]]], linewidth=1.2)
    ax.set_xscale("log")  # log-log axes: counts span several orders of magnitude
    ax.set_yscale("log")
    ax.set_xlabel("Peptide rank")
    ax.set_ylabel("Read-pair count")
    ax.set_title("Rank-abundance (own-library peptides)" + group["title_suffix"])
    line_legend(ax, group, colors, linestyles)
    # bbox_inches="tight" grows the saved image so the outside legend isn't clipped.
    fig.savefig(group["fig_dir"] / "rank_abundance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)  # free the figure's memory (many figures are made in a loop)


# plot_lorenz: Lorenz curve per sample (library evenness), own-library peptides.
#   In:  as plot_rank_abundance. Out: writes <group fig_dir>/lorenz_curves.png.
#   The diagonal = perfect evenness; area between it and a curve is half the Gini.
def plot_lorenz(mat, lib_of_ref, group, colors):
    linestyles = style_map(group, LINESTYLE_CYCLE)  # style value -> line style
    fig, ax = plt.subplots(figsize=(6, 6))  # square, since both axes run 0..1
    for _, srow in group["samples"].iterrows():
        # Tuple-unpack the (x, y) arrays returned by lorenz_curve.
        frac_pop, frac_reads = lorenz_curve(own_library_counts(mat, lib_of_ref, srow).values)
        ax.plot(frac_pop, frac_reads, color=colors[srow["treatment"]],
                linestyle=linestyles[srow[group["style_col"]]], linewidth=1.2)
    # Reference diagonal: "k--" is matplotlib shorthand for black (k) dashed (--).
    ax.plot([0, 1], [0, 1], "k--", linewidth=1)
    # Label placed along the diagonal (data coordinates, rotated 45 degrees).
    ax.text(0.55, 0.5, "perfect evenness", rotation=45, fontsize=8, ha="center", va="bottom")
    ax.set_xlabel("Cumulative fraction of peptides")
    ax.set_ylabel("Cumulative fraction of reads")
    ax.set_title("Lorenz curves (library evenness)" + group["title_suffix"])
    line_legend(ax, group, colors, linestyles)
    fig.savefig(group["fig_dir"] / "lorenz_curves.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_batch_correlation: sample-vs-sample correlation heatmap, one per library.
#   In:  mat, lib_of_ref, group, colors.
#   Out: writes <group fig_dir>/batch_correlation_<library>.png for each library
#        present in the group. Called once per figure group.
# High correlation between replicates = reproducible IPs; low = noisy/batch effects.
def plot_batch_correlation(mat, lib_of_ref, group, colors):
    """Pearson correlation heatmap of log1p counts, per library, over the
    peptides passing the matrix filter within this group's samples, with
    treatment colour strips along the top and left edges."""
    # groupby yields (library name, sub-DataFrame of that library's samples);
    # sort=False keeps the libraries in order of first appearance.
    for library, sub_samples in group["samples"].groupby("library", sort=False):
        # This library's peptide IDs, intersected with the matrix rows to be safe
        # (they are the same set, since the matrix is built from the metadata).
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        tc_ids = sub_samples["tc_id"].tolist()
        # List comprehension: one tick label per sample (ip-prefixed in the "all" group).
        labels = [sample_label(srow, group) for _, srow in sub_samples.iterrows()]
        n = len(tc_ids)
        sub = mat.loc[own_refs, tc_ids]  # peptides x these samples
        # Keep only peptides passing the matrix filter within these samples
        # (boolean Series used as a row mask in .loc).
        sub = sub.loc[passes_filter(sub)]
        # log1p(x) = ln(1 + x): tames huge counts, keeps 0 at 0. DataFrame.corr
        # computes the Pearson correlation between every pair of COLUMNS (samples),
        # giving an n x n matrix.
        corr = np.log1p(sub).corr(method="pearson")

        # Figure grows with the number of samples so the cells stay readable.
        fig, ax = plt.subplots(figsize=(1.2 * n + 2, 1.2 * n + 2))
        # Draw the matrix as an image; colour scale fixed to 0..1 (negative
        # correlations would be clipped to the bottom colour).
        im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="viridis")
        # Write each correlation value in its cell. imshow puts column j at x=j and
        # row i at y=i. f"{...:.2f}" formats to 2 decimal places.
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)

        # Attach thin extra axes to the heatmap's edges: colour strips on top and
        # left (sharing the heatmap's x / y axis so cells line up exactly), and a
        # colorbar on the right. size/pad are in inches.
        divider = make_axes_locatable(ax)
        ax_top = divider.append_axes("top", size=0.3, pad=0.05, sharex=ax)
        ax_left = divider.append_axes("left", size=0.3, pad=0.05, sharey=ax)
        cax = divider.append_axes("right", size=0.3, pad=0.15)
        # One RGBA colour per sample (by treatment) -> array of shape (n, 4).
        strip = np.array([mcolors.to_rgba(colors[t]) for t in sub_samples["treatment"]])
        # imshow of an RGBA image: np.newaxis adds a dimension, making the strip
        # 1 x n x 4 (a row, for the top) or n x 1 x 4 (a column, for the left).
        # aspect="auto" lets the cells stretch to fill the thin axes.
        ax_top.imshow(strip[np.newaxis, :, :], aspect="auto")
        ax_left.imshow(strip[:, np.newaxis, :], aspect="auto")

        # Sample labels along the bottom of the heatmap, angled for readability
        # (ha="right" anchors each rotated label's end at its tick).
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=HEATMAP_TICK_FONTSIZE)
        ax.set_yticks(range(n))
        # Hide the heatmap's own y labels; they're drawn on the left strip instead,
        # so they sit outside the strip rather than between strip and heatmap.
        ax.tick_params(axis="y", left=False, labelleft=False)
        # ax_left shares the y axis, so it already has the n tick positions.
        ax_left.set_yticklabels(labels, fontsize=HEATMAP_TICK_FONTSIZE)
        ax_left.set_xticks([])  # no ticks along the strips' short sides
        ax_top.set_yticks([])
        # The top strip shares the x axis; hide its bottom tick labels so the sample
        # names appear only once (under the heatmap).
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        # Title goes on the top strip (the topmost axes), so it sits above everything.
        # Two adjacent f-strings are joined automatically; \n starts a second line.
        ax_top.set_title(f"{library} batch correlation{group['title_suffix']}\n"
                         f"(log1p read-pair counts; {len(sub)} peptides with {GROUP_FILTER_DESCRIPTION})")
        fig.colorbar(im, cax=cax)  # colour scale drawn into the right-hand axes
        # anchored to the colorbar (in its axes units) so it sits just past the colorbar's tick labels
        cax.legend(handles=treatment_legend_handles(sub_samples, colors), loc="center left",
                   bbox_to_anchor=(3.5, 0.5), borderaxespad=0, fontsize=HEATMAP_TICK_FONTSIZE)
        fig.savefig(group["fig_dir"] / f"batch_correlation_{library}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


# plot_shannon: bar chart of normalised Shannon entropy per sample.
#   In:  qc -- the qc_summary DataFrame built in main(); group; colors.
#   Out: writes <group fig_dir>/shannon_entropy.png. Called once per figure group.
#   1.0 = reads spread perfectly evenly over the whole designed library.
def plot_shannon(qc, group, colors):
    # Index the QC table by sample ID, then .loc with the group's tc_ids both
    # subsets it to this group and puts the rows in the group's plotting order.
    group_qc = qc.set_index("tc_id").loc[group["samples"]["tc_id"]]
    labels = [sample_label(srow, group) for _, srow in group["samples"].iterrows()]
    # Width scales with the number of bars.
    fig, ax = plt.subplots(figsize=(1.2 * len(labels) + 2, 4.5))
    # Bars at x = 0..n-1; list comprehension gives each bar its treatment colour.
    ax.bar(range(len(labels)), group_qc["normalized_shannon_entropy"],
           color=[colors[t] for t in group["samples"]["treatment"]], edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=XTICK_ROTATION, ha="center")
    ax.set_ylim(0, 1.05)  # the normalised score is 0..1; small headroom above
    ax.set_ylabel("Normalized Shannon entropy\n(H / log2(designed library size))")
    ax.set_title("Library diversity: how evenly are reads spread across peptides?" + group["title_suffix"])
    legend_outside(ax, treatment_legend_handles(group["samples"], colors), fontsize=8)
    fig.savefig(group["fig_dir"] / "shannon_entropy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# write_pca_scores: save PCA coordinates for one figure group.
#   In:  scores, explained_var (from compute_pca), group, n_pcs (max PCs to write).
#   Out: writes <group fig_dir>/pca_scores.csv with sample info + PC columns
#        named like "PC1 (34.2%)". Called by main() after a successful PCA.
def write_pca_scores(scores, explained_var, group, n_pcs=5):
    """Per-sample PC scores (first n_pcs) for this group, with each PC's
    variance explained in the column name."""
    samples = group["samples"]
    # Can't write more PCs than exist (scores has min(n_samples, n_peptides) columns).
    n_pcs = min(n_pcs, scores.shape[1])
    # Double brackets select a list of columns -> DataFrame; .copy() so adding
    # columns below doesn't modify (or warn about) the original sample table.
    out = samples[["tc_id", "ip", "treatment", "replicate"]].copy()
    for k in range(n_pcs):
        # Format spec ".1%" multiplies by 100 and adds a % sign, 1 decimal place.
        # scores[:, k] = column k (all samples' scores on PC k+1); rows match samples' order.
        out[f"PC{k + 1} ({explained_var[k]:.1%})"] = scores[:, k]
    out.to_csv(group["fig_dir"] / "pca_scores.csv", index=False)  # no row-number column


# plot_pca: scatter of PC1 vs PC2 for one figure group.
#   In:  scores, explained_var, n_eligible (from compute_pca), group, colors.
#   Out: writes <group fig_dir>/pca_samples.png. Colour = treatment,
#        marker shape = ip (all-samples figure) or replicate (per-IP figures).
def plot_pca(scores, explained_var, n_eligible, group, colors):
    """PCA computed on this group's samples only (see compute_pca)."""
    samples = group["samples"]
    markers = style_map(group, MARKER_CYCLE)  # style value -> marker shape
    fig, ax = plt.subplots(figsize=(7.5, 6))
    # One scatter call per marker shape (matplotlib can't vary marker within a call).
    for value, marker in markers.items():
        # np.where on a boolean array returns a tuple of index arrays; [0] takes the
        # positions (row numbers) of samples whose style column equals this value.
        sel = np.where(samples[group["style_col"]].values == value)[0]
        # scores[sel, 0] / scores[sel, 1] = those samples' PC1 / PC2. Per-point colours
        # come from a list comprehension; .iloc looks the treatment up by position,
        # which matches the row order of scores.
        ax.scatter(scores[sel, 0], scores[sel, 1], c=[colors[samples["treatment"].iloc[i]] for i in sel],
                   marker=marker, s=70, edgecolor="black", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({explained_var[0]:.1%} variance)")
    # Conditional expression: only show PC2's % if a second PC exists.
    ax.set_ylabel(f"PC2 ({explained_var[1]:.1%} variance)" if len(explained_var) > 1 else "PC2")
    which = "all samples" if group["name"] == "all" else f"{group['name']} samples only"
    # Three-line title; min(...) because fewer than N peptides may be eligible.
    ax.set_title(f"PCA of samples{group['title_suffix']}\n(computed on {which}: log10(CPM + 1), top "
                 f"{min(PCA_N_TOP_VARIABLE_PEPTIDES, n_eligible)} most variable of {n_eligible} peptides\n"
                 f"with {GROUP_FILTER_DESCRIPTION})",
                 fontsize=10)
    # kind="marker" makes legend entries dots/shapes rather than lines.
    legend_outside(ax, treatment_legend_handles(group["samples"], colors, kind="marker")
                   + style_legend_handles(group, markers, kind="marker"), fontsize=8)
    fig.savefig(group["fig_dir"] / "pca_samples.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# main: the whole run, top to bottom.
#   1. build + save the raw count matrix; 2. apply/save the all-sample matrix filter;
#   3. compute per-sample QC -> qc_summary.csv; 4. draw every figure for each
#   figure group (+ PCA if all samples share one library); 5. write summary.md.
def main():
    # Create results/figures/ (and parents) if needed; exist_ok avoids an error if present.
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    # Merged samples.tsv + samples_meta.csv (CLI args 1/2 override the paths),
    # sorted ip -> treatment -> sheet order, plus the treatment -> colour dict.
    samples, colors = load_samples()
    # dtype=str reads every column as text (stops pandas guessing types).
    meta = pd.read_csv(META_CSV, dtype=str)
    meta["ref_id"] = meta["ref_id"].astype(str)  # already str via dtype=str; harmless safeguard

    print("Building count matrix...", file=sys.stderr)  # progress messages go to stderr
    mat = build_matrix(samples, meta)
    mat.to_csv(RESULTS_DIR / "count_matrix.csv")  # row index (ref_id) written as first column

    # Series mapping ref_id -> library ("CoV"/"Vir3"): set_index makes ref_id the
    # index, then ["library"] picks that column.
    lib_of_ref = meta.set_index("ref_id")["library"]

    # The all-sample filtered matrix is written for sharing (and feeds
    # top_peptides.py's tables); the between-sample figures apply the same
    # filter within each figure group instead (passes_filter); library-
    # completeness metrics below use the raw matrix.
    # filter_matrix returns 3 things: the kept rows, a per-sample table of reads
    # kept/removed, and a per-library table of peptides kept/removed.
    mat_filtered, filter_report, filter_library_report = filter_matrix(mat, lib_of_ref)
    mat_filtered.to_csv(RESULTS_DIR / "count_matrix_filtered.csv")
    filter_report.to_csv(RESULTS_DIR / "filter_report.csv", index=False)
    print(f"Matrix filter ({FILTER_DESCRIPTION}):\n{filter_library_report.to_string()}", file=sys.stderr)

    # --- Per-sample QC: build a list of dicts (one per sample), then one DataFrame ---
    qc_rows = []
    for _, srow in samples.iterrows():
        tc_id, library = srow["tc_id"], srow["library"]  # tuple unpacking of two values
        # Peptides of the sample's own library vs. those of the other library.
        own_refs = lib_of_ref[lib_of_ref == library].index
        other_refs = lib_of_ref[lib_of_ref != library].index

        col = mat[tc_id]  # this sample's counts (Series indexed by ref_id)
        own_counts = col.loc[own_refs]
        other_counts = col.loc[other_refs]

        stats = sample_readstats(tc_id)  # values parsed from the cutadapt/bowtie2 logs

        # All counted read pairs (own + other library = every ref in the matrix).
        total_mapped = own_counts.sum() + other_counts.sum()
        n_designed = len(own_refs)  # designed library size = denominator for coverage
        entropy = shannon_entropy(own_counts.values)
        qc_rows.append({
            "tc_id": tc_id,
            "library": library,
            "ip": srow["ip"],
            "treatment": srow["treatment"],
            "raw_read_pairs": stats["raw_read_pairs"],
            "trim_rate": stats["trim_rate"],
            "concordant_map_rate": stats["concordant_map_rate"],
            "concordant_multimap_rate": stats["concordant_multimap_rate"],
            "overall_map_rate": stats["overall_map_rate"],
            # Cross-library mapping: fraction of counted pairs landing on the OTHER
            # library's peptides -- a contamination / sample-mix-up check. (Stored as a
            # fraction 0-1 despite the "pct" name.)
            "pct_mapped_to_other_library": other_counts.sum() / total_mapped if total_mapped else float("nan"),
            "n_designed_peptides": n_designed,
            # (own_counts > 0) is a boolean Series; .sum() counts the Trues and
            # .mean() gives the fraction True (here: fraction of library detected).
            "n_peptides_detected": int((own_counts > 0).sum()),
            "pct_library_covered": (own_counts > 0).mean(),  # fraction, not %
            "n_dropouts": int((own_counts == 0).sum()),  # designed peptides never seen
            # Same coverage metrics at the stricter >= coverage_min_count (10) reads
            # threshold; the f-string builds the column name, e.g. "..._min10".
            f"n_peptides_detected_min{COVERAGE_MIN_COUNT}": int((own_counts >= COVERAGE_MIN_COUNT).sum()),
            f"pct_library_covered_min{COVERAGE_MIN_COUNT}": (own_counts >= COVERAGE_MIN_COUNT).mean(),
            # Typical count among detected peptides (boolean mask drops zeros first).
            "median_count_detected": own_counts[own_counts > 0].median(),
            "gini_own_library": gini(own_counts.values),
            "shannon_entropy": entropy,  # in bits
            # Normalise by the maximum possible entropy, log2(n_designed) (all designed
            # peptides equally abundant), giving 0..1. Undefined for a 1-peptide library.
            "normalized_shannon_entropy": entropy / np.log2(n_designed) if n_designed > 1 else float("nan"),
        })

    # List of dicts -> DataFrame: dict keys become columns, one row per dict.
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(RESULTS_DIR / "qc_summary.csv", index=False)
    print(qc.to_string(index=False), file=sys.stderr)  # full table to the job log

    # --- Figures: all samples, then one set per IP condition ---
    # List of group dicts (name, samples, fig_dir, title_suffix, style_col);
    # also creates the per-IP figure folders.
    groups = figure_groups(samples, FIG_DIR)
    libraries = samples["library"].unique()  # numpy array of distinct library names
    # PCA only makes sense when all samples share one library (peptide sets differ
    # between libraries), so it is skipped entirely for mixed-library runs.
    pca_written = len(libraries) == 1
    if not pca_written:
        print(f"Skipping PCA: samples span several libraries ({', '.join(libraries)}).", file=sys.stderr)

    for group in groups:
        plot_rank_abundance(mat, lib_of_ref, group, colors)
        plot_lorenz(mat, lib_of_ref, group, colors)
        plot_batch_correlation(mat, lib_of_ref, group, colors)
        plot_shannon(qc, group, colors)
        if pca_written:
            # Peptides of the single library present (libraries[0]).
            own_refs = lib_of_ref[lib_of_ref == libraries[0]].index.intersection(mat.index)
            pca = compute_pca(mat, own_refs, group["samples"]["tc_id"].tolist())
            if pca is None:
                print(f"Skipping PCA{group['title_suffix']}: too few samples or peptides.", file=sys.stderr)
            else:
                scores, explained_var, n_eligible = pca  # unpack the 3-tuple
                plot_pca(scores, explained_var, n_eligible, group, colors)
                write_pca_scores(scores, explained_var, group)
                print(f"PCA{group['title_suffix']}: {n_eligible} eligible peptides, PC1 {explained_var[0]:.1%}, "
                      f"PC2 {explained_var[1]:.1%}", file=sys.stderr)

    # --- Markdown summary ---
    # Build the report as a list of lines, joined with newlines at the end.
    lines = ["# PhIP-seq library QC summary", ""]
    lines.append("## Per-sample QC")
    lines.append("")
    # DataFrame -> Markdown table (needs the `tabulate` package); floats to 3 d.p.
    lines.append(qc.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    # Links to files written by later scripts (top_peptides.py, differential_enrichment.py).
    lines.append("Top most abundant peptides per sample: [top_abundant_peptides.csv](top_abundant_peptides.csv)")
    lines.append("")
    lines.append("Differential enrichment vs the background library (written separately by "
                 "differential_enrichment.py): [differential/summary.md](differential/summary.md)")
    lines.append("")
    # Plain-language explanation of the low-count filters, with the configured
    # thresholds filled in by the f-strings.
    lines.append("## Low-count filtering")
    lines.append("")
    lines.append(f"- Coverage: `*_min{COVERAGE_MIN_COUNT}` QC columns and the solid coverage bars count a peptide as "
                 f"detected only with >= {COVERAGE_MIN_COUNT} reads; the denominator is always the full designed "
                 f"library. Rank-abundance, Lorenz, Shannon entropy and evenness use the raw matrix.")
    lines.append(f"- Matrix filter (>= {MATRIX_FILTER_MIN_COUNT} reads in >= {MATRIX_FILTER_MIN_SAMPLES} samples), "
                 f"applied within each figure group (all samples, each IP) on that group's samples only for the PCA, "
                 f"batch correlation and top-variable heatmaps -- so the IP4 figures only use peptides seen in IP4 "
                 f"samples. Across all samples it gives [count_matrix_filtered.csv](count_matrix_filtered.csv) "
                 f"(also used for the top-abundant tables; per-sample reads removed in "
                 f"[filter_report.csv](filter_report.csv)), summarised below.")
    lines.append(f"- PCA: computed separately for each figure group on log10(CPM + 1) of that group's filtered "
                 f"peptides (top {PCA_N_TOP_VARIABLE_PEPTIDES} most variable); scores in each figures folder's "
                 f"pca_scores.csv. Axes are not comparable between groups.")
    lines.append("")
    # Per-library peptides raw/kept/removed table; reset_index turns the "library"
    # index into an ordinary column so it shows up in the Markdown table.
    lines.append(filter_library_report.reset_index().to_markdown(index=False))
    lines.append("")
    # Names (without .png) of every figure to embed, in display order. Includes
    # figures made by the other scripts run after this one (completeness_figures.py,
    # evenness_summary_figure.py, overview_figure.py, top_peptides.py).
    figure_names = ["rank_abundance", "lorenz_curves"]
    # One batch-correlation figure per library (list comprehension over unique libraries).
    figure_names += [f"batch_correlation_{lib}" for lib in samples["library"].unique()]
    # Add the PCA figure only when PCA ran (list + empty-or-one-item list).
    figure_names += ["shannon_entropy"] + (["pca_samples"] if pca_written else [])
    figure_names += ["coverage_bar", "rarefaction_curves", f"rarefaction_curves_min{COVERAGE_MIN_COUNT}",
                     "detection_histogram", "evenness_summary",
                     "overview_all_samples", "top_variable_heatmap"]
    # One "## Figures (...)" section per figure group, embedding each image.
    for group in groups:
        lines.append(f"## Figures{group['title_suffix']}")
        lines.append("")
        # Folder path relative to results/ (where summary.md lives), e.g. "figures/IP4",
        # with forward slashes (as_posix) so Markdown links work on any OS.
        rel_dir = group["fig_dir"].relative_to(RESULTS_DIR).as_posix()
        for name in figure_names:
            lines.append(f"![{name}]({rel_dir}/{name}.png)")  # Markdown image syntax
            lines.append("")
    (RESULTS_DIR / "summary.md").write_text("\n".join(lines))
    print(f"Wrote {RESULTS_DIR / 'summary.md'}", file=sys.stderr)


# Run main() only when executed as a script, not when imported as a module.
if __name__ == "__main__":
    main()
