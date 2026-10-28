"""One figure summarizing library evenness per the plain-language stats
already discussed: for each sample, the P10-P90 "typical range" of
per-peptide read counts, the median, and any near-dropout outliers (<=10
reads) called out explicitly rather than buried in a max/min ratio.
Coloured by treatment and drawn once for all samples (results/figures/)
and once per IP condition (results/figures/<ip>/, see sample_groups.py).

Usage: evenness_summary_figure.py [samples_tsv] [samples_meta_csv]
  (default to reference/samples.tsv and reference/samples_meta.csv)
"""
import json  # reads the JSON config files under reference/
import sys  # sys.stderr for progress messages (CLI args are read in sample_groups.py)
from pathlib import Path  # object-oriented file paths ("/" joins path parts)

import numpy as np  # used for percentiles and building x-position arrays
import pandas as pd  # tables: count matrix and peptide metadata
import matplotlib

# Select the non-interactive "Agg" backend (renders straight to PNG files, no
# window/display needed -- important on HPC nodes). Must happen before pyplot
# is imported, which is why this import order looks unusual.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # the plotting interface (figures, axes)
# Line2D lets us build "fake" line/marker artists purely to use as legend entries.
from matplotlib.lines import Line2D

# Shared helpers (scripts/sample_groups.py): which sample subsets to draw, the
# outside-the-axes legend, tick labels, the treatment colour legend, and sample loading.
from sample_groups import figure_groups, legend_outside, load_samples, sample_label, treatment_legend_handles

# Repo checkout root: this file is scripts/<name>.py, so go up two levels.
PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())  # parse JSON -> dict
RUN_DIR = Path(_paths["run_dir"])  # where pipeline outputs are read from/written to
RESULTS_DIR = RUN_DIR / "results"  # holds count_matrix.csv (written by summarize.py)
FIG_DIR = RESULTS_DIR / "figures"  # top-level figure folder (per-IP subfolders go inside)
# Per-peptide metadata (one row per designed peptide: ref_id, library, ...)
# written by build_reference.py; used to know which library each peptide belongs to.
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

# Figure thresholds live in reference/figure_params.json (shared with
# overview_figure.py) rather than being duplicated as a literal here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
OUTLIER_THRESHOLD = _fig_params["outlier_threshold"]  # reads; peptides at or below this are called out individually
XTICK_ROTATION = _fig_params["xtick_rotation"]  # degrees to rotate x tick labels (90 = vertical)

# Near-dropout markers are black so they stay distinct from the treatment palette.
OUTLIER_COLOR = "black"

# Room below the axes for the rotated "<ip> <tc_id> (N ND, N miss)" tick labels.
TICK_LABEL_HEIGHT_IN = 3.2


# plot_evenness: draw and save one evenness figure for one sample group.
#   In:  mat        -- raw count matrix (rows = peptide ref_ids, columns = sample tc_ids)
#        lib_of_ref -- Series mapping ref_id -> library name ("CoV"/"Vir3")
#        group      -- dict from figure_groups(): samples (DataFrame), fig_dir, title_suffix, name
#        colors     -- dict treatment -> matplotlib colour
#   Out: nothing returned; writes <group fig_dir>/evenness_summary.png.
#   Per sample (one x position each) it draws the P10-P90 "typical range" bar, the
#   median dot and every near-dropout peptide (1..OUTLIER_THRESHOLD reads) as a marker.
def plot_evenness(mat, lib_of_ref, group, colors):
    samples = group["samples"]  # the sample-sheet rows in this group, already in plotting order
    # One figure + one axes; width grows with the number of samples (1.4 in each + 2 in
    # for margins), height = 5.5 in of plot plus the space reserved for tick labels.
    fig, ax = plt.subplots(figsize=(1.4 * len(samples) + 2, 5.5 + TICK_LABEL_HEIGHT_IN))
    xtick_labels = []  # filled with one label per sample, applied after the loop

    # iterrows() yields (index, row) pairs; enumerate adds i = 0, 1, 2... used as the x position.
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]  # sample ID and the library it was grown from
        # "Own-library" peptides: boolean mask picks ref_ids whose library matches this
        # sample's, .index gives those ref_ids, and intersection() keeps only ones that are
        # actually rows of the count matrix (so the row lookup in .loc below cannot raise a KeyError).
        # Peptides from the other library are ignored -- they are not expected in this sample.
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]  # Series: read count of each own-library peptide
        detected = counts[counts > 0]  # boolean mask: peptides seen at least once
        # Dropouts: designed peptides with zero reads. (counts == 0) is a boolean Series;
        # .sum() counts the Trues; int() turns the numpy integer into a plain int.
        n_missing = int((counts == 0).sum())
        # Near-dropouts: detected (>0) but only at a very low count (<= threshold reads).
        outliers = detected[detected <= OUTLIER_THRESHOLD]
        # Only draw if something was detected (percentiles of an empty array would fail).
        if len(detected):
            # 10th, 50th (median) and 90th percentiles of the detected counts in one call.
            # P10-P90 = the range the middle 80% of peptides fall in: a robust evenness
            # measure that, unlike min/max, is not dominated by a few extreme peptides.
            p10, p50, p90 = np.percentile(detected, [10, 50, 90])

            # the tight "typical" band most peptides fall in
            # Drawn as a thick vertical line from (i, p10) to (i, p90) in the treatment colour;
            # "butt" caps make the bar end exactly at P10/P90; alpha = slight transparency.
            ax.plot([i, i], [p10, p90], color=colors[srow["treatment"]], linewidth=6, solid_capstyle="butt",
                    alpha=0.85)
            # Median as a white circle with a black edge; zorder=5 draws it on top of the bar.
            ax.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=8, zorder=5)

            # near-dropout outliers, called out individually rather than folded into a min/max ratio
            # np.full(n, i) = an array of n copies of i, so every outlier sits at this sample's x;
            # "v" = downward triangle marker, linestyle="" = markers only (no connecting line).
            ax.plot(np.full(len(outliers), i), outliers.values, "v", color=OUTLIER_COLOR, markersize=9,
                    linestyle="", zorder=6)
        # Build the tick-label suffix, e.g. " (3 ND, 12 miss)", from whichever counts are non-zero.
        label_bits = []
        if len(outliers):  # any near-dropouts?
            label_bits.append(f"{len(outliers)} ND")
        if n_missing:  # any peptides with zero reads?
            label_bits.append(f"{n_missing} miss")
        # Join the bits inside brackets; a sample with neither is labelled " (clean)".
        extra = " (" + ", ".join(label_bits) + ")" if label_bits else " (clean)"
        # Kept as a single line: rotated multi-line tick labels lay out very
        # differently from rotated single-line ones (each line becomes its
        # own column) and need much more room than is practical here.
        # sample_label() gives "<ip> <tc_id>" in the all-samples figure, just tc_id per IP.
        xtick_labels.append(f"{sample_label(srow, group)}{extra}")

    # Log y-axis: read counts span several orders of magnitude (a few reads to thousands).
    ax.set_yscale("log")
    # Half a slot of padding either side of the first (x=0) and last (x=n-1) sample.
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_xticks(range(len(samples)))  # one tick per sample at x = 0, 1, 2, ...
    # Replace the numeric ticks with the sample labels, rotated so long labels don't overlap.
    ax.set_xticklabels(xtick_labels, rotation=XTICK_ROTATION, ha="center", fontsize=8)
    ax.set_ylabel("Reads per peptide (log scale)")
    # Two-line title: question + group suffix (e.g. " (IP4)"), then the ND definition.
    ax.set_title("Library evenness: how tightly clustered is representation?" + group["title_suffix"]
                 + f"\n(ND = near-dropout, <={OUTLIER_THRESHOLD} reads)")
    # Legend entries for the plot elements, built from empty Line2D "proxy" artists that are
    # never drawn on the axes -- they only carry the style + label shown in the legend.
    handles = [Line2D([], [], color="gray", linewidth=6, label="typical range (10th-90th percentile)"),
               Line2D([], [], marker="o", linestyle="", color="white", markeredgecolor="black", markersize=8,
                      label="median"),
               Line2D([], [], marker="v", linestyle="", color=OUTLIER_COLOR, markersize=9,
                      label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)")]
    # One legend to the right of the axes (bbox_to_anchor in sample_groups.py): the plot-element
    # entries above plus a "treatment" header and one coloured square per treatment present.
    legend_outside(ax, handles + treatment_legend_handles(samples, colors), fontsize=8)
    # Fixed widths/heights (in inches, converted to figure fractions) are
    # reserved for the outside legend and the rotated tick labels. Not using
    # bbox_inches="tight" here: with this matplotlib version it can blow up
    # the canvas to a pathological size when combined with a log-scaled axis
    # and many plotted markers.
    # subplots_adjust takes fractions of the figure: the axes end 3.2 in from the right edge
    # (legend space) and start TICK_LABEL_HEIGHT_IN in above the bottom (tick label space).
    fig.subplots_adjust(right=1 - 3.2 / fig.get_figwidth(), bottom=TICK_LABEL_HEIGHT_IN / fig.get_figheight())
    fig.savefig(group["fig_dir"] / "evenness_summary.png", dpi=150)  # write the PNG at 150 dpi
    plt.close(fig)  # free the figure's memory (several groups are drawn in one run)
    # Progress message on stderr so it doesn't mix with any stdout output.
    print(f"Wrote {group['fig_dir'] / 'evenness_summary.png'}", file=sys.stderr)


# main: load inputs once, then draw the figure for every sample group.
#   In:  optional CLI args [samples_tsv] [samples_meta_csv] (read inside load_samples()),
#        plus <run_dir>/reference/combined_metadata.csv and <run_dir>/results/count_matrix.csv.
#   Out: evenness_summary.png in results/figures/ and in each results/figures/<ip>/.
def main():
    # Merged, sorted sample sheet (tc_id, library, ip, treatment, replicate) + treatment colours.
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)  # dtype=str: keep every column as text (IDs stay exact)
    # Series indexed by ref_id whose values are the library name -- a ref_id -> library lookup.
    lib_of_ref = meta.set_index("ref_id")["library"]
    # Raw (unfiltered) peptide x sample count matrix; index_col=0 makes ref_id the row index.
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)

    # figure_groups() returns the "all samples" group plus one group per IP (and creates
    # their output folders); draw one figure per group.
    for group in figure_groups(samples, FIG_DIR):
        plot_evenness(mat, lib_of_ref, group, colors)


# Only run main() when executed as a script, not when imported as a module.
if __name__ == "__main__":
    main()
