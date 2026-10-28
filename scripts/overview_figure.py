"""One combined figure comparing samples together: coverage (%% of
designed peptides detected at >=1 read and at >= coverage_min_count reads,
see count_filters.py) and evenness (typical range + near-dropout
outliers) side by side, coloured by treatment. Drawn once for all samples
(results/figures/) and once per IP condition (results/figures/<ip>/, see
sample_groups.py).

Usage: overview_figure.py [samples_tsv] [samples_meta_csv]
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
# Line2D / Patch let us build "fake" artists purely to use as legend entries
# (Line2D for lines/markers, Patch for filled rectangles matching the bars).
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Coverage threshold from reference/filter_params.json: a peptide counts as
# "detected" for the stricter coverage metric only with >= this many reads.
from count_filters import COVERAGE_MIN_COUNT
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

# Figure thresholds/styling live in reference/figure_params.json (shared
# with evenness_summary_figure.py) rather than being duplicated as literals
# here.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
OUTLIER_THRESHOLD = _fig_params["outlier_threshold"]  # reads
XTICK_ROTATION = _fig_params["xtick_rotation"]  # degrees to rotate x tick labels (90 = vertical)

# Near-dropout markers are black so they stay distinct from the treatment palette.
OUTLIER_COLOR = "black"

# Room below the evenness panel for the rotated "<ip> <tc_id> (N ND, N miss)" tick labels.
TICK_LABEL_HEIGHT_IN = 3.2


# plot_overview: draw and save one two-panel overview figure for one sample group.
#   In:  mat        -- raw count matrix (rows = peptide ref_ids, columns = sample tc_ids)
#        lib_of_ref -- Series mapping ref_id -> library name ("CoV"/"Vir3")
#        group      -- dict from figure_groups(): samples (DataFrame), fig_dir, title_suffix, name
#        colors     -- dict treatment -> matplotlib colour
#   Out: nothing returned; writes <group fig_dir>/overview_all_samples.png.
#   Top panel: coverage bars (% of own-library peptides with >=1 read, light, and with
#   >= COVERAGE_MIN_COUNT reads, solid). Bottom panel: the same evenness plot as
#   evenness_summary_figure.py (P10-P90 bar, median dot, near-dropout markers).
def plot_overview(mat, lib_of_ref, group, colors):
    samples = group["samples"]  # the sample-sheet rows in this group, already in plotting order
    # Two axes stacked vertically (2 rows, 1 column): unpacked into ax_cov (top) and ax_even
    # (bottom). Width grows 0.7 in per sample; height = 8.5 in of panels + tick label space.
    # Note the x axes are not shared (no sharex) -- each panel sets its own ticks below.
    fig, (ax_cov, ax_even) = plt.subplots(2, 1, figsize=(0.7 * len(samples) + 5, 8.5 + TICK_LABEL_HEIGHT_IN))

    cov_labels, even_labels = [], []  # tick labels for each panel, filled in the loop
    # iterrows() yields (index, row) pairs; enumerate adds i = 0, 1, 2... used as the x position.
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]  # sample ID and the library it was grown from
        color = colors[srow["treatment"]]  # this sample's treatment colour (same in both panels)
        # "Own-library" peptides: boolean mask picks ref_ids whose library matches this
        # sample's, .index gives those ref_ids, and intersection() keeps only ones that are
        # actually rows of the count matrix (so the row lookup in .loc below cannot raise a KeyError).
        # Peptides from the other library are ignored -- they are not expected in this sample.
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]  # Series: read count of each own-library peptide
        total = len(own_refs)  # denominator: number of designed peptides in this library
        detected = counts[counts > 0]  # boolean mask: peptides seen at least once
        # Peptides passing the stricter coverage threshold: True/False per peptide, summed.
        n_min = int((counts >= COVERAGE_MIN_COUNT).sum())
        # Coverage percentages; the "if total else 0" guards against dividing by zero
        # when no own-library peptides are found in the matrix.
        pct = 100 * len(detected) / total if total else 0
        pct_min = 100 * n_min / total if total else 0

        # --- coverage panel: light bar >=1 read, solid bar >= COVERAGE_MIN_COUNT reads ---
        # Both bars sit at the same x, so the shorter solid bar overlays the taller light one
        # (alpha=0.35 makes the >=1-read bar pale); thin black edges outline each bar.
        ax_cov.bar(i, pct, color=color, alpha=0.35, edgecolor="black", linewidth=0.5)
        ax_cov.bar(i, pct_min, color=color, edgecolor="black", linewidth=0.5)
        # Write the >=COVERAGE_MIN_COUNT peptide count just above the top of the (taller) light
        # bar, rotated vertical; va="bottom" anchors the text's bottom edge at that y.
        ax_cov.text(i, pct + 1.5, f"{n_min}", ha="center", va="bottom", rotation=90, fontsize=7)

        # --- evenness panel ---
        # Near-dropouts: detected (>0) but only at a very low count (<= threshold reads).
        outliers = detected[detected <= OUTLIER_THRESHOLD]
        # Only draw if something was detected (percentiles of an empty array would fail).
        if len(detected):
            # 10th, 50th (median) and 90th percentiles of the detected counts in one call;
            # P10-P90 is the "typical range" covering the middle 80% of peptides.
            p10, p50, p90 = np.percentile(detected, [10, 50, 90])
            # Thick vertical bar from P10 to P90; "butt" caps end it exactly at those values.
            ax_even.plot([i, i], [p10, p90], color=color, linewidth=6, solid_capstyle="butt", alpha=0.85)
            # Median as a white circle with a black edge, drawn on top of the bar (zorder=5).
            ax_even.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=7, zorder=5)
            # Each near-dropout peptide as a black down-triangle at this sample's x
            # (np.full repeats i once per outlier); linestyle="" = markers only.
            ax_even.plot(np.full(len(outliers), i), outliers.values, "v", color=OUTLIER_COLOR, markersize=8,
                         linestyle="", zorder=6)

        # Dropouts: designed peptides with zero reads in this sample.
        n_missing = int((counts == 0).sum())
        # Kept as a single line: rotated multi-line tick labels lay out very
        # differently from rotated single-line ones (each line becomes its
        # own column) and need much more room than is practical here.
        # List comprehension builds "N ND" / "N miss" strings (or "" when zero) and the
        # trailing "if b" keeps only the non-empty ones.
        extra_bits = [b for b in [f"{len(outliers)} ND" if len(outliers) else "",
                                  f"{n_missing} miss" if n_missing else ""] if b]
        label = sample_label(srow, group)  # "<ip> <tc_id>" in the all-samples figure, else tc_id
        cov_labels.append(label)  # coverage panel: plain label
        # Evenness panel: label plus " (N ND, N miss)" when there is anything to report
        # (unlike evenness_summary_figure.py, no " (clean)" suffix otherwise).
        even_labels.append(label + (" (" + ", ".join(extra_bits) + ")" if extra_bits else ""))

    # --- coverage panel axes ---
    ax_cov.set_ylim(0, 125)  # y up to 125% (not 100) leaves room for the count labels above bars
    ax_cov.set_xticks(range(len(samples)))  # one tick per sample at x = 0, 1, 2, ...
    ax_cov.set_xticklabels(cov_labels, rotation=XTICK_ROTATION, ha="center")  # sample names, rotated
    ax_cov.set_ylabel("% of designed peptides detected")
    # NOTE: `total` here is whatever the loop left behind, i.e. the LAST sample's library size.
    ax_cov.set_title(f"Coverage (bar labels: peptides with >= {COVERAGE_MIN_COUNT} reads, of {total} designed)")
    # Legend right of the coverage panel: treatment header + colour squares, then a
    # "detected at" header (an invisible Line2D whose only purpose is its label) and two
    # gray Patch swatches mimicking the light (>=1 read) and solid (>=min count) bars.
    legend_outside(ax_cov, treatment_legend_handles(samples, colors) + [
        Line2D([], [], linestyle="", label="— detected at —"),
        Patch(facecolor="gray", alpha=0.35, edgecolor="black", label=">= 1 read"),
        Patch(facecolor="gray", edgecolor="black", label=f">= {COVERAGE_MIN_COUNT} reads"),
    ], fontsize=8)

    # --- evenness panel axes ---
    # Log y-axis: read counts span several orders of magnitude (a few reads to thousands).
    ax_even.set_yscale("log")
    # Half a slot of padding either side of the first (x=0) and last (x=n-1) sample.
    ax_even.set_xlim(-0.6, len(samples) - 0.4)
    ax_even.set_xticks(range(len(samples)))  # one tick per sample
    # Sample labels with ND/miss counts, rotated so long labels don't overlap.
    ax_even.set_xticklabels(even_labels, fontsize=8, rotation=XTICK_ROTATION, ha="center")
    ax_even.set_ylabel("Reads per peptide (log scale)")
    ax_even.set_title(f"Evenness (10th-90th percentile + outliers)\n(ND = near-dropout, <={OUTLIER_THRESHOLD} reads)")
    # Legend right of the evenness panel, built from empty Line2D "proxy" artists that only
    # carry a style + label (treatment colours are already explained in the top legend).
    legend_outside(ax_even, [
        Line2D([], [], color="gray", linewidth=6, label="typical range (10th-90th percentile)"),
        Line2D([], [], marker="o", linestyle="", color="white", markeredgecolor="black", markersize=7,
               label="median"),
        Line2D([], [], marker="v", linestyle="", color=OUTLIER_COLOR, markersize=8,
               label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)"),
    ], fontsize=8)

    # Overall figure title above both panels, with the group suffix (e.g. " (IP4)").
    fig.suptitle("Library completeness and evenness overview" + group["title_suffix"], fontsize=13)
    # Fixed widths/heights (in inches, converted to figure fractions) are
    # reserved for the outside legends, the rotated tick labels, and a gap
    # between the two stacked panels (their rotated labels/title would
    # otherwise collide). Not using
    # bbox_inches="tight" here: with this matplotlib version it can blow up
    # the canvas to a pathological size when combined with a log-scaled axis
    # and many plotted markers.
    # right: leave 3.2 in for legends; bottom: tick label space; top: 0.9 in for the suptitle;
    # hspace=0.9: vertical gap between panels as a fraction of the average panel height.
    fig.subplots_adjust(right=1 - 3.2 / fig.get_figwidth(), bottom=TICK_LABEL_HEIGHT_IN / fig.get_figheight(),
                        top=1 - 0.9 / fig.get_figheight(), hspace=0.9)
    # Same file name in every group folder (the per-IP copy covers only that IP's samples).
    fig.savefig(group["fig_dir"] / "overview_all_samples.png", dpi=150)
    plt.close(fig)  # free the figure's memory (several groups are drawn in one run)
    # Progress message on stderr so it doesn't mix with any stdout output.
    print(f"Wrote {group['fig_dir'] / 'overview_all_samples.png'}", file=sys.stderr)


# main: load inputs once, then draw the overview figure for every sample group.
#   In:  optional CLI args [samples_tsv] [samples_meta_csv] (read inside load_samples()),
#        plus <run_dir>/reference/combined_metadata.csv and <run_dir>/results/count_matrix.csv.
#   Out: overview_all_samples.png in results/figures/ and in each results/figures/<ip>/.
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
        plot_overview(mat, lib_of_ref, group, colors)


# Only run main() when executed as a script, not when imported as a module.
if __name__ == "__main__":
    main()
