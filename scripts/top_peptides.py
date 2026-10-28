"""Peptide-level views of the count matrix:
  1. top_variable_heatmap.png -- the N most variable peptides (by variance
     of log10(CPM + 1) across the group's samples), rows clustered,
     columns in treatment order with a treatment colour strip. Drawn once
     for all samples (results/figures/) and once per IP condition
     (results/figures/<ip>/, see sample_groups.py), each group picking its
     own top N.
  2. results/top_abundant_peptides.csv -- the N most abundant peptides in
     each sample (long format: one row per sample x rank), with read-pair
     count, CPM, fraction of the sample's reads and peptide annotation;
     results/top_abundant_peptides_wide.csv -- the same as a rank x sample
     grid of short peptide labels, for scanning by eye.

The heatmap applies the matrix filter within each group, on that group's
samples only (count_filters.passes_filter); the tables use the
all-sample filtered matrix (count_matrix_filtered.csv). Counts are CPM-normalised first -- per sample, against
its full own-library read total from the raw matrix, so CPM means the same
as before filtering -- since sequencing depth differs several-fold between
samples.

Usage: top_peptides.py [samples_tsv] [samples_meta_csv]
  (default to reference/samples.tsv and reference/samples_meta.csv)
"""
# --- Imports -----------------------------------------------------------------
import json  # to read the reference/*.json config files
import sys  # sys.exit for fatal errors, sys.stderr for progress messages
from pathlib import Path  # object-oriented file paths ("/" joins path parts)

import numpy as np  # vectorised maths (log10, arrays for the colour strip)
import pandas as pd  # DataFrames: count matrices, metadata, output tables
import matplotlib

# Select the non-interactive "Agg" backend BEFORE importing pyplot so figures
# render straight to PNG without a display (e.g. on a cluster node).
matplotlib.use("Agg")
import matplotlib.colors as mcolors  # colour-name/tuple -> RGBA conversion
import matplotlib.pyplot as plt  # the plotting interface (figures, axes)
# make_axes_locatable lets us attach extra thin axes (colour strip, colorbar)
# to the sides of the heatmap axes, sized in inches and aligned to it.
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Hierarchical clustering: linkage() builds the cluster tree; leaves_list()
# returns the row order of the tree's leaves (similar rows end up adjacent).
from scipy.cluster.hierarchy import leaves_list, linkage

# passes_filter(counts) -> boolean per peptide: >= matrix_filter_min_count reads
# in >= matrix_filter_min_samples of the given samples (filter_params.json).
# GROUP_FILTER_DESCRIPTION is the human-readable form used in the title.
from count_filters import GROUP_FILTER_DESCRIPTION, passes_filter
# Shared sample-sheet loading, figure grouping (all samples / per IP) and
# label/legend helpers, so every figure script looks the same.
from sample_groups import figure_groups, load_samples, sample_label, treatment_legend_handles

# Repo root: this file is scripts/top_peptides.py, so go up two levels.
PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())  # dict of paths
RUN_DIR = Path(_paths["run_dir"])  # the data directory (may live anywhere)
RESULTS_DIR = RUN_DIR / "results"  # count matrices are read, tables written, here
FIG_DIR = RESULTS_DIR / "figures"  # all-samples heatmap here; per-IP ones in FIG_DIR/<ip>/
# Per-peptide annotation (ref_id -> library, organism, protein, ...) from build_reference.py.
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

# Figure settings shared with the other figure scripts (reference/figure_params.json).
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
N_TOP_VARIABLE = _fig_params["heatmap_n_top_variable_peptides"]  # heatmap rows per group
N_TOP_ABUNDANT = _fig_params["top_abundant_n_peptides"]  # table rows per sample
HEATMAP_TICK_FONTSIZE = _fig_params["heatmap_tick_fontsize"]  # base size; x labels use 75% of it

# Heatmap geometry: figure size is derived from the number of rows/columns so
# that each cell gets a fixed physical size and labels stay readable.
ROW_LABEL_FONTSIZE = 7  # font size of the peptide labels on the y axis
ROW_HEIGHT_IN = 0.14  # per peptide row: keeps ROW_LABEL_FONTSIZE labels from overlapping
COL_WIDTH_IN = 0.35  # inches per sample column


# ---------------------------------------------------------------------------
# short_name(text, max_len, strip_parens=True) -> str
#   In:  an annotation string (organism/protein name; may be NaN), max length.
#   Out: a tidied, length-capped label ("unknown" if missing).
#   Used by peptide_label() here and imported by differential_enrichment.py
#   (volcano annotations, organism bar-chart labels).
# ---------------------------------------------------------------------------
def short_name(text, max_len, strip_parens=True):
    """Annotation trimmed for an axis label: underscores -> spaces, then
    truncate. strip_parens drops parenthesised synonyms / EC numbers --
    turn it off where they are the only thing distinguishing two entries
    (e.g. virus strains sharing a species name)."""
    # Missing pandas values arrive as float NaN (not str), and "" is falsy.
    if not isinstance(text, str) or not text:
        return "unknown"
    text = text.replace("_", " ")  # underscores -> spaces for readability
    if strip_parens:
        # Keep only the part before the first " (" e.g. "Spike (S protein)" -> "Spike".
        text = text.split(" (")[0]
    text = text.strip()  # trim leading/trailing whitespace
    # If too long, cut to max_len - 1 characters and add an ellipsis (total max_len).
    return text if len(text) <= max_len else text[:max_len - 1] + "…"


# ---------------------------------------------------------------------------
# peptide_label(ann) -> str like "source_id | organism | protein start-end"
#   In:  one peptide's annotation row (a Series from combined_metadata.csv).
#   Used for heatmap row labels and table columns here, and imported by
#   differential_enrichment.py for its per-contrast tables.
# ---------------------------------------------------------------------------
def peptide_label(ann):
    """Short one-line peptide label, shared by the heatmap rows and the tables."""
    # Organism capped at 32 chars, protein name at 24; start/end are the
    # peptide's position within the protein.
    return (f"{ann['source_id']} | {short_name(ann['organism'], 32)} | "
            f"{short_name(ann['protein_name'], 24)} {ann['start']}-{ann['end']}")


# peptide_labels(meta, ref_ids) -> list of labels, one per ref_id, in order.
# meta must be indexed by ref_id. Called by draw_peptide_heatmap().
def peptide_labels(meta, ref_ids):
    return [peptide_label(meta.loc[ref_id]) for ref_id in ref_ids]  # meta.loc[id] = one row


# ---------------------------------------------------------------------------
# cpm_matrix(mat_filtered, mat_raw, lib_of_ref, samples) -> DataFrame
#   In:  mat_filtered = the counts to normalise (raw or filtered matrix);
#        mat_raw = raw matrix used for each sample's library size (denominator);
#        lib_of_ref = ref_id -> library Series; samples = sample-sheet rows.
#   Out: own-library peptides x samples matrix of counts per million (CPM).
#   Called by main() twice (raw matrix for the heatmap, filtered for the tables).
#
# CPM = count / (sample's total reads) * 1,000,000 -- puts samples sequenced
# to different depths on a common scale. The denominator always comes from the
# RAW matrix so a peptide's CPM is the same whether or not filtering ran.
# ---------------------------------------------------------------------------
def cpm_matrix(mat_filtered, mat_raw, lib_of_ref, samples):
    """Counts per million for the filtered matrix's own-library peptides,
    scaled by each sample's raw own-library read total; all samples here
    share one library, so this is one peptides x samples frame
    (other-library peptides dropped)."""
    libraries = samples["library"].unique()  # distinct libraries in the sample sheet
    if len(libraries) != 1:
        # One shared denominator set only makes sense for a single library.
        sys.exit(f"ERROR: top_peptides.py expects samples from one library, got {', '.join(libraries)}")
    own = lib_of_ref[lib_of_ref == libraries[0]].index  # ref_ids of that library (boolean mask)
    # Per-sample library size: sum down the rows (axis=0 -> one total per column)
    # of the raw counts over own-library peptides present in the raw matrix.
    lib_sizes = mat_raw.loc[own.intersection(mat_raw.index), samples["tc_id"]].sum(axis=0)
    # DataFrame / Series divides each column by the Series value with the
    # matching label (pandas aligns the Series index to the columns), then x 1e6.
    return mat_filtered.loc[own.intersection(mat_filtered.index), samples["tc_id"]] / lib_sizes * 1e6


# ---------------------------------------------------------------------------
# draw_peptide_heatmap(values, meta, group, colors, title, out_path,
#                      cbar_label="log10(CPM + 1)") -> None (writes a PNG)
#   In:  values = peptides x samples DataFrame of values to colour (normally
#        log10(CPM + 1)); meta = annotation indexed by ref_id; group = a
#        figure-group dict (needs "samples" and "name"); colors = treatment ->
#        colour; title; out_path = PNG to write.
#   Out: the heatmap PNG -- rows clustered, columns in sample order, a
#        treatment colour strip on top, a colorbar on the left, peptide
#        labels on the right.
#   Called by plot_top_variable_heatmap() here and by
#   differential_enrichment.plot_enriched_heatmap().
# ---------------------------------------------------------------------------
def draw_peptide_heatmap(values, meta, group, colors, title, out_path, cbar_label="log10(CPM + 1)"):
    """Peptides x samples heatmap (rows clustered, columns in the group's
    sample order under a treatment colour strip, rows labelled with
    peptide_label). Shared with differential_enrichment.py."""
    samples = group["samples"]  # sample-sheet rows for this group, in plotting order
    values = values[samples["tc_id"]]  # select + reorder columns to match that order
    # Cluster rows only when there are enough to make ordering meaningful.
    # linkage(): agglomerative clustering on Euclidean distance between row
    # profiles, merging clusters by "average" linkage (UPGMA: mean pairwise
    # distance); optimal_ordering=True flips branches so neighbouring leaves
    # are as similar as possible. leaves_list() gives the resulting row order,
    # and .iloc reorders the rows by position.
    if len(values) > 2:
        values = values.iloc[leaves_list(linkage(values.to_numpy(), method="average", metric="euclidean",
                                                 optimal_ordering=True))]
    n_rows, n_cols = values.shape  # peptides, samples

    # Figure size scales with the grid; the +7 / +3 inches leave room for labels,
    # colorbar and legend.
    fig, ax = plt.subplots(figsize=(COL_WIDTH_IN * n_cols + 7, ROW_HEIGHT_IN * n_rows + 3))
    # imshow draws the 2-D array as coloured cells (row 0 at the top).
    # aspect="auto" lets cells be non-square to fill the axes; "nearest" keeps
    # sharp cell edges (no smoothing). The returned image feeds the colorbar.
    im = ax.imshow(values.to_numpy(), aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_xticks(range(n_cols))  # one tick per sample column
    # Sample labels, vertical, at 75% of the configured heatmap tick font size.
    ax.set_xticklabels([sample_label(srow, group) for _, srow in samples.iterrows()],
                       rotation=90, fontsize=HEATMAP_TICK_FONTSIZE * 0.75)
    ax.set_yticks(range(n_rows))  # one tick per peptide row
    # Peptide labels in the (clustered) row order.
    ax.set_yticklabels(peptide_labels(meta, values.index), fontsize=ROW_LABEL_FONTSIZE)
    ax.yaxis.tick_right()  # put the (long) peptide labels on the right-hand side

    # Divider carves extra axes off the heatmap's edges, aligned with it.
    divider = make_axes_locatable(ax)
    # Thin (0.25 in) strip above the heatmap; sharex keeps its columns aligned.
    ax_top = divider.append_axes("top", size=0.25, pad=0.05, sharex=ax)
    # One RGBA colour per sample (its treatment colour) -> array of shape (n_cols, 4).
    strip = np.array([mcolors.to_rgba(colors[t]) for t in samples["treatment"]])
    # np.newaxis adds a row dimension -> (1, n_cols, 4): a 1-pixel-high RGBA
    # image that imshow shows directly (RGBA arrays need no colormap).
    ax_top.imshow(strip[np.newaxis, :, :], aspect="auto")
    ax_top.set_yticks([])  # no y ticks on the strip
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)  # nor x ticks/labels
    ax_top.set_title(title, fontsize=11)  # the title sits above the strip

    # Narrow axes on the left for the colorbar (0.6 in gap, clear of the heatmap).
    cax = divider.append_axes("left", size=0.2, pad=0.6)
    fig.colorbar(im, cax=cax, label=cbar_label)  # colour scale for the heatmap values
    cax.yaxis.set_ticks_position("left")  # colorbar ticks on its outer (left) side
    cax.yaxis.set_label_position("left")  # and its label too
    # below the sample tick labels, clear of the long peptide labels on the right
    # [1:] drops the "— treatment —" header entry (legend title used instead).
    # bbox_to_anchor y is in axes fractions: -0.02 minus roughly 2.2 inches
    # expressed as a fraction of figure height (approximates the label space).
    ax.legend(handles=treatment_legend_handles(samples, colors)[1:], title="treatment", loc="upper center",
              bbox_to_anchor=(0.5, -0.02 - 2.2 / fig.get_figheight()), ncol=3, fontsize=9)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")  # "tight" includes outside legend/labels
    plt.close(fig)  # release memory
    print(f"Wrote {out_path}", file=sys.stderr)  # progress message


# ---------------------------------------------------------------------------
# plot_top_variable_heatmap(log_cpm, mat_raw, meta, group, colors) -> None
#   In:  log_cpm = log10(CPM + 1) for all own-library peptides x all samples;
#        mat_raw = raw counts (for the within-group filter); meta, group, colors.
#   Out: <group fig_dir>/top_variable_heatmap.png via draw_peptide_heatmap().
#   Called by main() once per figure group.
#
# Picks the N peptides whose log10(CPM + 1) varies most across THIS group's
# samples -- the ones that best distinguish samples. log10(x + 1) compresses
# the huge dynamic range (so a few very abundant peptides don't dominate the
# variance) and the +1 avoids log10(0).
# ---------------------------------------------------------------------------
def plot_top_variable_heatmap(log_cpm, mat_raw, meta, group, colors):
    samples = group["samples"]  # this group's samples
    sub = log_cpm[samples["tc_id"]]  # keep only this group's columns
    # Within-group low-count filter: passes_filter returns a boolean Series
    # (one per peptide) computed on raw counts of just these samples; .loc
    # with that mask keeps the passing rows.
    sub = sub.loc[passes_filter(mat_raw.loc[sub.index, samples["tc_id"]])]
    n_eligible = len(sub)  # number of peptides passing (reported in the title)
    # Variance across samples for each row (axis=1 = along columns), largest
    # first; take the top N row labels and select those rows.
    top = sub.loc[sub.var(axis=1).sort_values(ascending=False).head(N_TOP_VARIABLE).index]
    # Two-line title: what is shown, then how the peptides were chosen.
    title = (f"Top {len(top)} most variable peptides{group['title_suffix']}\n"
             f"(variance of log10(CPM + 1) among {n_eligible} peptides with {GROUP_FILTER_DESCRIPTION}; "
             "rows clustered)")
    draw_peptide_heatmap(top, meta, group, colors, title, group["fig_dir"] / "top_variable_heatmap.png")


# ---------------------------------------------------------------------------
# top_abundant_table(mat, cpm, meta, samples) -> long-format DataFrame
#   In:  mat = filtered count matrix (for raw read-pair counts); cpm = CPM of
#        that matrix (from cpm_matrix); meta = annotation indexed by ref_id;
#        samples = sample-sheet rows.
#   Out: one row per (sample, rank) for the N_TOP_ABUNDANT highest-CPM
#        peptides of each sample, with counts and annotation.
#   Called by main().
# ---------------------------------------------------------------------------
def top_abundant_table(mat, cpm, meta, samples):
    rows = []  # list of dicts -> DataFrame at the end (faster than appending to a frame)
    for _, srow in samples.iterrows():
        tc_id = srow["tc_id"]  # this sample's column name
        # This sample's CPM column, sorted high to low, first N entries.
        top = cpm[tc_id].sort_values(ascending=False).head(N_TOP_ABUNDANT)
        # .items() yields (ref_id, cpm) pairs; enumerate(start=1) numbers them 1..N.
        for rank, (ref_id, sample_cpm) in enumerate(top.items(), start=1):
            ann = meta.loc[ref_id]  # this peptide's annotation row
            rows.append({
                "tc_id": tc_id,
                "ip": srow["ip"],
                "treatment": srow["treatment"],
                "rank": rank,
                "ref_id": ref_id,
                "peptide_label": peptide_label(ann),  # short human-readable label
                "read_pairs": int(mat.at[ref_id, tc_id]),  # .at = fast single-cell lookup
                "cpm": sample_cpm,
                # CPM / 1e6 = fraction of the sample's (raw own-library) reads.
                "fraction_of_sample_reads": sample_cpm / 1e6,
                # Full (untrimmed) annotation for reference.
                "organism": ann["organism"],
                "protein_name": ann["protein_name"],
                "start": ann["start"],
                "end": ann["end"],
                "peptide_aa": ann["peptide_aa"],  # amino-acid sequence
            })
    return pd.DataFrame(rows)  # each dict becomes one row; keys become columns


# ---------------------------------------------------------------------------
# main() -- entry point.
#   In:  optional CLI args (sample sheet, sample metadata; read by
#        sample_groups.load_samples), combined_metadata.csv,
#        results/count_matrix.csv and results/count_matrix_filtered.csv
#        (both written by summarize.py).
#   Out: top_variable_heatmap.png per figure group, plus
#        results/top_abundant_peptides.csv and top_abundant_peptides_wide.csv.
# ---------------------------------------------------------------------------
def main():
    # samples: merged, plot-ordered sample sheet; colors: treatment -> colour.
    samples, colors = load_samples()
    # Annotation table, all columns as text, indexed by ref_id for fast lookups.
    meta = pd.read_csv(META_CSV, dtype=str).set_index("ref_id")
    lib_of_ref = meta["library"]  # Series: ref_id -> library name
    # Raw counts (every peptide) and the all-sample filtered counts; the first
    # column (ref_id) becomes the row index in both.
    mat_raw = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)
    mat = pd.read_csv(RESULTS_DIR / "count_matrix_filtered.csv", index_col=0)

    # heatmap: raw matrix, filtered within each group by plot_top_variable_heatmap
    # log10(CPM + 1), element-wise over the whole DataFrame.
    log_cpm_raw = np.log10(cpm_matrix(mat_raw, mat_raw, lib_of_ref, samples) + 1)
    # One heatmap for all samples, then one per IP (each in its own folder).
    for group in figure_groups(samples, FIG_DIR):
        plot_top_variable_heatmap(log_cpm_raw, mat_raw, meta, group, colors)

    # CPM of the filtered matrix (denominators still from the raw matrix).
    cpm = cpm_matrix(mat, mat_raw, lib_of_ref, samples)

    table = top_abundant_table(mat, cpm, meta, samples)  # long format
    out = RESULTS_DIR / "top_abundant_peptides.csv"
    table.to_csv(out, index=False)  # index=False: don't write the 0..n row numbers
    # Reshape long -> wide: one row per rank, one column per sample, cells =
    # peptide labels; [samples["tc_id"]] puts columns in plotting order.
    wide = table.pivot(index="rank", columns="tc_id", values="peptide_label")[samples["tc_id"]]
    # Friendlier column headers: "<ip> <treatment> <tc_id>" (same order as above).
    wide.columns = [f"{srow['ip']} {srow['treatment']} {srow['tc_id']}" for _, srow in samples.iterrows()]
    wide.to_csv(RESULTS_DIR / "top_abundant_peptides_wide.csv")  # keeps "rank" as the first column
    print(f"Wrote {out} and {RESULTS_DIR / 'top_abundant_peptides_wide.csv'}", file=sys.stderr)


# Run main() only when executed as a script, not when imported (e.g. by
# differential_enrichment.py, which imports the heatmap/label helpers).
if __name__ == "__main__":
    main()
