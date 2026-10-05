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
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import leaves_list, linkage

from count_filters import GROUP_FILTER_DESCRIPTION, passes_filter
from sample_groups import figure_groups, load_samples, sample_label, treatment_legend_handles

PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())
RUN_DIR = Path(_paths["run_dir"])
RESULTS_DIR = RUN_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
N_TOP_VARIABLE = _fig_params["heatmap_n_top_variable_peptides"]
N_TOP_ABUNDANT = _fig_params["top_abundant_n_peptides"]
HEATMAP_TICK_FONTSIZE = _fig_params["heatmap_tick_fontsize"]

ROW_LABEL_FONTSIZE = 7
ROW_HEIGHT_IN = 0.14  # per peptide row: keeps ROW_LABEL_FONTSIZE labels from overlapping
COL_WIDTH_IN = 0.35


def short_name(text, max_len, strip_parens=True):
    """Annotation trimmed for an axis label: underscores -> spaces, then
    truncate. strip_parens drops parenthesised synonyms / EC numbers --
    turn it off where they are the only thing distinguishing two entries
    (e.g. virus strains sharing a species name)."""
    if not isinstance(text, str) or not text:
        return "unknown"
    text = text.replace("_", " ")
    if strip_parens:
        text = text.split(" (")[0]
    text = text.strip()
    return text if len(text) <= max_len else text[:max_len - 1] + "…"


def peptide_label(ann):
    """Short one-line peptide label, shared by the heatmap rows and the tables."""
    return (f"{ann['source_id']} | {short_name(ann['organism'], 32)} | "
            f"{short_name(ann['protein_name'], 24)} {ann['start']}-{ann['end']}")


def peptide_labels(meta, ref_ids):
    return [peptide_label(meta.loc[ref_id]) for ref_id in ref_ids]


def cpm_matrix(mat_filtered, mat_raw, lib_of_ref, samples):
    """Counts per million for the filtered matrix's own-library peptides,
    scaled by each sample's raw own-library read total; all samples here
    share one library, so this is one peptides x samples frame
    (other-library peptides dropped)."""
    libraries = samples["library"].unique()
    if len(libraries) != 1:
        sys.exit(f"ERROR: top_peptides.py expects samples from one library, got {', '.join(libraries)}")
    own = lib_of_ref[lib_of_ref == libraries[0]].index
    lib_sizes = mat_raw.loc[own.intersection(mat_raw.index), samples["tc_id"]].sum(axis=0)
    return mat_filtered.loc[own.intersection(mat_filtered.index), samples["tc_id"]] / lib_sizes * 1e6


def draw_peptide_heatmap(values, meta, group, colors, title, out_path, cbar_label="log10(CPM + 1)"):
    """Peptides x samples heatmap (rows clustered, columns in the group's
    sample order under a treatment colour strip, rows labelled with
    peptide_label). Shared with differential_enrichment.py."""
    samples = group["samples"]
    values = values[samples["tc_id"]]
    if len(values) > 2:
        values = values.iloc[leaves_list(linkage(values.to_numpy(), method="average", metric="euclidean",
                                                 optimal_ordering=True))]
    n_rows, n_cols = values.shape

    fig, ax = plt.subplots(figsize=(COL_WIDTH_IN * n_cols + 7, ROW_HEIGHT_IN * n_rows + 3))
    im = ax.imshow(values.to_numpy(), aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([sample_label(srow, group) for _, srow in samples.iterrows()],
                       rotation=90, fontsize=HEATMAP_TICK_FONTSIZE * 0.75)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(peptide_labels(meta, values.index), fontsize=ROW_LABEL_FONTSIZE)
    ax.yaxis.tick_right()

    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size=0.25, pad=0.05, sharex=ax)
    strip = np.array([mcolors.to_rgba(colors[t]) for t in samples["treatment"]])
    ax_top.imshow(strip[np.newaxis, :, :], aspect="auto")
    ax_top.set_yticks([])
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_top.set_title(title, fontsize=11)

    cax = divider.append_axes("left", size=0.2, pad=0.6)
    fig.colorbar(im, cax=cax, label=cbar_label)
    cax.yaxis.set_ticks_position("left")
    cax.yaxis.set_label_position("left")
    # below the sample tick labels, clear of the long peptide labels on the right
    ax.legend(handles=treatment_legend_handles(samples, colors)[1:], title="treatment", loc="upper center",
              bbox_to_anchor=(0.5, -0.02 - 2.2 / fig.get_figheight()), ncol=3, fontsize=9)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def plot_top_variable_heatmap(log_cpm, mat_raw, meta, group, colors):
    samples = group["samples"]
    sub = log_cpm[samples["tc_id"]]
    sub = sub.loc[passes_filter(mat_raw.loc[sub.index, samples["tc_id"]])]
    n_eligible = len(sub)
    top = sub.loc[sub.var(axis=1).sort_values(ascending=False).head(N_TOP_VARIABLE).index]
    title = (f"Top {len(top)} most variable peptides{group['title_suffix']}\n"
             f"(variance of log10(CPM + 1) among {n_eligible} peptides with {GROUP_FILTER_DESCRIPTION}; "
             "rows clustered)")
    draw_peptide_heatmap(top, meta, group, colors, title, group["fig_dir"] / "top_variable_heatmap.png")


def top_abundant_table(mat, cpm, meta, samples):
    rows = []
    for _, srow in samples.iterrows():
        tc_id = srow["tc_id"]
        top = cpm[tc_id].sort_values(ascending=False).head(N_TOP_ABUNDANT)
        for rank, (ref_id, sample_cpm) in enumerate(top.items(), start=1):
            ann = meta.loc[ref_id]
            rows.append({
                "tc_id": tc_id,
                "ip": srow["ip"],
                "treatment": srow["treatment"],
                "rank": rank,
                "ref_id": ref_id,
                "peptide_label": peptide_label(ann),
                "read_pairs": int(mat.at[ref_id, tc_id]),
                "cpm": sample_cpm,
                "fraction_of_sample_reads": sample_cpm / 1e6,
                "organism": ann["organism"],
                "protein_name": ann["protein_name"],
                "start": ann["start"],
                "end": ann["end"],
                "peptide_aa": ann["peptide_aa"],
            })
    return pd.DataFrame(rows)


def main():
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str).set_index("ref_id")
    lib_of_ref = meta["library"]
    mat_raw = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)
    mat = pd.read_csv(RESULTS_DIR / "count_matrix_filtered.csv", index_col=0)

    # heatmap: raw matrix, filtered within each group by plot_top_variable_heatmap
    log_cpm_raw = np.log10(cpm_matrix(mat_raw, mat_raw, lib_of_ref, samples) + 1)
    for group in figure_groups(samples, FIG_DIR):
        plot_top_variable_heatmap(log_cpm_raw, mat_raw, meta, group, colors)

    cpm = cpm_matrix(mat, mat_raw, lib_of_ref, samples)

    table = top_abundant_table(mat, cpm, meta, samples)
    out = RESULTS_DIR / "top_abundant_peptides.csv"
    table.to_csv(out, index=False)
    wide = table.pivot(index="rank", columns="tc_id", values="peptide_label")[samples["tc_id"]]
    wide.columns = [f"{srow['ip']} {srow['treatment']} {srow['tc_id']}" for _, srow in samples.iterrows()]
    wide.to_csv(RESULTS_DIR / "top_abundant_peptides_wide.csv")
    print(f"Wrote {out} and {RESULTS_DIR / 'top_abundant_peptides_wide.csv'}", file=sys.stderr)


if __name__ == "__main__":
    main()
