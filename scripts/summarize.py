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
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from count_filters import (COVERAGE_MIN_COUNT, FILTER_DESCRIPTION, GROUP_FILTER_DESCRIPTION,
                           MATRIX_FILTER_MIN_COUNT, MATRIX_FILTER_MIN_SAMPLES, filter_matrix, passes_filter)
from sample_groups import (LINESTYLE_CYCLE, MARKER_CYCLE, figure_groups, legend_outside, load_samples,
                           sample_label, style_legend_handles, style_map, treatment_legend_handles)

PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())
RUN_DIR = Path(_paths["run_dir"])
RESULTS_DIR = RUN_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

# Figure styling shared across all figure-generating scripts.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
XTICK_ROTATION = _fig_params["xtick_rotation"]
HEATMAP_TICK_FONTSIZE = _fig_params["heatmap_tick_fontsize"]
PCA_N_TOP_VARIABLE_PEPTIDES = _fig_params["pca_n_top_variable_peptides"]


def load_counts(tc_id):
    path = RESULTS_DIR / tc_id / "counts.tsv"
    if not path.exists():
        return {}
    counts = {}
    with open(path) as f:
        for line in f:
            ref_id, n = line.rstrip("\n").split("\t")
            counts[ref_id] = int(n)
    return counts


def build_matrix(samples, meta):
    """One peptide x sample matrix of paired-end fragment counts (one
    count per aligned read pair -- see run_sample.sh's counting logic)."""
    cols = {}
    for tc_id in samples["tc_id"]:
        cols[tc_id] = load_counts(tc_id)
    df = pd.DataFrame(index=meta["ref_id"], columns=list(cols.keys()), dtype="int64")
    for tc_id, counts in cols.items():
        df[tc_id] = df.index.map(lambda r: counts.get(r, 0)).astype("int64")
    return df


def parse_cutadapt_paired_log(path):
    """Parse a paired-mode cutadapt summary (pair-level fields, not the
    single-end 'Total reads processed'/'Reads written' field names)."""
    text = path.read_text()
    total = int(re.search(r"Total read pairs processed:\s*([\d,]+)", text).group(1).replace(",", ""))
    written = int(re.search(r"Pairs written \(passing filters\):\s*([\d,]+)", text).group(1).replace(",", ""))
    return total, written


def parse_bowtie2_paired_log(path):
    """Parse a paired-mode bowtie2 summary: concordant-pair alignment
    counts plus the overall (mate-level) alignment rate it reports on its
    last line."""
    text = path.read_text()
    total_pairs = int(re.search(r"(\d+) reads; of these:", text).group(1))
    concordant_0 = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly 0 times", text).group(1))
    concordant_1 = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", text).group(1))
    concordant_multi = int(re.search(r"(\d+) \([\d.]+%\) aligned concordantly >1 times", text).group(1))
    overall_pct = float(re.search(r"([\d.]+)% overall alignment rate", text).group(1))
    return total_pairs, concordant_0, concordant_1, concordant_multi, overall_pct


def sample_readstats(tc_id):
    """cutadapt/bowtie2 pair-level QC stats for one sample. One log file
    per tool per sample now (lanes already merged and R1/R2 already
    trimmed/aligned together upstream in run_sample.sh), so there's nothing
    left to sum across lanes/mates the way this used to."""
    stats = {
        "raw_read_pairs": float("nan"),
        "trim_rate": float("nan"),
        "concordant_map_rate": float("nan"),
        "concordant_multimap_rate": float("nan"),
        "overall_map_rate": float("nan"),
    }
    cpath = RESULTS_DIR / tc_id / "cutadapt.log"
    if cpath.exists():
        raw, written = parse_cutadapt_paired_log(cpath)
        stats["raw_read_pairs"] = raw
        stats["trim_rate"] = written / raw if raw else float("nan")
    bpath = RESULTS_DIR / tc_id / "bowtie2.log"
    if bpath.exists():
        total_pairs, _concordant_0, concordant_1, concordant_multi, overall_pct = parse_bowtie2_paired_log(bpath)
        stats["concordant_map_rate"] = (concordant_1 + concordant_multi) / total_pairs if total_pairs else float("nan")
        stats["concordant_multimap_rate"] = concordant_multi / total_pairs if total_pairs else float("nan")
        stats["overall_map_rate"] = overall_pct / 100.0
    return stats


def gini(counts):
    x = np.sort(np.asarray(counts, dtype=float))
    x = x[x >= 0]
    n = len(x)
    if n == 0 or x.sum() == 0:
        return float("nan")
    cum = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cum) / cum[-1]) / n


def lorenz_curve(counts):
    x = np.sort(np.asarray(counts, dtype=float))
    cum = np.cumsum(x)
    cum = np.insert(cum, 0, 0) / cum[-1]
    frac_pop = np.linspace(0, 1, len(cum))
    return frac_pop, cum


def shannon_entropy(counts, base=2):
    """Shannon diversity (in bits, base 2) of a sample's read counts across
    its own peptide library: how evenly are reads spread across peptides,
    vs. concentrated on a few?"""
    x = np.asarray(counts, dtype=float)
    total = x.sum()
    if total <= 0:
        return float("nan")
    p = x[x > 0] / total
    return float(-np.sum(p * np.log(p)) / np.log(base))


def compute_pca(mat, own_refs, tc_ids):
    """PCA of just these samples (one figure group). Peptides are the
    group's own library peptides passing the matrix filter *within the
    group* (count_filters.passes_filter on these samples; zeros allowed),
    in log10(CPM + 1) -- CPM scaled by
    each sample's raw own-library read total, since depth varies several-
    fold -- then the top N most variable, mean-centered per peptide.
    Returns (scores, explained_var, n_peptides_eligible), or None if there
    isn't enough signal to run PCA on."""
    sub = mat.loc[own_refs, tc_ids]
    cpm = sub / sub.sum(axis=0) * 1e6
    keep = passes_filter(sub)
    if keep.sum() < 2 or len(tc_ids) < 3:
        return None
    log_cpm = np.log10(cpm.loc[keep] + 1)
    top_ids = log_cpm.var(axis=1).sort_values(ascending=False).head(PCA_N_TOP_VARIABLE_PEPTIDES).index

    X = log_cpm.loc[top_ids].T.to_numpy()  # samples x peptides
    X = X - X.mean(axis=0, keepdims=True)
    U, S, _Vt = np.linalg.svd(X, full_matrices=False)
    scores = U * S
    explained_var = (S ** 2) / np.sum(S ** 2)
    return scores, explained_var, int(keep.sum())

def own_library_counts(mat, lib_of_ref, srow):
    """One sample's counts over its own library's peptides."""
    own_refs = lib_of_ref[lib_of_ref == srow["library"]].index
    return mat[srow["tc_id"]].loc[own_refs]


def line_legend(ax, group, colors, linestyles):
    legend_outside(ax, treatment_legend_handles(group["samples"], colors, kind="line")
                   + style_legend_handles(group, linestyles, kind="line"), fontsize=8)


def plot_rank_abundance(mat, lib_of_ref, group, colors):
    """Rank-abundance (log-log) per sample, restricted to each sample's own library."""
    linestyles = style_map(group, LINESTYLE_CYCLE)
    fig, ax = plt.subplots(figsize=(7, 5))
    for _, srow in group["samples"].iterrows():
        vals = np.sort(own_library_counts(mat, lib_of_ref, srow).values)[::-1]
        vals = vals[vals > 0]
        ax.plot(np.arange(1, len(vals) + 1), vals, color=colors[srow["treatment"]],
                linestyle=linestyles[srow[group["style_col"]]], linewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Peptide rank")
    ax.set_ylabel("Read-pair count")
    ax.set_title("Rank-abundance (own-library peptides)" + group["title_suffix"])
    line_legend(ax, group, colors, linestyles)
    fig.savefig(group["fig_dir"] / "rank_abundance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_lorenz(mat, lib_of_ref, group, colors):
    linestyles = style_map(group, LINESTYLE_CYCLE)
    fig, ax = plt.subplots(figsize=(6, 6))
    for _, srow in group["samples"].iterrows():
        frac_pop, frac_reads = lorenz_curve(own_library_counts(mat, lib_of_ref, srow).values)
        ax.plot(frac_pop, frac_reads, color=colors[srow["treatment"]],
                linestyle=linestyles[srow[group["style_col"]]], linewidth=1.2)
    ax.plot([0, 1], [0, 1], "k--", linewidth=1)
    ax.text(0.55, 0.5, "perfect evenness", rotation=45, fontsize=8, ha="center", va="bottom")
    ax.set_xlabel("Cumulative fraction of peptides")
    ax.set_ylabel("Cumulative fraction of reads")
    ax.set_title("Lorenz curves (library evenness)" + group["title_suffix"])
    line_legend(ax, group, colors, linestyles)
    fig.savefig(group["fig_dir"] / "lorenz_curves.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_batch_correlation(mat, lib_of_ref, group, colors):
    """Pearson correlation heatmap of log1p counts, per library, over the
    peptides passing the matrix filter within this group's samples, with
    treatment colour strips along the top and left edges."""
    for library, sub_samples in group["samples"].groupby("library", sort=False):
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        tc_ids = sub_samples["tc_id"].tolist()
        labels = [sample_label(srow, group) for _, srow in sub_samples.iterrows()]
        n = len(tc_ids)
        sub = mat.loc[own_refs, tc_ids]
        sub = sub.loc[passes_filter(sub)]
        corr = np.log1p(sub).corr(method="pearson")

        fig, ax = plt.subplots(figsize=(1.2 * n + 2, 1.2 * n + 2))
        im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="viridis")
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)

        divider = make_axes_locatable(ax)
        ax_top = divider.append_axes("top", size=0.3, pad=0.05, sharex=ax)
        ax_left = divider.append_axes("left", size=0.3, pad=0.05, sharey=ax)
        cax = divider.append_axes("right", size=0.3, pad=0.15)
        strip = np.array([mcolors.to_rgba(colors[t]) for t in sub_samples["treatment"]])
        ax_top.imshow(strip[np.newaxis, :, :], aspect="auto")
        ax_left.imshow(strip[:, np.newaxis, :], aspect="auto")

        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=HEATMAP_TICK_FONTSIZE)
        ax.set_yticks(range(n))
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax_left.set_yticklabels(labels, fontsize=HEATMAP_TICK_FONTSIZE)
        ax_left.set_xticks([])
        ax_top.set_yticks([])
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_top.set_title(f"{library} batch correlation{group['title_suffix']}\n"
                         f"(log1p read-pair counts; {len(sub)} peptides with {GROUP_FILTER_DESCRIPTION})")
        fig.colorbar(im, cax=cax)
        # anchored to the colorbar (in its axes units) so it sits just past the colorbar's tick labels
        cax.legend(handles=treatment_legend_handles(sub_samples, colors), loc="center left",
                   bbox_to_anchor=(3.5, 0.5), borderaxespad=0, fontsize=HEATMAP_TICK_FONTSIZE)
        fig.savefig(group["fig_dir"] / f"batch_correlation_{library}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_shannon(qc, group, colors):
    group_qc = qc.set_index("tc_id").loc[group["samples"]["tc_id"]]
    labels = [sample_label(srow, group) for _, srow in group["samples"].iterrows()]
    fig, ax = plt.subplots(figsize=(1.2 * len(labels) + 2, 4.5))
    ax.bar(range(len(labels)), group_qc["normalized_shannon_entropy"],
           color=[colors[t] for t in group["samples"]["treatment"]], edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=XTICK_ROTATION, ha="center")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Normalized Shannon entropy\n(H / log2(designed library size))")
    ax.set_title("Library diversity: how evenly are reads spread across peptides?" + group["title_suffix"])
    legend_outside(ax, treatment_legend_handles(group["samples"], colors), fontsize=8)
    fig.savefig(group["fig_dir"] / "shannon_entropy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_pca_scores(scores, explained_var, group, n_pcs=5):
    """Per-sample PC scores (first n_pcs) for this group, with each PC's
    variance explained in the column name."""
    samples = group["samples"]
    n_pcs = min(n_pcs, scores.shape[1])
    out = samples[["tc_id", "ip", "treatment", "replicate"]].copy()
    for k in range(n_pcs):
        out[f"PC{k + 1} ({explained_var[k]:.1%})"] = scores[:, k]
    out.to_csv(group["fig_dir"] / "pca_scores.csv", index=False)


def plot_pca(scores, explained_var, n_eligible, group, colors):
    """PCA computed on this group's samples only (see compute_pca)."""
    samples = group["samples"]
    markers = style_map(group, MARKER_CYCLE)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for value, marker in markers.items():
        sel = np.where(samples[group["style_col"]].values == value)[0]
        ax.scatter(scores[sel, 0], scores[sel, 1], c=[colors[samples["treatment"].iloc[i]] for i in sel],
                   marker=marker, s=70, edgecolor="black", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({explained_var[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({explained_var[1]:.1%} variance)" if len(explained_var) > 1 else "PC2")
    which = "all samples" if group["name"] == "all" else f"{group['name']} samples only"
    ax.set_title(f"PCA of samples{group['title_suffix']}\n(computed on {which}: log10(CPM + 1), top "
                 f"{min(PCA_N_TOP_VARIABLE_PEPTIDES, n_eligible)} most variable of {n_eligible} peptides\n"
                 f"with {GROUP_FILTER_DESCRIPTION})",
                 fontsize=10)
    legend_outside(ax, treatment_legend_handles(group["samples"], colors, kind="marker")
                   + style_legend_handles(group, markers, kind="marker"), fontsize=8)
    fig.savefig(group["fig_dir"] / "pca_samples.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)
    meta["ref_id"] = meta["ref_id"].astype(str)

    print("Building count matrix...", file=sys.stderr)
    mat = build_matrix(samples, meta)
    mat.to_csv(RESULTS_DIR / "count_matrix.csv")

    lib_of_ref = meta.set_index("ref_id")["library"]

    # The all-sample filtered matrix is written for sharing (and feeds
    # top_peptides.py's tables); the between-sample figures apply the same
    # filter within each figure group instead (passes_filter); library-
    # completeness metrics below use the raw matrix.
    mat_filtered, filter_report, filter_library_report = filter_matrix(mat, lib_of_ref)
    mat_filtered.to_csv(RESULTS_DIR / "count_matrix_filtered.csv")
    filter_report.to_csv(RESULTS_DIR / "filter_report.csv", index=False)
    print(f"Matrix filter ({FILTER_DESCRIPTION}):\n{filter_library_report.to_string()}", file=sys.stderr)

    qc_rows = []
    for _, srow in samples.iterrows():
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index
        other_refs = lib_of_ref[lib_of_ref != library].index

        col = mat[tc_id]
        own_counts = col.loc[own_refs]
        other_counts = col.loc[other_refs]

        stats = sample_readstats(tc_id)

        total_mapped = own_counts.sum() + other_counts.sum()
        n_designed = len(own_refs)
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
            "pct_mapped_to_other_library": other_counts.sum() / total_mapped if total_mapped else float("nan"),
            "n_designed_peptides": n_designed,
            "n_peptides_detected": int((own_counts > 0).sum()),
            "pct_library_covered": (own_counts > 0).mean(),
            "n_dropouts": int((own_counts == 0).sum()),
            f"n_peptides_detected_min{COVERAGE_MIN_COUNT}": int((own_counts >= COVERAGE_MIN_COUNT).sum()),
            f"pct_library_covered_min{COVERAGE_MIN_COUNT}": (own_counts >= COVERAGE_MIN_COUNT).mean(),
            "median_count_detected": own_counts[own_counts > 0].median(),
            "gini_own_library": gini(own_counts.values),
            "shannon_entropy": entropy,
            "normalized_shannon_entropy": entropy / np.log2(n_designed) if n_designed > 1 else float("nan"),
        })

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(RESULTS_DIR / "qc_summary.csv", index=False)
    print(qc.to_string(index=False), file=sys.stderr)

    # --- Figures: all samples, then one set per IP condition ---
    groups = figure_groups(samples, FIG_DIR)
    libraries = samples["library"].unique()
    pca_written = len(libraries) == 1
    if not pca_written:
        print(f"Skipping PCA: samples span several libraries ({', '.join(libraries)}).", file=sys.stderr)

    for group in groups:
        plot_rank_abundance(mat, lib_of_ref, group, colors)
        plot_lorenz(mat, lib_of_ref, group, colors)
        plot_batch_correlation(mat, lib_of_ref, group, colors)
        plot_shannon(qc, group, colors)
        if pca_written:
            own_refs = lib_of_ref[lib_of_ref == libraries[0]].index.intersection(mat.index)
            pca = compute_pca(mat, own_refs, group["samples"]["tc_id"].tolist())
            if pca is None:
                print(f"Skipping PCA{group['title_suffix']}: too few samples or peptides.", file=sys.stderr)
            else:
                scores, explained_var, n_eligible = pca
                plot_pca(scores, explained_var, n_eligible, group, colors)
                write_pca_scores(scores, explained_var, group)
                print(f"PCA{group['title_suffix']}: {n_eligible} eligible peptides, PC1 {explained_var[0]:.1%}, "
                      f"PC2 {explained_var[1]:.1%}", file=sys.stderr)

    # --- Markdown summary ---
    lines = ["# PhIP-seq library QC summary", ""]
    lines.append("## Per-sample QC")
    lines.append("")
    lines.append(qc.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    lines.append("Top most abundant peptides per sample: [top_abundant_peptides.csv](top_abundant_peptides.csv)")
    lines.append("")
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
    lines.append(filter_library_report.reset_index().to_markdown(index=False))
    lines.append("")
    figure_names = ["rank_abundance", "lorenz_curves"]
    figure_names += [f"batch_correlation_{lib}" for lib in samples["library"].unique()]
    figure_names += ["shannon_entropy"] + (["pca_samples"] if pca_written else [])
    figure_names += ["coverage_bar", "rarefaction_curves", f"rarefaction_curves_min{COVERAGE_MIN_COUNT}",
                     "detection_histogram", "evenness_summary",
                     "overview_all_samples", "top_variable_heatmap"]
    for group in groups:
        lines.append(f"## Figures{group['title_suffix']}")
        lines.append("")
        rel_dir = group["fig_dir"].relative_to(RESULTS_DIR).as_posix()
        for name in figure_names:
            lines.append(f"![{name}]({rel_dir}/{name}.png)")
            lines.append("")
    (RESULTS_DIR / "summary.md").write_text("\n".join(lines))
    print(f"Wrote {RESULTS_DIR / 'summary.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
