"""Aggregate per-sample counts into a peptide x sample matrix and a QC report:
coverage, dropouts, evenness, cross-library mapping, batch reproducibility.

Reads <run_dir>/results/<TC_ID>/counts.tsv (one paired-end fragment count
per peptide, produced by run_sample.sh's merged-lane, paired-mode
cutadapt/bowtie2 pipeline) + its cutadapt.log/bowtie2.log,
<run_dir>/reference/combined_metadata.csv, and (optionally) a per-sample
metadata sheet (reference/samples_meta.csv by default) used to colour/shape
the PCA plot. Writes, all under <run_dir>:
  results/count_matrix.csv
  results/qc_summary.csv
  results/figures/*.png
  results/summary.md
<run_dir> is configured in reference/paths.json, separate from this repo
checkout.

Usage: summarize.py [samples_tsv] [samples_meta_csv]
  samples_tsv       defaults to reference/samples.tsv
  samples_meta_csv  defaults to reference/samples_meta.csv (skipped if missing)
"""
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PROJECT_DIR = Path(__file__).resolve().parent.parent
SAMPLES_TSV = PROJECT_DIR / "reference" / "samples.tsv"
SAMPLES_META_CSV = PROJECT_DIR / "reference" / "samples_meta.csv"

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
LIBRARY_COLORS = _fig_params["library_colors"]
LEGEND_BBOX_ANCHOR = tuple(_fig_params["legend_bbox_anchor"])
XTICK_ROTATION = _fig_params["xtick_rotation"]
HEATMAP_TICK_FONTSIZE = _fig_params["heatmap_tick_fontsize"]
PCA_N_TOP_VARIABLE_PEPTIDES = _fig_params["pca_n_top_variable_peptides"]

MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]


def legend_outside(ax, **kwargs):
    """Place a legend to the right of the axes rather than inside the
    plotting area."""
    return ax.legend(loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR, borderaxespad=0, **kwargs)


def load_samples():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    return pd.read_csv(path, sep="\t")


def load_sample_meta():
    """Optional per-sample metadata (e.g. ip, treatment) used only to
    colour/shape the PCA plot. Returns None if not available."""
    path = Path(sys.argv[2]) if len(sys.argv) > 2 else SAMPLES_META_CSV
    if not path.exists():
        return None
    return pd.read_csv(path, dtype=str)


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


def pick_color_shape_columns(sample_meta_cols):
    """From metadata columns (excluding the sample id column), use the
    lowest-cardinality column for point shape (matplotlib only has a
    handful of practical marker styles) and the highest-cardinality column
    for point colour (a colormap handles many categories better)."""
    candidates = [(sample_meta_cols[c].nunique(dropna=True), c)
                  for c in sample_meta_cols.columns if sample_meta_cols[c].nunique(dropna=True) > 1]
    candidates.sort()
    if not candidates:
        return None, None
    if len(candidates) == 1:
        return candidates[0][1], None
    return candidates[-1][1], candidates[0][1]


def compute_pca(mat, tc_ids):
    """PCA of samples, restricted to peptides detected (>0 reads) in every
    sample and the top N most variable of those (by variance of log1p
    counts), mean-centered per peptide. Returns (scores, explained_var) or
    (None, None) if there isn't enough shared signal to run PCA on."""
    sub = mat[tc_ids]
    detected_in_all = (sub > 0).all(axis=1)
    sub = sub.loc[detected_in_all]
    if sub.shape[0] < 2 or sub.shape[1] < 2:
        return None, None
    log_sub = np.log1p(sub)
    top_ids = log_sub.var(axis=1).sort_values(ascending=False).head(PCA_N_TOP_VARIABLE_PEPTIDES).index
    log_sub = log_sub.loc[top_ids]

    X = log_sub.T.to_numpy()  # samples x peptides
    X = X - X.mean(axis=0, keepdims=True)
    U, S, _Vt = np.linalg.svd(X, full_matrices=False)
    scores = U * S
    explained_var = (S ** 2) / np.sum(S ** 2)
    return scores, explained_var


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    samples = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)
    meta["ref_id"] = meta["ref_id"].astype(str)

    sample_meta = load_sample_meta()
    if sample_meta is not None:
        samples = samples.merge(sample_meta, left_on="tc_id", right_on="sample", how="left")

    print("Building count matrix...", file=sys.stderr)
    mat = build_matrix(samples, meta)
    mat.to_csv(RESULTS_DIR / "count_matrix.csv")

    lib_of_ref = meta.set_index("ref_id")["library"]

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
            "median_count_detected": own_counts[own_counts > 0].median(),
            "gini_own_library": gini(own_counts.values),
            "shannon_entropy": entropy,
            "normalized_shannon_entropy": entropy / np.log2(n_designed) if n_designed > 1 else float("nan"),
        })

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(RESULTS_DIR / "qc_summary.csv", index=False)
    print(qc.to_string(index=False), file=sys.stderr)

    # --- Figures ---
    # Rank-abundance (log-log) per sample, restricted to each sample's own library
    fig, ax = plt.subplots(figsize=(7, 5))
    for _, srow in samples.iterrows():
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index
        vals = np.sort(mat[tc_id].loc[own_refs].values)[::-1]
        vals = vals[vals > 0]
        ax.plot(np.arange(1, len(vals) + 1), vals, label=f"{tc_id} ({library})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Peptide rank")
    ax.set_ylabel("Read-pair count")
    ax.set_title("Rank-abundance (own-library peptides)")
    legend_outside(ax, fontsize=8)
    fig.savefig(FIG_DIR / "rank_abundance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Lorenz curves
    fig, ax = plt.subplots(figsize=(6, 6))
    for _, srow in samples.iterrows():
        tc_id, library = srow["tc_id"], srow["library"]
        own_refs = lib_of_ref[lib_of_ref == library].index
        frac_pop, frac_reads = lorenz_curve(mat[tc_id].loc[own_refs].values)
        ax.plot(frac_pop, frac_reads, label=tc_id)
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="perfect evenness")
    ax.set_xlabel("Cumulative fraction of peptides")
    ax.set_ylabel("Cumulative fraction of reads")
    ax.set_title("Lorenz curves (library evenness)")
    legend_outside(ax, fontsize=8)
    fig.savefig(FIG_DIR / "lorenz_curves.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Batch reproducibility heatmaps, per library
    for library, tc_ids in samples.groupby("library")["tc_id"].apply(list).items():
        own_refs = lib_of_ref[lib_of_ref == library].index
        sub = np.log1p(mat.loc[own_refs, tc_ids])
        corr = sub.corr(method="pearson")
        fig, ax = plt.subplots(figsize=(1.2 * len(tc_ids) + 2, 1.2 * len(tc_ids) + 2))
        im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(tc_ids)))
        ax.set_xticklabels(tc_ids, rotation=45, ha="right", fontsize=HEATMAP_TICK_FONTSIZE)
        ax.set_yticks(range(len(tc_ids)))
        ax.set_yticklabels(tc_ids, fontsize=HEATMAP_TICK_FONTSIZE)
        for i in range(len(tc_ids)):
            for j in range(len(tc_ids)):
                ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)
        ax.set_title(f"{library} batch correlation (log1p read-pair counts)")
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"batch_correlation_{library}.png", dpi=150)
        plt.close(fig)

    # Shannon entropy per sample
    fig, ax = plt.subplots(figsize=(1.2 * len(samples) + 2, 4.5))
    tc_ids_ordered = qc["tc_id"].tolist()
    bar_colors = [LIBRARY_COLORS.get(lib, "gray") for lib in qc["library"]]
    ax.bar(tc_ids_ordered, qc["normalized_shannon_entropy"], color=bar_colors)
    ax.set_xticks(range(len(tc_ids_ordered)))
    ax.set_xticklabels(tc_ids_ordered, rotation=XTICK_ROTATION, ha="center")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Normalized Shannon entropy\n(H / log2(designed library size))")
    ax.set_title("Library diversity: how evenly are reads spread across peptides?")
    handles = [Patch(color=c, label=lib) for lib, c in LIBRARY_COLORS.items() if lib in set(qc["library"])]
    if handles:
        legend_outside(ax, handles=handles, fontsize=8)
    fig.savefig(FIG_DIR / "shannon_entropy.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # PCA of samples
    pca_written = False
    scores, explained_var = compute_pca(mat, samples["tc_id"].tolist())
    if scores is not None:
        color_col, shape_col = (None, None)
        if sample_meta is not None:
            color_col, shape_col = pick_color_shape_columns(sample_meta.drop(columns=["sample"]))

        fig, ax = plt.subplots(figsize=(7.5, 6))

        if color_col:
            categories = samples[color_col].fillna("unknown").astype(str)
            uniq_colors = sorted(categories.unique())
            cmap = plt.get_cmap("tab20", max(len(uniq_colors), 1))
            color_map = {c: cmap(i) for i, c in enumerate(uniq_colors)}
            point_colors = categories.map(color_map)
        else:
            categories = samples["library"].astype(str)
            point_colors = categories.map(lambda lib: LIBRARY_COLORS.get(lib, "gray"))

        if shape_col:
            shape_categories = samples[shape_col].fillna("unknown").astype(str)
            uniq_shapes = sorted(shape_categories.unique())
            shape_map = {c: MARKER_CYCLE[i % len(MARKER_CYCLE)] for i, c in enumerate(uniq_shapes)}
            markers = shape_categories.map(shape_map)
        else:
            markers = pd.Series(["o"] * len(samples))

        for marker in markers.unique():
            idx = np.where(markers.values == marker)[0]
            ax.scatter(scores[idx, 0], scores[idx, 1], c=[point_colors.iloc[i] for i in idx],
                       marker=marker, s=70, edgecolor="black", linewidth=0.5)

        ax.set_xlabel(f"PC1 ({explained_var[0]:.1%} variance)")
        ax.set_ylabel(f"PC2 ({explained_var[1]:.1%} variance)" if len(explained_var) > 1 else "PC2")
        title_suffix = f" (top {PCA_N_TOP_VARIABLE_PEPTIDES} variable, shared-detected peptides)"
        ax.set_title("PCA of samples" + title_suffix, fontsize=10)

        legend_handles = []
        color_label = color_col or "library"
        legend_handles.append(Line2D([], [], linestyle="", label=f"— {color_label} —"))
        for cat, color in (color_map.items() if color_col else LIBRARY_COLORS.items()):
            legend_handles.append(Line2D([], [], marker="o", linestyle="", color=color, markeredgecolor="black",
                                          markersize=8, label=str(cat)))
        if shape_col:
            legend_handles.append(Line2D([], [], linestyle="", label=f"— {shape_col} —"))
            for cat, marker in shape_map.items():
                legend_handles.append(Line2D([], [], marker=marker, linestyle="", color="gray",
                                              markeredgecolor="black", markersize=8, label=str(cat)))
        ax.legend(handles=legend_handles, loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR,
                  borderaxespad=0, fontsize=8)
        fig.savefig(FIG_DIR / "pca_samples.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        pca_written = True
    else:
        print("Skipping PCA figure: fewer than 2 peptides detected in every sample, or fewer than 2 samples.",
              file=sys.stderr)

    # --- Markdown summary ---
    lines = ["# PhIP-seq library QC summary", ""]
    lines.append("## Per-sample QC")
    lines.append("")
    lines.append(qc.to_markdown(index=False, floatfmt=".3f"))
    lines.append("")
    lines.append("## Figures")
    lines.append("")
    lines.append("![Rank abundance](figures/rank_abundance.png)")
    lines.append("")
    lines.append("![Lorenz curves](figures/lorenz_curves.png)")
    lines.append("")
    for library in samples["library"].unique():
        lines.append(f"![{library} batch correlation](figures/batch_correlation_{library}.png)")
        lines.append("")
    lines.append("![Shannon entropy](figures/shannon_entropy.png)")
    lines.append("")
    if pca_written:
        lines.append("![PCA of samples](figures/pca_samples.png)")
        lines.append("")
    (RESULTS_DIR / "summary.md").write_text("\n".join(lines))
    print(f"Wrote {RESULTS_DIR / 'summary.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
