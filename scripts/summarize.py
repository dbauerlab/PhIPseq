"""Aggregate per-sample counts into a peptide x sample matrix and a QC report:
coverage, dropouts, evenness, cross-library mapping, batch reproducibility.

Reads results/<TC_ID>/counts.tsv (one paired-end fragment count per
peptide, produced by run_sample.sh's merged-lane, paired-mode
cutadapt/bowtie2 pipeline) + its cutadapt.log/bowtie2.log, and
reference/combined_metadata.csv. Writes:
  results/count_matrix.csv
  results/qc_summary.csv
  results/figures/*.png
  results/summary.md
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
SAMPLES_TSV = PROJECT_DIR / "scripts" / "samples.tsv"
META_CSV = PROJECT_DIR / "reference" / "combined_metadata.csv"


def load_samples():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    return pd.read_csv(path, sep="\t")


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


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    samples = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str)
    meta["ref_id"] = meta["ref_id"].astype(str)

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
        qc_rows.append({
            "tc_id": tc_id,
            "library": library,
            "raw_read_pairs": stats["raw_read_pairs"],
            "trim_rate": stats["trim_rate"],
            "concordant_map_rate": stats["concordant_map_rate"],
            "concordant_multimap_rate": stats["concordant_multimap_rate"],
            "overall_map_rate": stats["overall_map_rate"],
            "pct_mapped_to_other_library": other_counts.sum() / total_mapped if total_mapped else float("nan"),
            "n_designed_peptides": len(own_refs),
            "n_peptides_detected": int((own_counts > 0).sum()),
            "pct_library_covered": (own_counts > 0).mean(),
            "n_dropouts": int((own_counts == 0).sum()),
            "median_count_detected": own_counts[own_counts > 0].median(),
            "gini_own_library": gini(own_counts.values),
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
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "rank_abundance.png", dpi=150)
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
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "lorenz_curves.png", dpi=150)
    plt.close(fig)

    # Batch reproducibility heatmaps, per library
    for library, tc_ids in samples.groupby("library")["tc_id"].apply(list).items():
        own_refs = lib_of_ref[lib_of_ref == library].index
        sub = np.log1p(mat.loc[own_refs, tc_ids])
        corr = sub.corr(method="pearson")
        fig, ax = plt.subplots(figsize=(1.2 * len(tc_ids) + 2, 1.2 * len(tc_ids) + 2))
        im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(tc_ids)))
        ax.set_xticklabels(tc_ids, rotation=45, ha="right")
        ax.set_yticks(range(len(tc_ids)))
        ax.set_yticklabels(tc_ids)
        for i in range(len(tc_ids)):
            for j in range(len(tc_ids)):
                ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)
        ax.set_title(f"{library} batch correlation (log1p read-pair counts)")
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"batch_correlation_{library}.png", dpi=150)
        plt.close(fig)

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
    (RESULTS_DIR / "summary.md").write_text("\n".join(lines))
    print(f"Wrote {RESULTS_DIR / 'summary.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
