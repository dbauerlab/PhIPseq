"""One combined figure comparing all samples together: coverage (%% of
designed peptides detected) and evenness (typical range + near-dropout
outliers) side by side, colored by library so the two T7CoV vs T7Vir3.2
batches are visually grouped.

Usage: overview_figure.py [samples_tsv]  (defaults to scripts/samples.tsv,
i.e. all 6 samples once the full run has completed)
"""
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
OUTLIER_THRESHOLD = 10  # reads


LIBRARY_COLORS = {"CoV": "tab:blue", "Vir3": "tab:orange"}


def main():
    samples_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    samples = pd.read_csv(samples_path, sep="\t")
    meta = pd.read_csv(META_CSV, dtype=str)
    lib_of_ref = meta.set_index("ref_id")["library"]
    mat = pd.read_csv(RESULTS_DIR / "count_matrix_R1.csv", index_col=0)

    fig, (ax_cov, ax_even) = plt.subplots(1, 2, figsize=(1.6 * len(samples) + 4, 5.5))

    xtick_labels = []
    seen_libs = []
    for i, (_, srow) in enumerate(samples.iterrows()):
        tc_id, library = srow["tc_id"], srow["library"]
        color = LIBRARY_COLORS.get(library, "gray")
        own_refs = lib_of_ref[lib_of_ref == library].index.intersection(mat.index)
        counts = mat.loc[own_refs, tc_id]
        total = len(own_refs)
        detected = counts[counts > 0]
        pct = 100 * len(detected) / total if total else 0

        # --- coverage panel ---
        ax_cov.bar(i, pct, color=color, label=library if library not in seen_libs else None)
        ax_cov.text(i, pct + 1.5, f"{len(detected)}/{total}", ha="center", va="bottom", fontsize=8)

        # --- evenness panel ---
        if len(detected):
            p10, p50, p90 = np.percentile(detected, [10, 50, 90])
            ax_even.plot([i, i], [p10, p90], color=color, linewidth=6, solid_capstyle="butt", alpha=0.85,
                         label=library if library not in seen_libs else None)
            ax_even.plot(i, p50, "o", color="white", markeredgecolor="black", markersize=7, zorder=5)
            outliers = detected[detected <= OUTLIER_THRESHOLD]
            for val in outliers.values:
                ax_even.plot(i, val, "v", color="firebrick", markersize=8, zorder=6,
                             label=f"near-dropout (<={OUTLIER_THRESHOLD} reads)" if "outlier" not in seen_libs else None)
                seen_libs.append("outlier")

        seen_libs.append(library)
        n_missing = int((counts == 0).sum())
        n_outlier = int((detected <= OUTLIER_THRESHOLD).sum()) if len(detected) else 0
        extra_bits = [b for b in [f"{n_outlier} near-dropout" if n_outlier else "", f"{n_missing} missing" if n_missing else ""] if b]
        xtick_labels.append(tc_id + ("\n" + ", ".join(extra_bits) if extra_bits else ""))

    ax_cov.set_ylim(0, 108)
    ax_cov.set_xticks(range(len(samples)))
    ax_cov.set_xticklabels(samples["tc_id"])
    ax_cov.set_ylabel("% of designed peptides detected")
    ax_cov.set_title("Coverage")
    ax_cov.legend(fontsize=8, loc="lower right")

    ax_even.set_yscale("log")
    ax_even.set_xlim(-0.6, len(samples) - 0.4)
    ax_even.set_xticks(range(len(samples)))
    ax_even.set_xticklabels(xtick_labels, fontsize=8)
    ax_even.set_ylabel("Reads per peptide (log scale)")
    ax_even.set_title("Evenness (10th-90th percentile + outliers)")
    ax_even.legend(fontsize=8, loc="upper right")

    fig.suptitle("All samples: library completeness and evenness overview", fontsize=13)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "overview_all_samples.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIG_DIR / 'overview_all_samples.png'}", file=sys.stderr)


if __name__ == "__main__":
    main()
