"""Low-count filtering, in two deliberately separate forms (thresholds in
reference/filter_params.json):

  A. Coverage detection threshold (`coverage_min_count`): for library-
     completeness metrics a peptide counts as "detected" in a sample only
     with >= this many reads. Reported alongside the >=1-read metrics; the
     count matrix itself is untouched, and the denominator stays the full
     designed library.

  C. Matrix filter (`matrix_filter_min_count`, `matrix_filter_min_samples`):
     keep a peptide (row) only if it has >= min_count reads in >=
     min_samples samples -- min_samples matching the replicate count, so a
     peptide consistently seen in one condition's replicates survives
     without needing group labels. Whole rows are dropped, so there are
     no NAs. Applied in two ways:
       - across all samples -> count_matrix_filtered.csv (the shareable
         filtered matrix; also used for the top-abundant-peptide tables);
       - within each figure group (all samples / one IP), via
         passes_filter() on just that group's samples, for the
         between-sample figures (PCA, batch correlation, top-variable
         heatmap) -- so e.g. the IP4 figures only use peptides seen in
         IP4 samples, not ones kept because of IP5.
     Never used for coverage: the dropped rows are by definition the
     undetected ones, so coverage over the kept rows would be inflated.
"""
import json
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent

_filter_params = json.loads((PROJECT_DIR / "reference" / "filter_params.json").read_text())
COVERAGE_MIN_COUNT = _filter_params["coverage_min_count"]
MATRIX_FILTER_MIN_COUNT = _filter_params["matrix_filter_min_count"]
MATRIX_FILTER_MIN_SAMPLES = _filter_params["matrix_filter_min_samples"]

FILTER_DESCRIPTION = f">= {MATRIX_FILTER_MIN_COUNT} reads in >= {MATRIX_FILTER_MIN_SAMPLES} samples"
GROUP_FILTER_DESCRIPTION = (f">= {MATRIX_FILTER_MIN_COUNT} reads in >= {MATRIX_FILTER_MIN_SAMPLES} "
                            f"of these samples")


def passes_filter(counts):
    """Boolean per row of a peptides x samples count frame: does the
    peptide pass the matrix filter within these samples?"""
    return (counts >= MATRIX_FILTER_MIN_COUNT).sum(axis=1) >= MATRIX_FILTER_MIN_SAMPLES


def filter_matrix(mat, lib_of_ref):
    """Apply the matrix filter (C). Returns (filtered matrix, per-sample
    report, per-library peptide counts before/after)."""
    keep = passes_filter(mat)
    filtered = mat.loc[keep]

    raw_reads = mat.sum(axis=0)
    kept_reads = filtered.sum(axis=0)
    sample_report = pd.DataFrame({
        "tc_id": mat.columns,
        "read_pairs_raw": raw_reads.values,
        "read_pairs_kept": kept_reads.values,
        "pct_read_pairs_removed": ((raw_reads - kept_reads) / raw_reads).values * 100,
    })

    libs = lib_of_ref.reindex(mat.index)
    library_report = pd.DataFrame({
        "peptides_raw": libs.value_counts(),
        "peptides_kept": libs[keep].value_counts(),
    }).fillna(0).astype(int)
    library_report["peptides_removed"] = library_report["peptides_raw"] - library_report["peptides_kept"]
    library_report.index.name = "library"
    return filtered, sample_report, library_report
