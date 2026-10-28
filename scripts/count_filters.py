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
# Shared module, not run directly. Imported (as `from count_filters import ...`)
# by summarize.py, completeness_figures.py, overview_figure.py, top_peptides.py
# and differential_enrichment.py; works because they all live in scripts/.
import json  # parses reference/filter_params.json
from pathlib import Path  # object-oriented file paths

import pandas as pd  # DataFrames for the count matrix and reports

# Repo root: this file is scripts/count_filters.py, so go up two levels.
PROJECT_DIR = Path(__file__).resolve().parent.parent

# Load the thresholds from reference/filter_params.json into a dict.
_filter_params = json.loads((PROJECT_DIR / "reference" / "filter_params.json").read_text())
# (A) Detection threshold: reads needed for a peptide to count as "detected"
# in coverage metrics. Only changes how coverage is *reported*; rows are kept.
COVERAGE_MIN_COUNT = _filter_params["coverage_min_count"]
# (C) Row filter: a peptide needs >= MATRIX_FILTER_MIN_COUNT reads...
MATRIX_FILTER_MIN_COUNT = _filter_params["matrix_filter_min_count"]
# ...in at least this many samples. Set to 2 = the number of replicates per
# condition, so a peptide reproducibly seen in one condition's replicates is
# kept, while a one-sample spike (e.g. PCR jackpot) is not.
MATRIX_FILTER_MIN_SAMPLES = _filter_params["matrix_filter_min_samples"]

# Human-readable filter text for logs/figure captions, e.g. ">= 10 reads in >= 2 samples".
FILTER_DESCRIPTION = f">= {MATRIX_FILTER_MIN_COUNT} reads in >= {MATRIX_FILTER_MIN_SAMPLES} samples"
# Same, for the per-group version, e.g. ">= 10 reads in >= 2 of these samples".
# (Two adjacent f-strings inside parentheses are joined into one string.)
GROUP_FILTER_DESCRIPTION = (f">= {MATRIX_FILTER_MIN_COUNT} reads in >= {MATRIX_FILTER_MIN_SAMPLES} "
                            f"of these samples")


# passes_filter(counts)
#   In:  peptides x samples count DataFrame (rows = ref_ids, columns = samples),
#        typically just the columns of one figure group / contrast.
#   Out: boolean Series indexed like counts' rows: True = keep the peptide.
#   Called by: filter_matrix() below, summarize.py, top_peptides.py,
#              differential_enrichment.py.
def passes_filter(counts):
    """Boolean per row of a peptides x samples count frame: does the
    peptide pass the matrix filter within these samples?"""
    # counts >= N gives a True/False frame; .sum(axis=1) counts the Trues per
    # row (samples with enough reads); then compare to the required number.
    return (counts >= MATRIX_FILTER_MIN_COUNT).sum(axis=1) >= MATRIX_FILTER_MIN_SAMPLES


# filter_matrix(mat, lib_of_ref)
#   In:  mat        - full peptides x samples raw count DataFrame (index = ref_id,
#                     columns = tc_id sample names);
#        lib_of_ref - Series mapping ref_id -> library ('CoV' / 'Vir3').
#   Out: (filtered matrix, per-sample read report, per-library peptide report).
#   Called by: summarize.py (writes count_matrix_filtered.csv / filter_report.csv).
def filter_matrix(mat, lib_of_ref):
    """Apply the matrix filter (C). Returns (filtered matrix, per-sample
    report, per-library peptide counts before/after)."""
    keep = passes_filter(mat)  # True/False per peptide, across *all* samples
    filtered = mat.loc[keep]  # boolean indexing: keep only the passing rows

    # Total reads per sample (column sums) before and after filtering.
    raw_reads = mat.sum(axis=0)
    kept_reads = filtered.sum(axis=0)
    # One row per sample. .values strips the index so the columns line up by
    # position (all are in mat.columns order).
    sample_report = pd.DataFrame({
        "tc_id": mat.columns,  # sample name
        "read_pairs_raw": raw_reads.values,  # reads before filtering
        "read_pairs_kept": kept_reads.values,  # reads on surviving peptides
        # Percentage of each sample's reads that sat on dropped peptides.
        # (A sample with 0 raw reads would give NaN here: 0/0.)
        "pct_read_pairs_removed": ((raw_reads - kept_reads) / raw_reads).values * 100,
    })

    # Library label for every peptide in the matrix, in matrix row order.
    # reindex gives NaN for any ref_id not found in lib_of_ref.
    libs = lib_of_ref.reindex(mat.index)
    # value_counts() tallies peptides per library (NaN labels are not counted);
    # the two Series align on library name. fillna(0) covers a library with no
    # surviving peptides, and astype(int) undoes the float conversion NaN forced.
    library_report = pd.DataFrame({
        "peptides_raw": libs.value_counts(),  # peptides per library before
        "peptides_kept": libs[keep].value_counts(),  # peptides per library after
    }).fillna(0).astype(int)
    # Number of peptides dropped per library.
    library_report["peptides_removed"] = library_report["peptides_raw"] - library_report["peptides_kept"]
    library_report.index.name = "library"  # label the index (becomes a CSV column header)
    return filtered, sample_report, library_report
