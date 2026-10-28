"""Build a combined bowtie2 reference (FASTA + metadata) from the two designed
PhIP-seq libraries (T7CoV and T7 Vir3.2/VirScan3).

Each source reference stores the full synthesized oligo (constant vector flanks
+ variable peptide-coding insert), with a different masking convention:
  - CoV Library Reference.xlsx: a single 'Nucleotide sequence' column, insert
    flanked by fixed anchors 'GAATTCGGAGCGGT' (5') and 'CACTGCACTCGAGA' (3').
    The unique per-clone ID is 'Barcode ID' (e.g. 'YP_009725301.9.BC1') --
    'Peptide ID' is the *parent* peptide shared by synonymous barcode
    replicates (BC1/BC2/...), a deliberate design to catch PCR jackpotting.
  - virscan3.peptide.metadata.csv: an 'oligo' column, and a 'source' column
    with three provenances (Vir2/Vir3/IEDB) that were masked inconsistently
    when this file was built:
      - source == 'Vir2': mixed case, insert = longest lowercase run
        (verified against peptide length on a sample).
      - source in ('Vir3', 'IEDB'): no case masking at all -- insert is
        found via the same constant anchors 'GGAATTCCGCTGCGT' (5') /
        'GAAGAGCTCGA' (3') confirmed present in 100% of these rows.
    All three provenances share the same vector anchors (confirmed directly
    against the Vir2 rows), i.e. they're one physical synthesized pool --
    the actual T7 Vir3.2 library -- so all are included.

Output: <run_dir>/reference/combined_peptides.fasta (namespaced IDs
'CoV|<Barcode ID>' / 'Vir3|<id>') and
<run_dir>/reference/combined_metadata.csv, where <run_dir> is configured in
reference/paths.json -- generated reference data is kept out of the repo
checkout, separate from the hand-maintained config that lives in
reference/*.json here.

Source-file locations and library-design constants (vector anchors, minimum
insert length) live in reference/paths.json and
reference/build_reference_params.json rather than being hardcoded here, so
redeploying to a new machine only means editing those files.
"""
# --- Imports -----------------------------------------------------------------
# Usage: `python scripts/build_reference.py` (no command-line arguments; all
# settings come from the JSON config files below). Run once per library
# design, before build_index.sbatch turns the FASTA into a bowtie2 index.
import csv  # stdlib CSV reader/writer: reads the Vir3 source, writes the metadata
import json  # parses the reference/*.json config files
import re  # regular expressions, used to find lowercase runs in Vir2 oligos
import sys  # gives us sys.stderr for progress/log messages
from pathlib import Path  # object-oriented file paths (`/` joins path parts)

import pandas as pd  # only used here to read the CoV Excel (.xlsx) file

# Hand-maintained config lives alongside the scripts, in the repo checkout.
# __file__ is this script; .resolve() makes it absolute, .parent is scripts/,
# .parent.parent is the repo root, so CONFIG_DIR = <repo>/reference.
CONFIG_DIR = Path(__file__).resolve().parent.parent / "reference"
# Read each JSON file as text and parse it into a plain Python dict.
# paths.json: source-file locations + run_dir (where generated data goes).
_paths = json.loads((CONFIG_DIR / "paths.json").read_text())
# build_reference_params.json: vector anchor sequences + minimum insert length.
_params = json.loads((CONFIG_DIR / "build_reference_params.json").read_text())

# Generated reference data (FASTA + metadata) is written to run_dir instead
# of the repo checkout, so cloning/redeploying the repo never has to carry
# run output around.
OUT_DIR = Path(_paths["run_dir"]) / "reference"  # <run_dir>/reference/
FASTA_OUT = OUT_DIR / "combined_peptides.fasta"  # insert DNA sequences, input to bowtie2-build
META_OUT = OUT_DIR / "combined_metadata.csv"  # one annotation row per FASTA record

# Absolute paths to the two source library design files (outside the repo).
COV_XLSX = _paths["cov_xlsx"]  # T7CoV design spreadsheet
VIR3_CSV = _paths["vir3_csv"]  # T7 Vir3.2 / VirScan3 design table

# CoV "anchors": constant vector DNA either side of the variable insert.
# Every synthesized oligo is [5' vector flank][peptide-coding insert][3' vector flank];
# only the insert differs between clones, so it is all bowtie2 should align to.
COV_ANCHOR_5 = _params["cov_anchor_5"]  # last bases of the 5' flank, e.g. 'GAATTCGGAGCGGT'
COV_ANCHOR_3 = _params["cov_anchor_3"]  # first bases of the 3' flank, e.g. 'CACTGCACTCGAGA'

# Vir3 library uses a different vector, hence different anchors.
VIR3_ANCHOR_5 = _params["vir3_anchor_5"]  # e.g. 'GGAATTCCGCTGCGT'
VIR3_ANCHOR_3 = _params["vir3_anchor_3"]  # e.g. 'GAAGAGCTCGA'
# Inserts shorter than this (in nt) are treated as parse failures and dropped.
MIN_VIR3_INSERT_LEN = _params["min_vir3_insert_len"]

# Columns (and their order) written to combined_metadata.csv. Each record dict
# also carries "insert_seq", which goes to the FASTA only, not to this CSV.
META_COLS = ["ref_id", "library", "source_id", "parent_peptide_id", "organism", "protein_name", "peptide_aa", "start", "end", "insert_len"]


# load_cov()
#   In:  nothing (reads the CoV .xlsx at COV_XLSX).
#   Out: list of record dicts, one per CoV clone that has both anchors, with
#        keys META_COLS + "insert_seq" (the uppercase insert DNA).
#   Called by: main() below only.
def load_cov():
    # Read the first sheet of the spreadsheet into a DataFrame (needs openpyxl).
    df = pd.read_excel(COV_XLSX)
    records = []  # collected output records
    n_bad = 0  # rows we could not extract an insert from
    # iterrows() yields (index, row) pairs; each row is a pandas Series keyed
    # by column name. Slow for big tables, but the CoV library is small.
    for _, row in df.iterrows():
        # Full synthesized oligo as text; str() guards against non-string cells
        # (an empty cell would become "nan" and fail the anchor search below).
        seq = str(row["Nucleotide sequence"]).strip()
        # Position of the first 5' anchor (-1 if absent). Case-sensitive search.
        i5 = seq.find(COV_ANCHOR_5)
        # Position of the *last* 3' anchor (rfind searches from the right), so
        # an anchor-like motif inside the insert cannot truncate it.
        i3 = seq.rfind(COV_ANCHOR_3)
        # Skip the row if either anchor is missing or they are in the wrong order.
        if i5 == -1 or i3 == -1 or i3 <= i5:
            n_bad += 1
            continue
        # Slice out the DNA strictly between the anchors (start just after the
        # 5' anchor, stop at the start of the 3' anchor) and normalise to uppercase.
        insert = seq[i5 + len(COV_ANCHOR_5): i3].upper()
        # Anchors overlapping or directly adjacent leave an empty insert: skip.
        if not insert:
            n_bad += 1
            continue
        # Unique per-clone ID, e.g. 'YP_009725301.9.BC1' (BCn = barcode replicate).
        barcode_id = str(row["Barcode ID"])
        # One record per clone. row.get(col, "") returns "" if a column is absent.
        records.append({
            # Namespaced ID used as the FASTA header / bowtie2 reference name;
            # the "CoV|" prefix keeps IDs unique across the two libraries.
            "ref_id": f"CoV|{barcode_id}",
            "library": "CoV",  # which library this peptide belongs to
            "source_id": barcode_id,  # the ID as it appears in the source file
            # Parent peptide shared by synonymous barcode replicates (BC1/BC2/...):
            # same amino acids, different codons, so their counts can be compared
            # to spot PCR jackpotting.
            "parent_peptide_id": row.get("Peptide ID", ""),
            # Note: 'Organsim' [sic] is the actual (misspelt) column name in the xlsx.
            "organism": row.get("Organsim", ""),
            "protein_name": row.get("Protein name", ""),
            "peptide_aa": row.get("Peptide sequence", ""),  # encoded amino-acid sequence
            "start": row.get("Start position", ""),  # peptide start within the protein
            "end": row.get("End position", ""),  # peptide end within the protein
            "insert_len": len(insert),  # insert length in nucleotides (3 per amino acid)
            "insert_seq": insert,  # the DNA written to the FASTA (not to the CSV)
        })
    # Progress summary to stderr, so stdout stays clean.
    print(f"[CoV] parsed {len(records)} peptides, skipped {n_bad} without both anchors", file=sys.stderr)
    return records


# Regex for one or more consecutive lowercase DNA bases. E.g. in
# 'GGAATTCCGCTGCGTatgcgtacc...gaGAAGAGCTCGA' findall() returns the lowercase
# stretches such as ['atgcgtacc...ga']. In the Vir2 rows of the Vir3 table the
# vector flanks are UPPERCASE and the variable insert is lowercase.
LOWER_RUN = re.compile(r"[acgt]+")


# load_vir3()
#   In:  nothing (reads the VirScan3 CSV at VIR3_CSV).
#   Out: list of record dicts (keys META_COLS + "insert_seq"), one per usable
#        row, with exact duplicate rows collapsed and clashing IDs renamed.
#   Called by: main() below only.
def load_vir3():
    records = []  # collected output records
    n_bad = 0  # rows with no usable insert (missing anchors or too short)
    n = 0  # total data rows read
    n_exact_dup = 0  # rows identical (same id + same insert) to an earlier row
    n_id_collision = 0  # rows whose id was already used by a *different* insert
    seen = {}  # vid -> insert_seq of first occurrence
    # newline="" is what the csv module expects; errors="replace" swaps any
    # undecodable bytes for a placeholder character instead of crashing.
    with open(VIR3_CSV, newline="", encoding="utf-8", errors="replace") as f:
        # DictReader uses the header line as keys, yielding one dict per row.
        reader = csv.DictReader(f)
        for row in reader:
            n += 1
            # Full oligo; `or ""` also turns a None (short row) into "".
            oligo_raw = row.get("oligo", "") or ""
            # Provenance of this peptide: 'Vir2', 'Vir3' or 'IEDB'.
            src = row.get("source", "")
            if src == "Vir2":
                # Vir2 rows are case-masked: the insert is lowercase, so take
                # every lowercase run...
                runs = LOWER_RUN.findall(oligo_raw)
                # ...and keep the longest one (max by length), uppercased.
                # If there are no lowercase runs at all, the insert is "".
                insert = max(runs, key=len).upper() if runs else ""
            else:
                # Vir3 / IEDB rows have no case masking: uppercase everything
                # and locate the insert between the constant vector anchors.
                oligo = oligo_raw.upper()
                i5 = oligo.find(VIR3_ANCHOR_5)  # first 5' anchor, -1 if absent
                i3 = oligo.rfind(VIR3_ANCHOR_3)  # last 3' anchor, -1 if absent
                # Conditional expression: slice between the anchors when both
                # are found and the 3' anchor starts after the 5' one ends;
                # otherwise "" (which the length check below then rejects).
                insert = oligo[i5 + len(VIR3_ANCHOR_5): i3] if (i5 != -1 and i3 != -1 and i3 > i5 + len(VIR3_ANCHOR_5)) else ""
            # Applies to all three provenances: too-short (or empty) inserts
            # are treated as unparseable and dropped.
            if len(insert) < MIN_VIR3_INSERT_LEN:
                n_bad += 1
                continue
            # The source file's peptide ID (not guaranteed unique in this file).
            vid = row.get("id", "")
            if vid in seen:
                # Same id seen before: either an exact repeat or a real clash.
                if seen[vid] == insert:
                    n_exact_dup += 1  # identical row repeated in the source file, skip silently
                    continue
                n_id_collision += 1
                # Rename, e.g. '1234' -> '1234.dup3'. The number is the running
                # total of collisions in the whole file, not per id.
                vid = f"{vid}.dup{n_id_collision}"  # distinct peptide sharing a non-unique id, disambiguate
            else:
                # First time we see this id: remember its insert for comparison.
                seen[vid] = insert
            records.append({
                "ref_id": f"Vir3|{vid}",  # namespaced FASTA header, e.g. 'Vir3|1234'
                "library": "Vir3",  # Vir2/Vir3/IEDB rows are all one physical Vir3 pool
                "source_id": vid,  # the (possibly .dupN-suffixed) source id
                "parent_peptide_id": "",  # Vir3 has no barcode replicates
                # Column names below are as spelt in virscan3.peptide.metadata.csv.
                "organism": row.get("Organism", ""),
                "protein_name": row.get("Protein.names", ""),
                "peptide_aa": row.get("peptide", ""),  # encoded amino-acid sequence
                "start": row.get("start", ""),  # peptide start within the protein
                "end": row.get("end", ""),  # peptide end within the protein
                "insert_len": len(insert),  # insert length in nucleotides
                "insert_seq": insert,  # the DNA written to the FASTA
            })
    # Summary of what was kept/skipped, to stderr. (n_bad covers both missing
    # anchors and inserts shorter than MIN_VIR3_INSERT_LEN.)
    print(f"[Vir3] parsed {len(records)} peptides out of {n} rows, skipped {n_bad} (no insert), "
          f"{n_exact_dup} exact-duplicate rows collapsed, {n_id_collision} id collisions disambiguated", file=sys.stderr)
    return records


# main()
#   In:  nothing (uses the module-level config above).
#   Out: writes FASTA_OUT and META_OUT; returns None.
#   Called by: the `if __name__ == "__main__"` block (i.e. running this script).
def main():
    # Create <run_dir>/reference/ (and any missing parents); no error if it exists.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # Parse both libraries and concatenate the two lists (CoV first, then Vir3).
    records = load_cov() + load_vir3()

    seen = set()  # ref_ids already written, for a final de-duplication pass
    n_dup = 0  # how many records were dropped as duplicate ref_ids
    # Open both outputs at once; the `with` block closes them automatically.
    # newline="" stops the csv module writing blank lines on some platforms.
    with open(FASTA_OUT, "w") as fasta, open(META_OUT, "w", newline="") as meta_f:
        # DictWriter writes dicts as CSV rows, with columns in META_COLS order.
        writer = csv.DictWriter(meta_f, fieldnames=META_COLS)
        writer.writeheader()  # first line: the column names
        for r in records:
            # Keep only the first record for any ref_id (bowtie2 needs unique names).
            if r["ref_id"] in seen:
                n_dup += 1
                continue
            seen.add(r["ref_id"])
            # FASTA format: a '>' header line with the sequence name, then the
            # sequence on the next line, e.g. '>CoV|YP_009725301.9.BC1\nCCAACT...\n'.
            fasta.write(f">{r['ref_id']}\n{r['insert_seq']}\n")
            # Dict comprehension: keep only the metadata columns (drops insert_seq).
            writer.writerow({k: r[k] for k in META_COLS})

    # Only warn if any duplicates were actually dropped.
    if n_dup:
        print(f"WARNING: {n_dup} duplicate ref_ids skipped", file=sys.stderr)
    # len(seen) == number of unique records written to both files.
    print(f"Wrote {len(seen)} references to {FASTA_OUT} and {META_OUT}", file=sys.stderr)


# Run main() only when executed as a script, not when imported as a module.
if __name__ == "__main__":
    main()
