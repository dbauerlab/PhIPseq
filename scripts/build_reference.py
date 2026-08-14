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

Output: reference/combined_peptides.fasta (namespaced IDs 'CoV|<Barcode ID>'
/ 'Vir3|<id>') and reference/combined_metadata.csv.

Source-file locations and library-design constants (vector anchors, minimum
insert length) live in reference/paths.json and
reference/build_reference_params.json rather than being hardcoded here, so
redeploying to a new machine only means editing those files.
"""
import csv
import json
import re
import sys
from pathlib import Path

import pandas as pd

OUT_DIR = Path(__file__).resolve().parent.parent / "reference"
FASTA_OUT = OUT_DIR / "combined_peptides.fasta"
META_OUT = OUT_DIR / "combined_metadata.csv"

_paths = json.loads((OUT_DIR / "paths.json").read_text())
_params = json.loads((OUT_DIR / "build_reference_params.json").read_text())

COV_XLSX = _paths["cov_xlsx"]
VIR3_CSV = _paths["vir3_csv"]

COV_ANCHOR_5 = _params["cov_anchor_5"]
COV_ANCHOR_3 = _params["cov_anchor_3"]

VIR3_ANCHOR_5 = _params["vir3_anchor_5"]
VIR3_ANCHOR_3 = _params["vir3_anchor_3"]
MIN_VIR3_INSERT_LEN = _params["min_vir3_insert_len"]

META_COLS = ["ref_id", "library", "source_id", "parent_peptide_id", "organism", "protein_name", "peptide_aa", "start", "end", "insert_len"]


def load_cov():
    df = pd.read_excel(COV_XLSX)
    records = []
    n_bad = 0
    for _, row in df.iterrows():
        seq = str(row["Nucleotide sequence"]).strip()
        i5 = seq.find(COV_ANCHOR_5)
        i3 = seq.rfind(COV_ANCHOR_3)
        if i5 == -1 or i3 == -1 or i3 <= i5:
            n_bad += 1
            continue
        insert = seq[i5 + len(COV_ANCHOR_5): i3].upper()
        if not insert:
            n_bad += 1
            continue
        barcode_id = str(row["Barcode ID"])
        records.append({
            "ref_id": f"CoV|{barcode_id}",
            "library": "CoV",
            "source_id": barcode_id,
            "parent_peptide_id": row.get("Peptide ID", ""),
            "organism": row.get("Organsim", ""),
            "protein_name": row.get("Protein name", ""),
            "peptide_aa": row.get("Peptide sequence", ""),
            "start": row.get("Start position", ""),
            "end": row.get("End position", ""),
            "insert_len": len(insert),
            "insert_seq": insert,
        })
    print(f"[CoV] parsed {len(records)} peptides, skipped {n_bad} without both anchors", file=sys.stderr)
    return records


LOWER_RUN = re.compile(r"[acgt]+")


def load_vir3():
    records = []
    n_bad = 0
    n = 0
    n_exact_dup = 0
    n_id_collision = 0
    seen = {}  # vid -> insert_seq of first occurrence
    with open(VIR3_CSV, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n += 1
            oligo_raw = row.get("oligo", "") or ""
            src = row.get("source", "")
            if src == "Vir2":
                runs = LOWER_RUN.findall(oligo_raw)
                insert = max(runs, key=len).upper() if runs else ""
            else:
                oligo = oligo_raw.upper()
                i5 = oligo.find(VIR3_ANCHOR_5)
                i3 = oligo.rfind(VIR3_ANCHOR_3)
                insert = oligo[i5 + len(VIR3_ANCHOR_5): i3] if (i5 != -1 and i3 != -1 and i3 > i5 + len(VIR3_ANCHOR_5)) else ""
            if len(insert) < MIN_VIR3_INSERT_LEN:
                n_bad += 1
                continue
            vid = row.get("id", "")
            if vid in seen:
                if seen[vid] == insert:
                    n_exact_dup += 1  # identical row repeated in the source file, skip silently
                    continue
                n_id_collision += 1
                vid = f"{vid}.dup{n_id_collision}"  # distinct peptide sharing a non-unique id, disambiguate
            else:
                seen[vid] = insert
            records.append({
                "ref_id": f"Vir3|{vid}",
                "library": "Vir3",
                "source_id": vid,
                "parent_peptide_id": "",
                "organism": row.get("Organism", ""),
                "protein_name": row.get("Protein.names", ""),
                "peptide_aa": row.get("peptide", ""),
                "start": row.get("start", ""),
                "end": row.get("end", ""),
                "insert_len": len(insert),
                "insert_seq": insert,
            })
    print(f"[Vir3] parsed {len(records)} peptides out of {n} rows, skipped {n_bad} (no insert), "
          f"{n_exact_dup} exact-duplicate rows collapsed, {n_id_collision} id collisions disambiguated", file=sys.stderr)
    return records


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = load_cov() + load_vir3()

    seen = set()
    n_dup = 0
    with open(FASTA_OUT, "w") as fasta, open(META_OUT, "w", newline="") as meta_f:
        writer = csv.DictWriter(meta_f, fieldnames=META_COLS)
        writer.writeheader()
        for r in records:
            if r["ref_id"] in seen:
                n_dup += 1
                continue
            seen.add(r["ref_id"])
            fasta.write(f">{r['ref_id']}\n{r['insert_seq']}\n")
            writer.writerow({k: r[k] for k in META_COLS})

    if n_dup:
        print(f"WARNING: {n_dup} duplicate ref_ids skipped", file=sys.stderr)
    print(f"Wrote {len(seen)} references to {FASTA_OUT} and {META_OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
