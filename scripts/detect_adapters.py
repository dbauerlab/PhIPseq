"""Determine the constant 5' adapter (vector sequence before the peptide
insert) for R1 and R2 of each library, from consensus base-calling over many
reads of one representative sample per library. Written to
reference/adapters.json for use by run_pipeline.sh's cutadapt calls.
"""
import gzip
import json
from collections import Counter
from pathlib import Path

FASTQ_DIR = "/Volumes/lab-bauerd/data/STPs/genomics-stp/inputs/tiphaine.cayol/PM26049/20260810_LH00442_0282_A23G2Y5LT3/fastq"
OUT = Path(__file__).resolve().parent.parent / "reference" / "adapters.json"

REPRESENTATIVE_SAMPLES = {
    "CoV": "CAY9116A7_S70_L007",
    "Vir3": "CAY9116A11_S74_L007",
}

WIDTH = 70
N_READS = 5000
DROP_THRESHOLD = 0.8  # consensus fraction below this marks the adapter/insert boundary


def consensus_prefix(fastq_path, n_reads=N_READS, width=WIDTH):
    counters = [Counter() for _ in range(width)]
    with gzip.open(fastq_path, "rt") as f:
        i = 0
        for lineno, line in enumerate(f):
            if lineno % 4 == 1:
                seq = line.strip()
                for pos in range(min(width, len(seq))):
                    counters[pos][seq[pos]] += 1
                i += 1
                if i >= n_reads:
                    break
    bases, fracs = [], []
    for c in counters:
        base, count = c.most_common(1)[0]
        bases.append(base)
        fracs.append(count / sum(c.values()))
    return bases, fracs


def find_adapter(fastq_path):
    bases, fracs = consensus_prefix(fastq_path)
    boundary = next((i for i, f in enumerate(fracs) if f < DROP_THRESHOLD), len(bases))
    return "".join(bases[:boundary])


def main():
    adapters = {}
    for lib, sample_prefix in REPRESENTATIVE_SAMPLES.items():
        r1 = find_adapter(f"{FASTQ_DIR}/{sample_prefix}_R1_001.fastq.gz")
        r2 = find_adapter(f"{FASTQ_DIR}/{sample_prefix}_R2_001.fastq.gz")
        adapters[lib] = {"R1": r1, "R2": r2}
        print(f"{lib}: R1 adapter ({len(r1)}nt) = {r1}")
        print(f"{lib}: R2 adapter ({len(r2)}nt) = {r2}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(adapters, indent=2))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
