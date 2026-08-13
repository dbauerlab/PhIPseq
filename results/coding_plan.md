# PhIP-Seq phage library QC: peptide representation

## Context

You grew 6 batches of two T7 phage-display libraries yourselves (TC_1–TC_4 = T7CoV library, TC_5–TC_6 = T7Vir3.2/VirScan3 library) and sequenced them **before any selection** — these are naive/input pools. The open question is basic QC: after cloning and growing the library ourselves, how much of the designed peptide set actually survived, and how evenly represented is it (dropouts, jackpotting, batch-to-batch consistency)? That's a prerequisite check before using these libraries in real serum-profiling PhIP-seq experiments.

Deliverable: for each of the 6 samples, a per-peptide read count against the *designed* library it belongs to, plus summary QC (coverage, evenness, dropouts, batch reproducibility) — answering "what peptides, at what representation."

This will live in a **new, separate project** (not the airway RNA-seq repo), so it is not bound by that repo's R/Quarto convention. Python is the natural fit here (pandas/openpyxl for the Excel/CSV references, bowtie2/samtools for alignment-based counting, matplotlib for plots) and keeps the pipeline transparent and inspectable rather than pulling in a heavyweight Nextflow framework for a single QC run.

## What we've confirmed so far

- **Raw data**: `/Volumes/lab-bauerd/data/STPs/genomics-stp/inputs/tiphaine.cayol/PM26049/20260810_LH00442_0282_A23G2Y5LT3/fastq/` — 6 samples (`CAY9116A7`…`CAY9116A12`, `S70`–`S75`), each split across lanes L007/L008, paired-end 101bp. Assumed mapping (tentatively confirmed, and self-verifying — see QC step below): `A7→TC_1, A8→TC_2, A9→TC_3, A10→TC_4` (T7CoV), `A11→TC_5, A12→TC_6` (T7Vir3.2).
- **Read architecture** (inspected directly from fastq): R1 = `[~30nt outer vector/primer seq][GAATTCGGAGCGGT][insert, forward strand, partial — read runs out before insert ends]`. R2 = `[~40nt outer vector/primer seq][NotI/HindIII vector junction][insert, reverse-complement, partial from the 3' end]`. The two reads cover the two ends of the insert independently, with a gap in the middle — plenty for unique peptide ID (each fragment is ~50–57nt, versus ~128k peptides max in the larger library).
- **Reference files** (already provided):
  - `.../PhIP-Seq/Elledge Lab Files/References/CoV Library Reference.xlsx` — 6,932 peptides (T7CoV, SARS-CoV-2 + related coronaviruses, 56mer tiles). Has a `Nucleotide sequence` column: `GAATTCGGAGCGGT` + insert + `CACTGCACTCGAGA` — i.e. the exact constant flanks seen in the real reads, confirming this is the right reference.
  - `.../References/VirScan3 .../virscan_annotations/virscan3.peptide.metadata.csv` — 128,257 peptides (T7 Vir3.2/VirScan3, pan-viral). Has an `oligo` column, mixed-case: uppercase = constant flank, lowercase = the insert — self-describing, easy to split.
  - **Gotcha for implementation**: my reconnaissance parse of the `.xlsx` via raw zipfile/XML looked column-shifted (values didn't line up with the header names I listed). Load it properly with `pandas.read_excel` (or `openpyxl` directly) in the real script — don't reuse my manual XML parser.
- **Local machine has no bioinformatics tools at all** (no R, no bowtie2/samtools/cutadapt) but does have Homebrew and Python 3.9 + pip3. Per your choice, we install what's needed locally via Homebrew and process on this laptop against the mounted share.

## Approach

1. **Environment setup**
   - `brew install bowtie2 samtools cutadapt`
   - Python venv in the new project with `pandas openpyxl biopython matplotlib`

2. **Build a combined reference** (one-off script)
   - Parse both reference files, strip constant flanks from each sequence to get the pure insert, and write one combined FASTA with namespaced IDs (e.g. `CoV|YP_009725301.9.BC1`, `Vir3|<id>`) plus one combined metadata table (organism/protein/peptide sequence/position) keyed by the same ID.
   - `bowtie2-build` the combined FASTA into a single index. One shared index (rather than a separate index per sample) is deliberate: it lets us measure how many reads from a "CoV" sample land on Vir3.2 references and vice versa — a built-in sanity check on the assumed TC_# mapping and on cross-contamination, at negligible extra cost.

3. **Per-sample trimming** — `cutadapt` on R1 and R2 independently, anchored 5' adapter = the constant flank sequence, `--discard-untrimmed`. The discard rate itself is a useful QC number (fraction of reads with the expected construct architecture at all).

4. **Alignment & counting** — `bowtie2` end-to-end against the combined index: R1 with `--norc` (forward strand only, matches insert orientation), R2 separately with `--nofw` (reverse strand only). Lanes (L007/L008) are aligned separately and counts summed per peptide, avoiding a wasteful pre-concatenation of the raw files. `samtools view -F4 | cut -f3 | sort | uniq -c` (or `idxstats`) gives per-peptide counts per read (R1-derived and R2-derived kept separate for a concordance check, R1 used as the primary count).

5. **Aggregate & QC report** (Python script + saved PNGs + a short markdown summary)
   - Peptide × sample count matrix, joined to metadata.
   - Per sample: total reads, % correct architecture (from cutadapt), % mapping, % mapping to the *other* library (contamination/mislabelling check), % of the designed library seen at all (coverage), dropout list, read-count distribution, an evenness metric (Gini coefficient / Lorenz curve) and top-N over-represented clones.
   - Batch reproducibility: pairwise per-peptide count correlation across TC_1–TC_4 (CoV batches) and separately TC_5 vs TC_6 (Vir3.2 batches).

6. **Test small, then run full** — run the whole pipeline on a truncated subsample (e.g. first ~200k reads) of one sample first to sanity-check trimming/mapping rates and tune anything before launching the full ~70GB run (background job given the size).

## New project location

`~/Library/CloudStorage/Dropbox-TheFrancisCrick/Tiphaine Cayol/Coding/phipseq-library-qc/` (sibling to the existing `claude_code_demo` repo) — a fresh git repo with roughly:
```
phipseq-library-qc/
  build_reference.py        # parses both reference files -> combined FASTA + metadata CSV
  run_pipeline.sh           # per-sample cutadapt -> bowtie2 -> counts, for all 6 samples
  summarize.py              # builds count matrix + QC report (figures + markdown)
  reference/                # generated FASTA/index/metadata (or .gitignore'd if large)
  results/
```

## Verification

- Subsample test run first; check cutadapt trim rate and bowtie2 mapping rate look sane (expect >90% for both if the adapter/reference sequences are right) before committing to the full run.
- Cross-library mapping numbers double as verification that the TC_# ↔ CAY9116A# ↔ library assignment is correct — if e.g. "TC_5" (assumed Vir3.2) actually maps overwhelmingly to the CoV reference, that tells us the mapping guess was wrong and needs correcting.
- Final sanity check: known VirScan/CoV tiles should show non-zero, roughly comparable representation across the naive library batches — wildly bimodal or all-zero results would flag a construct/adapter mismatch to debug before trusting the counts.
