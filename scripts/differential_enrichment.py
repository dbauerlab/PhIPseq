"""Differential enrichment of peptides vs a reference treatment -- the
"Phage & beads" mock IP by default (reference_treatment in
reference/de_params.json) -- DESeq2-style via PyDESeq2, separately for
each IP condition.

The mock IP goes through every step except the antibody, so comparing
against it controls for phage that bind the beads / plastic: hits reflect
the antibody. (An earlier version compared against the "Vir3 only" input
library; set reference_treatment and output_dir to rerun that way.)

One negative-binomial GLM per contrast (design ~ condition, treatment vs
reference, 2 vs 2 samples), two-sided Wald test and Benjamini-Hochberg
adjusted p-values with DESeq2 independent filtering. Every non-reference
treatment in the IP is tested, including the "Vir3 only" input library,
whose contrast is a QC view (peptides higher in the input than after a
mock IP) rather than an IP result.

- apeGLM LFC shrinkage is off by default (`shrink_lfc`): PyDESeq2's
  implementation is unstable on this data -- for a minority of peptides
  the shrunken value collapsed to ~0.3 while the raw LFC was 13-14, so
  it cannot be trusted for ranking. Raw log2FoldChange is reported and
  plotted; enabling the flag adds log2FoldChange_shrunk as a
  supplementary column.
- `frac_p_below_alpha` in enrichment_summary.csv is a calibration check:
  under the null it should sit near alpha. Much higher means the model
  understates variability and the p-values are anti-conservative -- which
  is the case for IP4 here, whose replicates barely correlate (r =
  0.01-0.16, see the QC batch-correlation figures). Contrasts above
  CALIBRATION_WARN_FRAC get calibration_flag=True in
  enrichment_summary.csv, a stderr warning, and a red volcano title.

- Each contrast is filtered and fitted on *its own* samples
  (count_filters.passes_filter: >= 10 reads in >= 2 of the treatment +
  reference samples). The filter is blind to which sample is in which
  group, which keeps the p-values valid. Fitting all of an IP's samples
  at once instead admitted peptides absent from the reference that hit
  PyDESeq2's optimiser bound (+-30/ln(2) = +-43.28), broke apeGLM and let
  2-8 CPM noise reach padj ~1e-17; fitting per contrast avoids this.
- Peptides with zero reads in every reference replicate are still
  tested -- strong true hits are often absent from the background (the
  12CA5 positive control's top HA peptide has 0 reads in both "Vir3 only"
  samples and in one beads sample) -- but their fold change is a lower
  bound, not an estimate. They are flagged rather than removed (see the
  flags below).
- Size factors: `poscounts` (DESeq2's estimator for sparse data). The
  default median-of-ratios needs peptides counted in every sample, which
  few are in IP4 given its low library coverage.
- Peptides failing a contrast's filter are still reported, with NA
  statistics and filtered=True, so every peptide appears for every
  contrast.
- Cook's-distance outlier handling is off: it needs >= 3 replicates per
  treatment (DESeq2 skips it below that too), and there are 2 here.
- A peptide is "enriched" with padj < alpha and log2FoldChange >
  lfc_threshold (unshrunk; shrunken LFCs are for ranking and plotting),
  "depleted" with padj < alpha and log2FoldChange < -lfc_threshold.
- Flags on enriched peptides, added after testing (they never change the
  calls, they say how much to trust them):
    reference_absent  0 reads in every reference replicate: the LFC is a
                      lower bound ("present vs absent").
    low_count_hit     mean DESeq2-normalised count in the treatment below
                      low_count_max (de_params.json), so both sides are
                      low (e.g. 1 -> 5): likely sampling noise.
    mock_enriched     only when the reference is NOT the mock IP: also
                      enriched in the same IP's mock_treatment contrast,
                      i.e. likely a bead / plastic binder. Off (column
                      absent) when the mock IP is the reference.

Reads <run_dir>/results/count_matrix.csv (raw integer counts),
reference/samples.tsv + samples_meta.csv, and
<run_dir>/reference/combined_metadata.csv. Writes, under
<run_dir>/results/<output_dir>/ (de_params.json; "differential_vs_beads"):
  differential_results.csv   all peptides x all contrasts (long format)
  enrichment_summary.csv     per contrast: n tested / enriched / depleted,
                             flag counts, calibration, top hit
  organism_summary.csv       per contrast x organism: n enriched peptides
  <ip>/<treatment>_vs_<reference>.csv   one contrast, sorted by padj
  <ip>/size_factors.csv      per contrast x sample
  <ip>/*.png                 volcano + MA (top hits labelled), p-value
                             histograms, dispersions, normalisation box
                             plots, enriched-peptide heatmap, enrichment
                             counts, organism dot plot
  <ip>/top_hits/<treatment>_vs_<reference>_top<N>.png
                             normalised counts of each contrast's top hits

Usage: differential_enrichment.py [samples_tsv] [samples_meta_csv]
  (default to reference/samples.tsv and reference/samples_meta.csv)
Threads: SLURM_CPUS_PER_TASK if set, else 1.
"""
# --- standard library ---
import json  # reads the JSON config files in reference/
import os  # os.environ: picks up the SLURM CPU allocation
import re  # regular expressions, used to make treatment names filename-safe
import sys  # sys.stderr for progress messages, sys.exit for fatal errors
from pathlib import Path  # object-oriented file paths ("dir / name" joins paths)

# --- third-party numerics / plotting ---
import numpy as np
import pandas as pd
import matplotlib

# "Agg" is a non-interactive backend that only renders to files: needed on
# cluster nodes with no display. It must be selected before pyplot is imported.
matplotlib.use("Agg")
import matplotlib.colors as mcolors  # colour-name -> RGBA conversion
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # invisible "proxy" markers for hand-built legends
from matplotlib.patches import Patch  # filled-square proxy handles for legends
# PyDESeq2 is a Python re-implementation of the R package DESeq2:
#  - DeseqDataSet holds the counts + sample metadata and fits the model
#    (size factors, dispersions, per-peptide negative-binomial GLM coefficients);
#  - DefaultInference is the numerical back-end (IRLS / optimisers), run in
#    parallel over peptides with joblib;
#  - DeseqStats takes a fitted DeseqDataSet and runs the Wald test for one
#    contrast, then computes adjusted p-values.
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Sibling modules in scripts/ (importable because the script is run as
# `python scripts/differential_enrichment.py`, which puts scripts/ on sys.path):
#  - passes_filter: ">= 10 reads in >= 2 samples" row filter (filter_params.json);
#    GROUP_FILTER_DESCRIPTION is the same rule as human-readable text.
from count_filters import GROUP_FILTER_DESCRIPTION, passes_filter
#  - load_samples: merged, ordered sample sheet + fixed treatment -> colour map;
#    legend_outside: places a legend to the right of an axes;
#    treatment_legend_handles: one colour swatch per treatment.
from sample_groups import legend_outside, load_samples, treatment_legend_handles
#  - draw_peptide_heatmap: clustered peptides x samples heatmap (shared with
#    the QC heatmaps); peptide_label / short_name: compact annotation strings.
from top_peptides import draw_peptide_heatmap, peptide_label, short_name

# Repo root: this file is scripts/differential_enrichment.py, so go up two levels.
PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())
RUN_DIR = Path(_paths["run_dir"])
RESULTS_DIR = RUN_DIR / "results"  # count_matrix.csv is read from here
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"  # per-peptide annotation (ref_id index)

# Analysis parameters from reference/de_params.json (keys must all be present,
# so a missing key fails loudly with a KeyError rather than silently defaulting).
_de = json.loads((PROJECT_DIR / "reference" / "de_params.json").read_text())
# Everything this script writes goes under results/<output_dir>/ (e.g.
# "differential_vs_beads"), so runs against different references can sit
# side by side.
OUT_DIR = RESULTS_DIR / _de["output_dir"]
REFERENCE = _de["reference_treatment"]  # baseline of every contrast ("Phage & beads", the mock IP)
ALPHA = _de["alpha"]  # FDR level: a peptide is significant if padj < ALPHA (0.05)
LFC_THRESHOLD = _de["lfc_threshold"]  # min |log2FC| for a call (1 = 2-fold)
# How DESeq2 normalises for sequencing depth ("poscounts" here, see fit_contrast).
SIZE_FACTORS_FIT_TYPE = _de["size_factors_fit_type"]
SHRINK_LFC = _de["shrink_lfc"]  # whether to also compute apeGLM-shrunken LFCs (false here)
N_LABEL_TOP = _de["n_label_top"]  # top enriched hits labelled on volcano/MA panels and plotted in top_hits/
# low_count_hit threshold: an enriched peptide whose mean normalised count in
# the treatment is below this has low counts on both sides of the comparison.
LOW_COUNT_MAX = _de["low_count_max"]
# The mock IP used for the mock_enriched flag, or null in the JSON (None here)
# to switch the flag off. It is also off when the mock IP is the reference
# itself, since every contrast is then already measured against the mock.
MOCK_TREATMENT = _de["mock_treatment"]
MOCK_ENABLED = MOCK_TREATMENT is not None and MOCK_TREATMENT != REFERENCE
# The post-test flag columns in use; a peptide with none of them set counts as "unflagged".
FLAG_COLS = ["reference_absent", "low_count_hit"] + (["mock_enriched"] if MOCK_ENABLED else [])
# Short tag shown after a peptide's plot label for each flag it carries.
FLAG_TAGS = {"reference_absent": "[ref 0]", "low_count_hit": "[low]", "mock_enriched": "[mock]"}
HEATMAP_MAX_PEPTIDES = _de["heatmap_max_peptides"]  # row cap for the enriched-peptide heatmap
ORGANISM_PLOT_MAX = _de["organism_plot_max_organisms"]  # row cap for the organism dot plot

# CPUs for the PyDESeq2 back-end: SLURM sets this env var (as a string) inside a
# job; default to 1 when run interactively.
N_CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
REF_CODE = "a_ref"  # model level names: formula-safe, and the reference sorts first (treatment coding baseline)
NOT_SIG_COLOR = "lightgray"  # plot colour for non-significant peptides
DEPLETED_COLOR = "dimgray"  # depleted peptides (enriched: treatment colour)
PANEL_NCOLS = 3  # multi-panel figures: at most 3 contrast panels per row
LABELLED_PANEL_WIDTH = 7.5  # inches per volcano / MA panel: the plot plus its label column
# Which LFC column the volcano / MA plots use on their x / y axis.
PLOT_LFC = "log2FoldChange_shrunk" if SHRINK_LFC else "log2FoldChange"
CALIBRATION_WARN_FRAC = 0.10  # fraction of raw p < alpha above which a contrast's p-values look anti-conservative


# Symbols that carry meaning in treatment names are spelled out before the
# generic clean-up below, so e.g. "Serum 2%" and "Serum 2" get different names.
SAFE_NAME_SYMBOLS = {"%": "pct", "&": "and", "+": "plus"}


# safe_name: in = any string (a treatment name such as "Serum 0.125%");
# out = a filename-safe version ("Serum_0.125pct"; "Phage & beads" ->
# "Phage_and_beads"). Used for the per-contrast CSV and figure names.
def safe_name(text):
    for symbol, word in SAFE_NAME_SYMBOLS.items():
        text = text.replace(symbol, word)
    # Replace every run of characters other than letters, digits and "." with a
    # single "_", then trim leading/trailing underscores.
    return re.sub(r"[^A-Za-z0-9.]+", "_", text).strip("_")


# check_safe_names: in = every treatment name in the run. Exits with an error
# if two different names would still map to the same file name (which would
# make one contrast's files silently overwrite the other's). Called at the
# start of main(), before any fitting, so the problem is caught in seconds.
def check_safe_names(treatments):
    by_name = {}  # safe name -> list of original treatment names
    for t in dict.fromkeys(treatments):  # unique, in order
        by_name.setdefault(safe_name(t), []).append(t)
    clashes = {name: ts for name, ts in by_name.items() if len(ts) > 1}
    if clashes:
        sys.exit("ERROR: these treatment names map to the same output file name -- rename one of them in "
                 "samples_meta.csv:\n" + "\n".join(f"  {' / '.join(ts)} -> {name}" for name, ts in clashes.items()))


# fit_contrast: in = the IP's raw count matrix (peptides x samples), that IP's
# sample sheet, and one treatment name. Out = (fitted DeseqDataSet, per-peptide
# results DataFrame, index of peptides that passed the filter and were tested).
# Called from main() once per (IP, non-reference treatment). Each contrast gets
# its own model fitted only on its own treatment + reference samples (see the
# module docstring for why this is per contrast rather than one fit per IP).
def fit_contrast(counts_all, ip_samples, treatment):
    """Fit one DESeq2 model for one contrast (treatment vs REFERENCE),
    on those samples only. Returns (dds, results DataFrame, tested index)."""
    # Sample IDs (tc_id) belonging to this treatment or to the reference
    # (normally 2 + 2 samples). .loc[boolean mask, column] gives a Series.
    cols = ip_samples.loc[ip_samples["treatment"].isin([treatment, REFERENCE]), "tc_id"]
    counts = counts_all[cols]  # peptides x (these samples) slice of raw counts
    # Independent pre-filter on just these samples (>= 10 reads in >= 2 of them);
    # it is blind to group labels, so it does not select on the comparison itself.
    tested = counts.index[passes_filter(counts)]
    # Sample metadata for the design formula: one "condition" column, coded
    # REF_CODE ("a_ref") for reference samples and "t01" for treatment samples.
    # The real treatment names (spaces, "%", "&") are not safe in a formula, and
    # "a_ref" sorts before "t01" so it is the baseline level of the factor.
    # np.where picks element-wise between the two codes; the index must be the
    # sample IDs so PyDESeq2 can match rows of metadata to rows of counts.
    metadata = pd.DataFrame(
        {"condition": np.where(ip_samples.set_index("tc_id").loc[cols, "treatment"] == REFERENCE,
                               REF_CODE, "t01")},
        index=cols.values)
    # One inference object shared by the dataset fit and the stats step below.
    inference = DefaultInference(n_cpus=N_CPUS)
    # Build the dataset. PyDESeq2 (like AnnData) wants samples as rows and
    # peptides as columns, hence the transpose .T, and integer counts.
    #  - design="~condition": log(mean count) = intercept + beta * [is treatment],
    #    so beta is the (natural-log) fold change of treatment over reference.
    #  - size_factors_fit_type="poscounts": per-sample depth normalisation that
    #    takes geometric means over the non-zero counts only, so it works on
    #    sparse data where few peptides are seen in every sample (the default
    #    "ratio" = median-of-ratios needs peptides with no zeros at all).
    #  - refit_cooks=False: do not replace Cook's-distance outlier counts and
    #    refit (pointless with 2 replicates per group).
    #  - quiet=True: suppress PyDESeq2's progress printing.
    dds = DeseqDataSet(counts=counts.loc[tested].T.astype(int), metadata=metadata, design="~condition",
                       size_factors_fit_type=SIZE_FACTORS_FIT_TYPE, refit_cooks=False,
                       inference=inference, quiet=True)
    # Run the DESeq2 fitting pipeline, in order: size factors -> peptide-wise
    # dispersion MLEs -> parametric dispersion-vs-mean trend -> prior variance
    # of dispersions around that trend -> MAP (empirical-Bayes) dispersions
    # shrunk towards the trend -> negative-binomial GLM coefficients (LFCs) by
    # IRLS -> Cook's distances (computed, but unused here).
    dds.deseq2()

    # Set up the statistical test for this contrast: ["factor", "level tested",
    # "baseline level"] = treatment vs reference.
    #  - alpha: the FDR target used by independent filtering when choosing its cutoff.
    #  - cooks_filter=False: do not set p-values to NA for Cook's outliers.
    #  - independent_filter=True: DESeq2's independent filtering -- try a range
    #    of baseMean cutoffs, BH-adjust only the peptides above each, and keep the
    #    cutoff that (roughly) maximises discoveries; peptides below it get padj NA.
    ds = DeseqStats(dds, contrast=["condition", "t01", REF_CODE], alpha=ALPHA, cooks_filter=False,
                    independent_filter=True, inference=inference, quiet=True)
    # summary() runs the two-sided Wald test (z = LFC / SE, normal p-value; null
    # LFC = 0), then the independent-filtering + Benjamini-Hochberg adjustment,
    # and stores the table in ds.results_df: baseMean (mean of normalised
    # counts), log2FoldChange, lfcSE (its standard error), stat (Wald z),
    # pvalue, padj; one row per tested peptide.
    ds.summary()
    res = ds.results_df.copy()  # copy, so a later lfc_shrink() cannot overwrite the raw LFCs
    if SHRINK_LFC:
        # apeGLM shrinkage: re-estimates the treatment coefficient with a
        # heavy-tailed (Cauchy) prior centred on 0, pulling noisy, low-count LFCs
        # towards 0. p-values are unchanged. "condition[T.t01]" is the
        # coefficient name the formula gives the t01-vs-baseline term.
        ds.lfc_shrink(coeff="condition[T.t01]")
        # lfc_shrink rewrote ds.results_df in place; keep its values as extra columns.
        res["log2FoldChange_shrunk"] = ds.results_df["log2FoldChange"]
        res["lfcSE_shrunk"] = ds.results_df["lfcSE"]
    # Count padj < ALPHA (NaN compares False, so untested peptides do not count).
    n_sig = int((res["padj"] < ALPHA).sum())
    # Progress line on stderr (stdout is left clean).
    print(f"  {treatment} vs {REFERENCE}: {len(tested)} tested ({len(cols)} samples), {n_sig} padj < {ALPHA}",
          file=sys.stderr)
    return dds, res, tested


# normalised_counts: in = a fitted DeseqDataSet. Out = DataFrame of DESeq2
# size-factor-normalised counts (raw count / sample size factor), peptides
# (rows) x samples (columns) -- the values the model actually compared.
# PyDESeq2 stores them samples x peptides in dds.layers["normed_counts"]
# (AnnData layout), hence the transpose.
def normalised_counts(dds):
    return pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T


# contrast_table: in = one contrast's PyDESeq2 results (`res`), treatment/IP
# names, the IP's sample sheet, all own-library peptide IDs, the peptides
# actually tested, a CPM matrix, this contrast's DESeq2-normalised counts
# (`norm`, tested peptides only), the IP's raw counts and peptide annotation.
# Out = one row per peptide (tested or not) with statistics, enriched/depleted
# calls, the reference_absent / low_count_hit flags (mock_enriched is filled
# in later by main(), once the IP's mock contrast is available), mean CPMs and
# normalised counts, and annotation, sorted by padj. main() writes it as
# <ip>/<treatment>_vs_<reference>.csv and stacks them into differential_results.csv.
def contrast_table(res, treatment, ip, ip_samples, all_refs, tested_refs, cpm, norm, counts_ip, meta):
    """Full per-peptide table for one contrast: every own-library peptide,
    untested (pre-filtered) ones with NA statistics."""
    trt_ids = ip_samples.loc[ip_samples["treatment"] == treatment, "tc_id"]  # treatment sample IDs
    ref_ids = ip_samples.loc[ip_samples["treatment"] == REFERENCE, "tc_id"]  # reference sample IDs
    # Empty frame indexed by every peptide; the index is named so that
    # reset_index() below turns it into a "ref_id" column.
    out = pd.DataFrame(index=pd.Index(all_refs, name="ref_id"))
    # Assigning a scalar fills the whole column (constant per contrast, but
    # needed once the tables are stacked into one long file).
    out["ip"] = ip
    out["treatment"] = treatment
    out["reference"] = REFERENCE
    out["n_treatment"] = len(trt_ids)  # number of replicates on each side
    out["n_reference"] = len(ref_ids)
    # reindex expands res to every peptide (NaN rows for untested ones), then
    # join attaches its columns by index.
    out = out.join(res.reindex(all_refs))
    out["filtered"] = ~out.index.isin(tested_refs)  # True = removed by the pre-filter, not tested
    # Mean CPM (reads per million, depth-normalised) per group, for an
    # intuitive read-out alongside the model's baseMean.
    out["mean_cpm_treatment"] = cpm.loc[all_refs, trt_ids].mean(axis=1)
    out["mean_cpm_reference"] = cpm.loc[all_refs, ref_ids].mean(axis=1)
    # Mean DESeq2-normalised count per group: the values the model compared
    # (baseMean is their mean over all 4 samples). reindex gives NaN for
    # untested peptides, which have no normalised counts.
    out["mean_norm_treatment"] = norm.reindex(all_refs)[trt_ids].mean(axis=1)
    out["mean_norm_reference"] = norm.reindex(all_refs)[ref_ids].mean(axis=1)
    sig = out["padj"] < ALPHA  # NaN padj -> False
    # Calls combine significance with a minimum effect size, using the raw
    # (unshrunken) LFC.
    out["enriched"] = sig & (out["log2FoldChange"] > LFC_THRESHOLD)
    out["depleted"] = sig & (out["log2FoldChange"] < -LFC_THRESHOLD)
    # --- Flags on enriched peptides (see the module docstring) ---
    # reference_absent: not a single read in any reference replicate, so the
    # fold change has no finite estimate (the optimiser just stops somewhere).
    # (counts_ip[ref_ids] == 0).all(axis=1) is True where every reference column is 0.
    ref_zero = (counts_ip.loc[all_refs, ref_ids] == 0).all(axis=1)
    out["reference_absent"] = out["enriched"] & ref_zero
    # low_count_hit: the treatment side is itself low (and, being enriched, the
    # reference is lower still), e.g. 1 -> 5 normalised counts.
    out["low_count_hit"] = out["enriched"] & (out["mean_norm_treatment"] < LOW_COUNT_MAX)
    if MOCK_ENABLED:
        out["mock_enriched"] = False  # placeholder; set by main() from the mock contrast
    ann = meta.loc[all_refs]  # annotation rows, in the same order as `out`
    # One compact label per peptide ("id | organism | protein start-end").
    out["peptide_label"] = [peptide_label(a) for _, a in ann.iterrows()]
    for col in ["organism", "protein_name", "start", "end", "peptide_aa"]:
        # .values copies positionally (same order), skipping index alignment.
        out[col] = ann[col].values
    # ref_id back to a column; most significant first (pvalue breaks the many
    # BH ties in padj), untested NaN rows last.
    return out.reset_index().sort_values(["padj", "pvalue"], na_position="last")


# calibration_fraction: in = a contrast table. Out = the share of tested
# peptides with raw p < ALPHA. Under the null this should be close to ALPHA;
# much higher means the model understates variability (anti-conservative
# p-values). max(..., 1) avoids division by zero when nothing was tested.
def calibration_fraction(tab):
    return float((tab["pvalue"] < ALPHA).sum() / max(tab["pvalue"].notna().sum(), 1))


# panel_grid: in = number of panels n and the width per panel in inches
# (the volcano / MA plots pass a wider value to fit their label column, see
# label_hits). Out = (figure, flat array of n or more axes) laid out up to
# PANEL_NCOLS wide, unused trailing axes hidden. Used by the volcano, MA,
# p-value histogram and dispersion figures (one panel per contrast).
def panel_grid(n, panel_width=5):
    ncols = min(PANEL_NCOLS, n)  # fewer columns if there are fewer panels
    nrows = int(np.ceil(n / ncols))  # enough rows to hold n panels
    # squeeze=False keeps `axes` a 2-D array even for a 1x1 grid.
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_width * ncols, 4.3 * nrows), squeeze=False)
    axes_flat = axes.flatten()  # 1-D, so it can be zipped against the contrasts
    for ax in axes_flat[n:]:
        ax.axis("off")  # hide the empty cells at the end of the last row
    return fig, axes_flat


# status_colors: in = a contrast table (needs boolean enriched/depleted
# columns) and the treatment's colour. Out = N x 4 RGBA array, one colour per
# row, for scatter(c=...). Used by the volcano and MA plots.
def status_colors(tab, color):
    """Nx4 RGBA: treatment colour if enriched, grey if depleted, pale otherwise.
    (Built as an array, not np.where -- an RGBA tuple would broadcast wrongly.)"""
    # Start with every row light grey: np.tile repeats the 4-value RGBA tuple len(tab) times.
    out = np.tile(mcolors.to_rgba(NOT_SIG_COLOR), (len(tab), 1))
    # Boolean-mask row assignment; .to_numpy(dtype=bool) drops the pandas index.
    out[tab["depleted"].to_numpy(dtype=bool)] = mcolors.to_rgba(DEPLETED_COLOR)
    out[tab["enriched"].to_numpy(dtype=bool)] = mcolors.to_rgba(color)
    return out


# status_legend: in = a figure. Adds one shared legend explaining the
# enriched / depleted / not-significant colours under the panels. Out = none.
# Used by the volcano and MA plots.
def status_legend(fig):
    # Proxy Patch handles (not tied to any plotted artist); "enriched" is shown
    # in neutral grey because each panel uses its own treatment colour.
    fig.legend(handles=[Patch(color="gray", label=f"enriched (padj < {ALPHA}, log2FC > {LFC_THRESHOLD}); "
                                                  "treatment colour"),
                        Patch(color=DEPLETED_COLOR, label=f"depleted (padj < {ALPHA}, log2FC < -{LFC_THRESHOLD})"),
                        Patch(color=NOT_SIG_COLOR, label="not significant"),
                        # hollow circle = the reference_absent flag (see draw_status_scatter)
                        Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor="gray",
                               label="enriched, 0 reference reads (LFC is a lower bound)"),
                        # invisible handle used only to explain the label tags
                        Line2D([], [], linestyle="", label=f"labels: top {N_LABEL_TOP} enriched by padj; "
                                                           + " / ".join(FLAG_TAGS[c] for c in FLAG_COLS)
                                                           + " = flags")],
               # anchored just below the bottom centre of the figure, in two rows
               loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol=3, fontsize=9)


# draw_status_scatter: in = an axes, a contrast table, the x and y values to
# plot (Series aligned with the table) and the treatment colour. Draws every
# peptide coloured by status (status_colors), except enriched peptides with
# the reference_absent flag, which are drawn as hollow circles on top so they
# stand out. Used by the volcano and MA plots. Out = none.
def draw_status_scatter(ax, tab, x, y, color):
    hollow = tab["reference_absent"].to_numpy(dtype=bool)  # numpy mask, one entry per row
    # s = marker size; linewidths=0 removes marker outlines; rasterized=True
    # stores the ~100k points as an image so the PNG renders quickly.
    ax.scatter(x[~hollow], y[~hollow], c=status_colors(tab.loc[~hollow], color), s=4, linewidths=0,
               rasterized=True)
    # facecolors="none" + a coloured edge = hollow circle.
    ax.scatter(x[hollow], y[hollow], facecolors="none", edgecolors=color, s=10, linewidths=0.5, rasterized=True)


# top_enriched: in = a contrast table (already sorted by padj, then pvalue).
# Out = its N_LABEL_TOP most significant enriched (up-regulated) rows -- the
# peptides labelled on the volcano and MA plots and drawn in top_hits/.
def top_enriched(tab):
    return tab.loc[tab["enriched"]].head(N_LABEL_TOP)


# hit_label: in = one contrast-table row. Out = a short label: the numeric part
# of the ref_id ("Vir3|81117" -> "81117"), a shortened organism name, and a
# tag for each flag set on the peptide, e.g. "81117 Influenza A virus [ref 0]".
def hit_label(r):
    tags = [FLAG_TAGS[col] for col in FLAG_COLS if r[col]]  # only the flags in use
    return " ".join([r["ref_id"].split("|")[-1], short_name(r["organism"], 22)] + tags)


# label_hits: in = an axes, the rows to label and their x / y positions (same
# order). Writes the labels as a column just outside the right edge of the
# panel, in rank order from the top, each joined to its point by a thin leader
# line -- so labels never overlap each other or the points, however tightly
# the hits cluster. Out = none. (Panels need extra width: see panel_grid.)
def label_hits(ax, top, xs, ys):
    for k, ((_, r), x, y) in enumerate(zip(top.iterrows(), xs, ys)):
        ax.annotate(f"{k + 1}. {hit_label(r)}",
                    (x, y),  # the point being labelled, in data coordinates
                    # text position in axes fractions (0-1 = inside the panel):
                    # x just right of the panel, y stepping down per label
                    xytext=(1.03, 0.98 - 0.095 * k), textcoords="axes fraction", fontsize=6,
                    # annotation_clip=False: draw the label even though it
                    # lies outside the axes
                    ha="left", va="center", annotation_clip=False,
                    # thin grey leader line from text to point
                    arrowprops=dict(arrowstyle="-", linewidth=0.4, color="gray", alpha=0.7,
                                    shrinkA=0, shrinkB=1))


# plot_volcano: in = {treatment: contrast table} for one IP, the IP name, the
# treatment -> colour map and the output dir. Out = <ip>/volcano.png: one
# panel per contrast, x = log2 fold change, y = -log10(padj), so strong and
# significant hits sit top-left/right. Titles are flagged in red for n=1
# treatments and for contrasts whose p-values look anti-conservative.
def plot_volcano(tables, ip, colors, out_dir):
    lfc = PLOT_LFC  # x-axis column (raw or shrunken LFC)
    fig, axes = panel_grid(len(tables), panel_width=LABELLED_PANEL_WIDTH)
    for ax, (treatment, tab) in zip(axes, tables.items()):
        # Calibration check, over every tested peptide's raw p-value (the same
        # value as enrichment_summary.csv's frac_p_below_alpha).
        frac = calibration_fraction(tab)
        # Only peptides with an adjusted p-value (tested and not removed by
        # independent filtering) can be placed on the y axis.
        tab = tab.loc[tab["padj"].notna()]
        # Clip padj at 1e-300 so a padj of exactly 0 does not give -log10(0) = inf.
        y = -np.log10(tab["padj"].clip(lower=1e-300))
        draw_status_scatter(ax, tab, tab[lfc], y, colors[treatment])
        # Dotted guide lines: the padj = ALPHA significance line...
        ax.axhline(-np.log10(ALPHA), color="black", linestyle=":", linewidth=0.8)
        for x in (-LFC_THRESHOLD, LFC_THRESHOLD):
            # ...and the +-LFC_THRESHOLD effect-size lines.
            ax.axvline(x, color="black", linestyle=":", linewidth=0.8)
        # Label the top N enriched (up-regulated) peptides by padj -- every one,
        # even if several come from the same organism, so recurring organisms show.
        top = top_enriched(tab)
        label_hits(ax, top, top[lfc], y.loc[top.index])  # y.loc[...] picks those rows' heights
        ax.margins(x=0.05)  # small horizontal padding (labels now sit outside the panel)
        # Hit counts among the plotted peptides for the title.
        n_up, n_down = int(tab["enriched"].sum()), int(tab["depleted"].sum())
        flags = []  # warnings to append to the panel title
        if tab["n_treatment"].iloc[0] < 2:
            flags.append("n=1")  # single treatment replicate: very weak test
        if frac > CALIBRATION_WARN_FRAC:
            flags.append(f"{frac:.0%} of raw p < {ALPHA}: anti-conservative")
        # Flags go on an extra title line in square brackets, if any.
        suffix = ("\n[" + "; ".join(flags) + "]") if flags else ""
        # Flagged panels get a red title.
        ax.set_title(f"{treatment} vs {REFERENCE}\n{n_up} enriched, {n_down} depleted{suffix}",
                     fontsize=10, color="firebrick" if flags else "black")
        ax.set_xlabel("log2 fold change" + (" (shrunken)" if SHRINK_LFC else ""))
        ax.set_ylabel("-log10(adjusted p)")
    fig.suptitle(f"Differential enrichment vs {REFERENCE} ({ip})", fontsize=13)
    status_legend(fig)
    # tight_layout within a rectangle (left, bottom, right, top in figure
    # fractions) that leaves room for the legend below and the suptitle above.
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    # bbox_inches="tight" grows the saved image to include the out-of-figure legend.
    fig.savefig(out_dir / "volcano.png", dpi=150, bbox_inches="tight")
    plt.close(fig)  # free the figure's memory (many figures are made in one run)


# plot_ma: same inputs as plot_volcano. Out = <ip>/ma_plots.png: one MA plot
# per contrast (x = baseMean, the mean of size-factor-normalised counts, on a
# log scale; y = log2 fold change). Shows how fold changes and calls depend on
# abundance -- low-count peptides have the noisiest LFCs.
def plot_ma(tables, ip, colors, out_dir):
    lfc = PLOT_LFC  # y-axis column (raw or shrunken LFC)
    fig, axes = panel_grid(len(tables), panel_width=LABELLED_PANEL_WIDTH)
    for ax, (treatment, tab) in zip(axes, tables.items()):
        tab = tab.loc[tab["baseMean"].notna()]  # tested peptides only (untested have NaN stats)
        draw_status_scatter(ax, tab, tab["baseMean"], tab[lfc], colors[treatment])
        # Same top-N enriched peptides as the volcano, at their (baseMean, LFC) positions.
        top = top_enriched(tab)
        label_hits(ax, top, top["baseMean"], top[lfc])
        ax.axhline(0, color="black", linewidth=0.8)  # no-change line
        ax.set_xscale("log")  # baseMean spans several orders of magnitude
        ax.margins(x=0.05)  # small horizontal padding (labels now sit outside the panel)
        ax.set_title(f"{treatment} vs {REFERENCE}", fontsize=10)
        ax.set_xlabel("mean of normalised counts")
        ax.set_ylabel("log2 fold change" + (" (shrunken)" if SHRINK_LFC else ""))
    fig.suptitle(f"MA plots vs {REFERENCE} ({ip})", fontsize=13)
    status_legend(fig)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))  # leave room for legend and suptitle
    fig.savefig(out_dir / "ma_plots.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_pvalue_histograms: same inputs as plot_volcano. Out =
# <ip>/pvalue_histograms.png, one histogram of raw p-values per contrast.
# A well-calibrated test gives a flat (uniform) histogram from the null
# peptides plus a spike near 0 from true hits; a broad excess at low p
# (or a hump) means the model under-estimates variance.
def plot_pvalue_histograms(tables, ip, colors, out_dir):
    fig, axes = panel_grid(len(tables))
    for ax, (treatment, tab) in zip(axes, tables.items()):
        # dropna removes untested peptides; 50 bins of width 0.02 over [0, 1].
        ax.hist(tab["pvalue"].dropna(), bins=50, color=colors[treatment], edgecolor="black", linewidth=0.3)
        ax.set_title(f"{treatment} vs {REFERENCE}", fontsize=10)
        ax.set_xlabel("raw p-value")
        ax.set_ylabel("peptides")
    fig.suptitle(f"p-value distributions ({ip}): expect flat + a spike near 0", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))  # leave room at the top for the suptitle
    fig.savefig(out_dir / "pvalue_histograms.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_dispersions: in = {treatment: fitted DeseqDataSet} for one IP. Out =
# <ip>/dispersions.png, DESeq2's standard dispersion diagnostic per contrast.
# Dispersion (alpha) is the negative-binomial over-dispersion: variance =
# mu + alpha * mu^2. The plot shows peptide-wise estimates (black), the
# fitted mean-dispersion trend (red) and the final shrunken values (blue)
# that the Wald test actually uses.
def plot_dispersions(dds_by_treatment, ip, out_dir):
    """One dispersion panel per contrast (each is its own fit)."""
    fig, axes = panel_grid(len(dds_by_treatment))
    for ax, (treatment, dds) in zip(axes, dds_by_treatment.items()):
        var = dds.var  # per-peptide (AnnData "var") table where PyDESeq2 stores its estimates
        # Mean of size-factor-normalised counts per peptide (a PyDESeq2 internal
        # column, hence the leading underscore).
        means = var["_normed_means"]
        # Peptide-wise maximum-likelihood dispersions: noisy with only 4 samples.
        ax.scatter(means, var["genewise_dispersions"], s=2, color="black", alpha=0.3, linewidths=0,
                   rasterized=True, label="peptide-wise estimate")
        # Final dispersions: MAP estimates shrunk towards the trend (PyDESeq2
        # keeps the peptide-wise value instead for strong high-side outliers).
        ax.scatter(means, var["dispersions"], s=2, color="tab:blue", alpha=0.3, linewidths=0, rasterized=True,
                   label="final (MAP)")
        # Sort by mean so the trend is drawn as a line left to right, not a scribble.
        order = np.argsort(means.values)
        # Fitted trend: dispersion = asymptDisp + extraPois / mean (DESeq2's
        # parametric fit), the prior centre for the MAP step.
        ax.plot(means.values[order], var["fitted_dispersions"].values[order], color="tab:red", linewidth=1.5,
                label="fitted trend")
        ax.set_xscale("log")  # both axes span orders of magnitude
        ax.set_yscale("log")
        ax.set_title(f"{treatment} vs {REFERENCE}", fontsize=10)
        ax.set_xlabel("mean of normalised counts")
        ax.set_ylabel("dispersion")
        ax.legend(markerscale=4, fontsize=7)  # markerscale enlarges the tiny dots in the legend
    fig.suptitle(f"Dispersion estimates ({ip})", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))  # leave room for the suptitle
    fig.savefig(out_dir / "dispersions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_enriched_heatmap: in = {treatment: contrast table}, IP name, IP sample
# sheet, CPM matrix for the IP's samples, annotation, colours, output dir.
# Out = <ip>/enriched_heatmap.png (via top_peptides.draw_peptide_heatmap), or
# nothing if the IP has no enriched peptides. Shows the raw abundance pattern
# of the hits in every sample of the IP, not just the contrast's own four.
def plot_enriched_heatmap(tables, ip, ip_samples, cpm, meta, colors, out_dir):
    """Union of enriched peptides across this IP's contrasts (capped at
    HEATMAP_MAX_PEPTIDES by best padj), log10(CPM + 1) across all of the
    IP's samples."""
    # Stack (ref_id, padj) of the enriched peptides from every contrast.
    best = pd.concat([t.loc[t["enriched"], ["ref_id", "padj"]] for t in tables.values()])
    if best.empty:
        print(f"  no enriched peptides in {ip}: skipping heatmap", file=sys.stderr)
        return
    # A peptide enriched in several contrasts keeps its smallest padj; sort best first.
    best = best.groupby("ref_id")["padj"].min().sort_values()
    n_total = len(best)  # union size, before the cap (reported in the title)
    ref_ids = best.head(HEATMAP_MAX_PEPTIDES).index  # cap the number of rows
    # draw_peptide_heatmap expects a "figure group" dict like the QC scripts use
    # (see sample_groups.figure_groups); build a minimal one for this IP.
    group = {"name": ip, "samples": ip_samples, "fig_dir": out_dir, "title_suffix": f" ({ip})"}
    title = (f"Enriched peptides vs {REFERENCE} ({ip}): {len(ref_ids)} of {n_total} "
             f"(union over contrasts, best padj first; rows clustered)")
    # log10(CPM + 1): compresses the huge dynamic range; +1 keeps zeros finite (log10(1) = 0).
    draw_peptide_heatmap(np.log10(cpm.loc[ref_ids] + 1), meta, group, colors, title,
                         out_dir / "enriched_heatmap.png")


# plot_enrichment_counts: in = the summary table (one row per contrast, all
# IPs), the IP name, colours, output dir. Out = <ip>/enrichment_counts.png: a
# bar chart with enriched counts as bars above zero and depleted counts as
# hatched bars below zero, one bar pair per contrast.
def plot_enrichment_counts(summary, ip, colors, out_dir):
    sub = summary.loc[summary["ip"] == ip]  # this IP's contrasts only
    # Width grows with the number of contrasts.
    fig, ax = plt.subplots(figsize=(0.8 * len(sub) + 3, 4.5))
    x = np.arange(len(sub))  # bar positions 0, 1, 2, ...
    bar_colors = [colors[t] for t in sub["treatment"]]
    ax.bar(x, sub["n_enriched"], color=bar_colors, edgecolor="black", linewidth=0.5)  # upward bars
    # Depleted counts drawn as negative heights: paler (alpha) and hatched.
    ax.bar(x, -sub["n_depleted"], color=bar_colors, edgecolor="black", linewidth=0.5, alpha=0.4, hatch="//")
    for xi, (up, down) in enumerate(zip(sub["n_enriched"], sub["n_depleted"])):
        ax.text(xi, up, f"{up}", ha="center", va="bottom", fontsize=8)  # count above the up bar
        ax.text(xi, -down, f"{down}", ha="center", va="top", fontsize=8)  # count below the down bar
    ax.axhline(0, color="black", linewidth=0.8)  # zero baseline between the two
    ax.set_xticks(x)
    ax.set_xticklabels(sub["treatment"], rotation=45, ha="right")  # angled treatment names
    ax.set_ylabel("peptides  (enriched up, depleted down)")
    ax.set_title(f"Significant peptides vs {REFERENCE} ({ip})\n(padj < {ALPHA}, |log2FC| > {LFC_THRESHOLD})")
    # 15% head-room above and below so the count labels are not clipped
    # (the ", 1" keeps the pad non-zero when every count is 0).
    pad = max(sub["n_enriched"].max(), sub["n_depleted"].max(), 1) * 0.15
    ax.set_ylim(-sub["n_depleted"].max() - pad, sub["n_enriched"].max() + pad)
    fig.savefig(out_dir / "enrichment_counts.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_organisms: in = organism summary (ip x treatment x organism ->
# n_enriched), the IP name, colours, output dir. Out =
# <ip>/organism_enrichment.png: a dot plot with organisms on y, contrasts on
# x and dot area proportional to the number of enriched peptides; limited to
# the ORGANISM_PLOT_MAX organisms with most enriched peptides across contrasts.
# Draws nothing if the IP has no enriched peptides at all.
def plot_organisms(org_summary, ip, colors, out_dir):
    sub = org_summary.loc[org_summary["ip"] == ip]
    if sub.empty:
        return
    # Organisms ranked by total enriched peptides over this IP's contrasts.
    top_orgs = sub.groupby("organism")["n_enriched"].sum().sort_values(ascending=False).head(ORGANISM_PLOT_MAX)
    # Contrasts with any hits, first-appearance order (dict.fromkeys de-duplicates, keeps order).
    treatments = list(dict.fromkeys(sub["treatment"]))
    sub = sub.loc[sub["organism"].isin(top_orgs.index)]  # keep only the top organisms' rows
    # Category -> integer position lookups. y is reversed so the top organism
    # is drawn at the top of the axis.
    y_of = {o: i for i, o in enumerate(top_orgs.index[::-1])}
    x_of = {t: i for i, t in enumerate(treatments)}
    # Figure size scales with the number of columns and rows.
    fig, ax = plt.subplots(figsize=(0.8 * len(treatments) + 6, 0.3 * len(top_orgs) + 2))
    scale = 400 / max(sub["n_enriched"].max(), 1)  # largest dot gets marker area 400 (points^2)
    # One dot per (contrast, organism); `s` is marker area, so area is proportional to count.
    ax.scatter([x_of[t] for t in sub["treatment"]], [y_of[o] for o in sub["organism"]],
               s=sub["n_enriched"] * scale, c=[colors[t] for t in sub["treatment"]], edgecolor="black",
               linewidth=0.5)
    ax.set_xticks(range(len(treatments)))
    ax.set_xticklabels(treatments, rotation=45, ha="right")
    ax.set_yticks(range(len(y_of)))
    # strains are often distinguished only by the parenthesised part, so keep it here
    ax.set_yticklabels([short_name(o, 52, strip_parens=False) for o in top_orgs.index[::-1]], fontsize=8)
    ax.set_xlim(-0.6, len(treatments) - 0.4)  # half a column of padding either side
    ax.grid(color="lightgray", linewidth=0.5)
    ax.set_axisbelow(True)  # grid lines behind the dots
    ax.set_title(f"Enriched peptides by organism vs {REFERENCE} ({ip})\n"
                 f"(top {len(top_orgs)} organisms by total enriched peptides)")
    # Size-legend reference values: median, 90th percentile and max dot count
    # (a set removes duplicates; sorted gives increasing sizes).
    sizes = sorted({int(v) for v in np.quantile(sub["n_enriched"], [0.5, 0.9, 1.0])})
    # Line2D markersize is a diameter in points, scatter `s` an area, hence sqrt.
    legend_outside(ax, [Line2D([], [], marker="o", linestyle="", color="gray", markeredgecolor="black",
                               markersize=np.sqrt(s * scale), label=f"{s} peptides") for s in sizes],
                   fontsize=8, labelspacing=1.5)
    fig.savefig(out_dir / "organism_enrichment.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_normalisation_boxplot: in = the IP's raw counts (own-library peptides x
# the IP's samples), its sample sheet, the IP name, colours and output dir.
# Out = <ip>/normalisation_boxplot.png: two stacked panels of box plots, one
# box per sample coloured by treatment -- top: raw log10(count + 1); bottom:
# log10(normalised count + 1) after DESeq2 size-factor normalisation. If the
# normalisation works, the bottom boxes line up much better than the top ones.
# The size factors are fitted once over all of the IP's samples, with the same
# method (poscounts) as the tests; each contrast's own fit uses only its 4
# samples, so its size factors (in size_factors.csv) differ slightly. Uses
# the peptides passing the filter within the IP, since with every peptide the
# mostly-zero IP4 samples would just show flat boxes at 0.
def plot_normalisation_boxplot(counts_ip, ip_samples, ip, colors, out_dir):
    kept = counts_ip.loc[passes_filter(counts_ip)]  # peptides passing the filter within this IP
    # PyDESeq2 needs a design even just to fit size factors; any per-sample
    # factor will do. pd.factorize numbers the treatments 0, 1, 2, ..., giving
    # formula-safe level names "t00", "t01", ...
    codes = pd.factorize(ip_samples["treatment"])[0]
    metadata = pd.DataFrame({"condition": [f"t{c:02d}" for c in codes]}, index=ip_samples["tc_id"].values)
    dds = DeseqDataSet(counts=kept.T.astype(int), metadata=metadata, design="~condition",
                       size_factors_fit_type=SIZE_FACTORS_FIT_TYPE, quiet=True)
    dds.fit_size_factors()  # only the normalisation step, no model fitting
    norm = normalised_counts(dds)
    size_factors = dds.obs["size_factors"]  # Series indexed by tc_id
    tc_ids = list(ip_samples["tc_id"])  # plotting order (ip -> treatment -> sheet order)
    box_colors = [colors[t] for t in ip_samples["treatment"]]

    # Two rows sharing the x axis; width grows with the number of samples.
    fig, axes = plt.subplots(2, 1, figsize=(0.45 * len(tc_ids) + 4, 9), sharex=True)
    panels = [(axes[0], kept, "raw counts"),
              (axes[1], norm, f"DESeq2-normalised counts ({SIZE_FACTORS_FIT_TYPE} size factors)")]
    for ax, data, what in panels:
        # One array of log10(count + 1) values per sample; +1 keeps zeros finite.
        values = [np.log10(data[t].to_numpy(dtype=float) + 1) for t in tc_ids]
        # patch_artist=True draws filled boxes we can colour; showfliers=False
        # hides the thousands of outlier points beyond the whiskers (1.5 x IQR).
        bp = ax.boxplot(values, patch_artist=True, showfliers=False, widths=0.6)
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
        for median in bp["medians"]:
            median.set_color("black")  # black median line, visible on any fill colour
        ax.set_ylabel("log10(count + 1)")
        ax.set_title(f"{what}, {len(kept):,} peptides passing the filter within {ip}", fontsize=10)
        ax.grid(axis="y", color="lightgray", linewidth=0.5)
        ax.set_axisbelow(True)  # grid behind the boxes
    # boxplot puts box i at x = i + 1; label each with its sample and size factor.
    axes[1].set_xticks(range(1, len(tc_ids) + 1))
    axes[1].set_xticklabels([f"{t}  (SF {size_factors[t]:.2f})" for t in tc_ids], rotation=90, fontsize=8)
    legend_outside(axes[0], treatment_legend_handles(ip_samples, colors), fontsize=8)
    fig.suptitle(f"Normalisation check ({ip}): per-sample count distributions before and after DESeq2 "
                 f"size factors", fontsize=12)
    fig.savefig(out_dir / "normalisation_boxplot.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# plot_top_hits: in = one contrast's table, its DESeq2-normalised counts
# (tested peptides x the contrast's samples), the IP's sample sheet, the
# treatment and IP names, colours and output dir. Out =
# <ip>/top_hits/<treatment>_vs_<reference>_top<N>.png: one small panel per
# top-N enriched peptide (as labelled on the volcano), showing each reference
# and treatment replicate's normalised count -- exactly the numbers DESeq2
# compared -- with a bar at each group mean. Draws nothing if there are no
# enriched peptides.
def plot_top_hits(tab, norm, ip_samples, treatment, ip, colors, out_dir):
    top = top_enriched(tab)
    if top.empty:
        return
    groups = [(REFERENCE, ip_samples.loc[ip_samples["treatment"] == REFERENCE, "tc_id"]),
              (treatment, ip_samples.loc[ip_samples["treatment"] == treatment, "tc_id"])]
    ncols = 5
    nrows = int(np.ceil(len(top) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.4 * nrows), squeeze=False)
    axes = axes.flatten()
    for ax in axes[len(top):]:
        ax.axis("off")  # hide unused panels when there are fewer than N hits
    for rank, (ax, (_, r)) in enumerate(zip(axes, top.iterrows()), start=1):
        for xi, (name, ids) in enumerate(groups):
            vals = norm.loc[r["ref_id"], ids]  # this peptide's normalised count in each replicate
            # Small horizontal spread so replicates with equal values do not overlap.
            jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) > 1 else np.zeros(1)
            ax.scatter(xi + jitter, vals, color=colors[name], edgecolor="black", s=45, zorder=3)
            ax.hlines(vals.mean(), xi - 0.22, xi + 0.22, color="black", linewidth=1.2)  # group mean
        # symlog: linear between 0 and 1, logarithmic above -- so zero counts
        # (common in the reference) still appear while large counts stay readable.
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(bottom=-0.3)  # a little room below 0 so zero points are not cut in half
        ax.set_xlim(-0.6, 1.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([short_name(name, 16) for name, _ in groups], fontsize=8)
        ax.set_title(f"{rank}. {hit_label(r)}\n{short_name(r['protein_name'], 30)} {r['start']}-{r['end']}\n"
                     f"log2FC {r['log2FoldChange']:.1f}, padj {r['padj']:.1e}", fontsize=7)
        if rank % ncols == 1:
            ax.set_ylabel("DESeq2-normalised count")  # only on the first panel of each row
    fig.suptitle(f"Top {len(top)} enriched peptides: {treatment} vs {REFERENCE} ({ip})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))  # leave room for the suptitle
    hits_dir = out_dir / "top_hits"
    hits_dir.mkdir(exist_ok=True)
    fig.savefig(hits_dir / f"{safe_name(treatment)}_vs_{safe_name(REFERENCE)}_top{N_LABEL_TOP}.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)


# main: loads samples, annotation and raw counts; for each IP that has
# reference samples, fits every treatment-vs-reference contrast, adds the
# mock_enriched flag (when in use) and writes its CSVs; then writes the combined tables,
# prints any calibration warning, and draws every figure per IP. Run as a
# script (see differential.sbatch).
def main():
    # Sample sheet (one row per sample: tc_id, library, ip, treatment, ...) in
    # plotting order, and the fixed treatment -> colour mapping.
    samples, colors = load_samples()
    check_safe_names(samples["treatment"])  # exits now if two treatments share a file name
    # Peptide annotation; dtype=str keeps IDs/coordinates as text; indexed by ref_id.
    meta = pd.read_csv(META_CSV, dtype=str).set_index("ref_id")
    # Raw integer counts, peptides (ref_id) x samples (tc_id); first column is the index.
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)
    libraries = samples["library"].unique()
    if len(libraries) != 1:
        # The analysis assumes every sample was made from the same phage library.
        sys.exit(f"ERROR: differential_enrichment.py expects samples from one library, got {', '.join(libraries)}")
    # Peptides that belong to this library and appear in the count matrix
    # (other libraries' peptides would be all-zero noise here).
    own_refs = meta.index[meta["library"] == libraries[0]].intersection(mat.index)
    counts_all = mat.loc[own_refs, samples["tc_id"]]  # own-library peptides x this run's samples
    # CPM: divide each column by its total (own-library reads) and scale to
    # 1e6; used for descriptive mean CPMs and the heatmap, not for the model.
    cpm = counts_all / counts_all.sum(axis=0) * 1e6
    OUT_DIR.mkdir(parents=True, exist_ok=True)  # create results/<output_dir>/ if missing

    all_tables, summary_rows = [], []  # every contrast table; one summary dict per contrast
    ip_tables = {}  # ip -> (tables, ip_samples, dds_by_treatment), kept for plotting later
    # Analyse each IP condition separately (sort=True: IPs in sorted order).
    for ip, ip_samples in samples.groupby("ip", sort=True):
        ip_samples = ip_samples.reset_index(drop=True)  # 0..n-1 row index within this IP
        if REFERENCE not in set(ip_samples["treatment"]):
            # Nothing to compare against in this IP.
            print(f"Skipping {ip}: no '{REFERENCE}' samples", file=sys.stderr)
            continue
        out_dir = OUT_DIR / ip  # per-IP output folder
        out_dir.mkdir(parents=True, exist_ok=True)
        counts_ip = counts_all[ip_samples["tc_id"]]  # this IP's samples only
        # Every non-reference treatment in this IP, de-duplicated in sample-sheet order.
        treatments = [t for t in dict.fromkeys(ip_samples["treatment"]) if t != REFERENCE]
        print(f"{ip}: {len(treatments)} contrasts vs {REFERENCE}, filtered per contrast "
              f"({GROUP_FILTER_DESCRIPTION}), size factors: {SIZE_FACTORS_FIT_TYPE}, {N_CPUS} CPUs",
              file=sys.stderr)

        # Per-IP collectors: contrast tables, fitted models (for the dispersion
        # plot), normalised counts (for the top-hit plots) and size-factor tables.
        tables, dds_by_treatment, norms, size_factor_rows = {}, {}, {}, []
        for treatment in treatments:
            # Fit treatment vs reference on just those samples.
            dds, res, tested = fit_contrast(counts_ip, ip_samples, treatment)
            dds_by_treatment[treatment] = dds
            norms[treatment] = normalised_counts(dds)  # the counts this fit compared
            # Size factors of this fit's samples (obs_names = sample IDs, in model order).
            size_factor_rows.append(pd.DataFrame({"treatment": treatment, "tc_id": dds.obs_names,
                                                  "size_factor": dds.obs["size_factors"].values}))
            # Expand to all own-library peptides, add calls, flags, CPMs and annotation.
            tables[treatment] = contrast_table(res, treatment, ip, ip_samples, own_refs, tested, cpm,
                                               norms[treatment], counts_ip, meta)

        # mock_enriched flag (only when the reference is not the mock IP): needs
        # the IP's mock contrast, so it is set only once every contrast has been
        # fitted. A peptide enriched in the mock IP ("Phage & beads", no
        # antibody) is likely binding the beads, not an antibody.
        if not MOCK_ENABLED:
            pass  # the mock IP is the reference (or no mock set): nothing to flag
        elif MOCK_TREATMENT in tables:
            mock = tables[MOCK_TREATMENT]
            mock_hits = set(mock.loc[mock["enriched"], "ref_id"])
            for treatment, tab in tables.items():
                if treatment != MOCK_TREATMENT:  # the mock contrast cannot flag itself
                    tab["mock_enriched"] = tab["enriched"] & tab["ref_id"].isin(mock_hits)
        else:
            print(f"  no '{MOCK_TREATMENT}' contrast in {ip}: mock_enriched not set", file=sys.stderr)

        for treatment, tab in tables.items():
            # One CSV per contrast, e.g. IP4/12CA5_2ug_vs_Phage_and_beads.csv.
            tab.to_csv(out_dir / f"{safe_name(treatment)}_vs_{safe_name(REFERENCE)}.csv", index=False)
            enriched = tab.loc[tab["enriched"]]  # sorted by padj, so row 0 is the top hit
            # Enriched peptides with none of the flags in use: the most trustworthy calls.
            flagged = enriched[FLAG_COLS].any(axis=1)
            frac = calibration_fraction(tab)
            # One summary row per contrast (int()/float() turn numpy scalars into
            # plain Python numbers for clean CSV output).
            summary_rows.append({
                "ip": ip, "treatment": treatment, "reference": REFERENCE,
                "n_treatment": int(tab["n_treatment"].iloc[0]), "n_reference": int(tab["n_reference"].iloc[0]),
                # peptides that passed the pre-filter and went into the model
                "n_tested": int((~tab["filtered"]).sum()),
                # calibration check over every tested peptide's raw p-value
                # (should be ~ALPHA if most peptides are null)
                "frac_p_below_alpha": frac,
                "calibration_flag": frac > CALIBRATION_WARN_FRAC,
                # tested peptides that survived independent filtering (got a padj)
                "n_with_padj": int(tab["padj"].notna().sum()),
                "n_enriched": int(tab["enriched"].sum()), "n_depleted": int(tab["depleted"].sum()),
                # how many enriched peptides carry each flag (a peptide can carry several)
                "n_enriched_reference_absent": int(enriched["reference_absent"].sum()),
                "n_enriched_low_count": int(enriched["low_count_hit"].sum()),
                # **{...} merges in the mock count only when that flag is in use
                **({"n_enriched_mock": int(enriched["mock_enriched"].sum())} if MOCK_ENABLED else {}),
                "n_enriched_unflagged": int((~flagged).sum()),
                # best enriched peptide (lowest padj), or blank/NaN if there is none
                "top_enriched_peptide": enriched["peptide_label"].iloc[0] if len(enriched) else "",
                "top_enriched_log2FC": enriched["log2FoldChange"].iloc[0] if len(enriched) else np.nan,
                "top_enriched_padj": enriched["padj"].iloc[0] if len(enriched) else np.nan,
            })
        # Stack the per-contrast size-factor tables into <ip>/size_factors.csv.
        pd.concat(size_factor_rows).to_csv(out_dir / "size_factors.csv", index=False)
        all_tables.extend(tables.values())
        ip_tables[ip] = (tables, ip_samples, dds_by_treatment, norms)

    # Per-contrast summary across all IPs -> enrichment_summary.csv.
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "enrichment_summary.csv", index=False)
    # Long table: every peptide x every contrast (ignore_index renumbers rows).
    combined = pd.concat(all_tables, ignore_index=True)
    combined.to_csv(OUT_DIR / "differential_results.csv", index=False)
    # Count enriched peptides per (ip, treatment, organism): keep enriched rows,
    # label missing organisms "unknown", groupby(...).size() counts rows per
    # group; then sort within each contrast by count, largest first.
    org_summary = (combined.loc[combined["enriched"]]
                   .assign(organism=lambda d: d["organism"].fillna("unknown"))
                   .groupby(["ip", "treatment", "organism"]).size().rename("n_enriched").reset_index()
                   .sort_values(["ip", "treatment", "n_enriched"], ascending=[True, True, False]))
    org_summary.to_csv(OUT_DIR / "organism_summary.csv", index=False)
    # Short console overview of the hit counts.
    print(summary[["ip", "treatment", "n_tested", "n_enriched", "n_enriched_unflagged", "n_depleted"]]
          .to_string(index=False), file=sys.stderr)
    # Calibration warning: these contrasts' p-values are anti-conservative,
    # so their hit counts are upper bounds.
    bad = summary.loc[summary["calibration_flag"]]
    if not bad.empty:
        print(f"WARNING: {len(bad)} of {len(summary)} contrasts have > {CALIBRATION_WARN_FRAC:.0%} of raw p < "
              f"{ALPHA} (expected ~{ALPHA:.0%}): p-values anti-conservative, hit counts are upper bounds -- "
              + ", ".join(f"{r.ip} {r.treatment} ({r.frac_p_below_alpha:.0%})" for r in bad.itertuples()),
              file=sys.stderr)

    # Figures, per IP (done after all fits so the summaries above are complete).
    for ip, (tables, ip_samples, dds_by_treatment, norms) in ip_tables.items():
        out_dir = OUT_DIR / ip
        plot_volcano(tables, ip, colors, out_dir)
        plot_ma(tables, ip, colors, out_dir)
        plot_pvalue_histograms(tables, ip, colors, out_dir)
        plot_dispersions(dds_by_treatment, ip, out_dir)
        plot_normalisation_boxplot(counts_all[ip_samples["tc_id"]], ip_samples, ip, colors, out_dir)
        for treatment, tab in tables.items():
            plot_top_hits(tab, norms[treatment], ip_samples, treatment, ip, colors, out_dir)
        # heatmap gets the CPM columns of this IP's samples only
        plot_enriched_heatmap(tables, ip, ip_samples, cpm[ip_samples["tc_id"]], meta, colors, out_dir)
        plot_enrichment_counts(summary, ip, colors, out_dir)
        plot_organisms(org_summary, ip, colors, out_dir)
        print(f"Wrote figures to {out_dir}", file=sys.stderr)
    print(f"Wrote {OUT_DIR / 'differential_results.csv'} ({len(combined)} rows)", file=sys.stderr)


# Run main() only when executed as a script, not when imported as a module.
if __name__ == "__main__":
    main()
