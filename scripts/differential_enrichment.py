"""Differential enrichment of peptides vs the background library
("Vir3 only" by default, reference/de_params.json), DESeq2-style via
PyDESeq2, separately for each IP condition.

One negative-binomial GLM per contrast (design ~ condition, treatment vs
reference, 2 vs 2 samples), two-sided Wald test and Benjamini-Hochberg
adjusted p-values with DESeq2 independent filtering. Every non-reference
treatment in the IP is tested, including the "Phage & beads" mock IP.

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
  CALIBRATION_WARN_FRAC are flagged in summary.md and in the figures.

- Each contrast is filtered and fitted on *its own* samples
  (count_filters.passes_filter: >= 10 reads in >= 2 of the treatment +
  reference samples). Fitting all of an IP's samples at once shares
  dispersion across 13-18 samples, but then admits peptides that are
  entirely absent from the 2 reference samples (they pass the per-IP
  filter on the strength of other treatments). Those have no finite MLE
  -- PyDESeq2 returns its optimiser bound, +-30/ln(2) = +-43.28 -- which
  both breaks apeGLM's prior variance (every shrunken LFC collapsed to
  ~0.2, regardless of effect size) and let 2-8 CPM noise reach padj
  ~1e-17, outranking the real signal. Fitting per contrast avoids this;
  dispersion is still stabilised by empirical Bayes across ~100k
  peptides, which matters more here than the sample count.
- Size factors: `poscounts` (DESeq2's estimator for sparse data). The
  default median-of-ratios needs peptides counted in every sample -- only
  47 in IP4, given its low coverage.
- Peptides failing a contrast's filter are still reported, with NA
  statistics and filtered=True, so every peptide appears for every
  contrast.
- Cook's-distance outlier handling is off: it needs >= 3 replicates per
  treatment (DESeq2 skips it below that too), and there are 2 here.
- A peptide is "enriched" with padj < alpha and log2FoldChange >
  lfc_threshold (unshrunk; shrunken LFCs are for ranking and plotting),
  "depleted" with padj < alpha and log2FoldChange < -lfc_threshold.

Reads <run_dir>/results/count_matrix.csv (raw integer counts),
reference/samples.tsv + samples_meta.csv, and
<run_dir>/reference/combined_metadata.csv. Writes, under
<run_dir>/results/differential/:
  differential_results.csv   all peptides x all contrasts (long format)
  enrichment_summary.csv     per contrast: n tested / enriched / depleted, top hit
  organism_summary.csv       per contrast x organism: n enriched peptides
  <ip>/<treatment>_vs_<reference>.csv   one contrast, sorted by padj
  <ip>/size_factors.csv      per contrast x sample
  <ip>/*.png                 volcano, MA, p-value histograms, dispersions,
                             enriched-peptide heatmap, enrichment counts,
                             organism dot plot
  summary.md                 the above, as a readable report

Usage: differential_enrichment.py [samples_tsv] [samples_meta_csv]
  (default to reference/samples.tsv and reference/samples_meta.csv)
Threads: SLURM_CPUS_PER_TASK if set, else 1.
"""
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from count_filters import GROUP_FILTER_DESCRIPTION, passes_filter
from sample_groups import legend_outside, load_samples
from top_peptides import draw_peptide_heatmap, peptide_label, short_name

PROJECT_DIR = Path(__file__).resolve().parent.parent

# run_dir (where results/ and generated reference data live) is configured
# in reference/paths.json, a hand-maintained config file kept in the repo
# checkout -- separate from run_dir itself.
_paths = json.loads((PROJECT_DIR / "reference" / "paths.json").read_text())
RUN_DIR = Path(_paths["run_dir"])
RESULTS_DIR = RUN_DIR / "results"
OUT_DIR = RESULTS_DIR / "differential"
META_CSV = RUN_DIR / "reference" / "combined_metadata.csv"

_de = json.loads((PROJECT_DIR / "reference" / "de_params.json").read_text())
REFERENCE = _de["reference_treatment"]
ALPHA = _de["alpha"]
LFC_THRESHOLD = _de["lfc_threshold"]
SIZE_FACTORS_FIT_TYPE = _de["size_factors_fit_type"]
SHRINK_LFC = _de["shrink_lfc"]
N_LABEL_TOP = _de["n_label_top"]
HEATMAP_MAX_PEPTIDES = _de["heatmap_max_peptides"]
ORGANISM_PLOT_MAX = _de["organism_plot_max_organisms"]

N_CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
REF_CODE = "a_ref"  # model level names: formula-safe, and the reference sorts first (treatment coding baseline)
LFC_BOUND = 30 / np.log(2)  # PyDESeq2 clamps |log fold change| at 30 natural-log units
NOT_SIG_COLOR = "lightgray"
DEPLETED_COLOR = "dimgray"
PANEL_NCOLS = 3
PLOT_LFC = "log2FoldChange_shrunk" if SHRINK_LFC else "log2FoldChange"
CALIBRATION_WARN_FRAC = 0.10  # fraction of raw p < alpha above which a contrast's p-values look anti-conservative


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9.]+", "_", text).strip("_")


def fit_contrast(counts_all, ip_samples, treatment):
    """Fit one DESeq2 model for one contrast (treatment vs REFERENCE),
    on those samples only. Returns (dds, results DataFrame, tested index)."""
    cols = ip_samples.loc[ip_samples["treatment"].isin([treatment, REFERENCE]), "tc_id"]
    counts = counts_all[cols]
    tested = counts.index[passes_filter(counts)]
    metadata = pd.DataFrame(
        {"condition": np.where(ip_samples.set_index("tc_id").loc[cols, "treatment"] == REFERENCE,
                               REF_CODE, "t01")},
        index=cols.values)
    inference = DefaultInference(n_cpus=N_CPUS)
    dds = DeseqDataSet(counts=counts.loc[tested].T.astype(int), metadata=metadata, design="~condition",
                       size_factors_fit_type=SIZE_FACTORS_FIT_TYPE, refit_cooks=False,
                       inference=inference, quiet=True)
    dds.deseq2()

    ds = DeseqStats(dds, contrast=["condition", "t01", REF_CODE], alpha=ALPHA, cooks_filter=False,
                    independent_filter=True, inference=inference, quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    # at PyDESeq2's optimiser bound: one group is all-zero, so the MLE is
    # really +-infinity -- flagged rather than reported as a real effect size
    res["lfc_unbounded"] = res["log2FoldChange"].abs() >= LFC_BOUND - 1e-6
    if SHRINK_LFC:
        ds.lfc_shrink(coeff="condition[T.t01]")
        res["log2FoldChange_shrunk"] = ds.results_df["log2FoldChange"]
        res["lfcSE_shrunk"] = ds.results_df["lfcSE"]
    n_sig = int((res["padj"] < ALPHA).sum())
    print(f"  {treatment} vs {REFERENCE}: {len(tested)} tested ({len(cols)} samples), {n_sig} padj < {ALPHA}, "
          f"{int(res['lfc_unbounded'].sum())} unbounded LFC", file=sys.stderr)
    return dds, res, tested


def contrast_table(res, treatment, ip, ip_samples, all_refs, tested_refs, cpm, meta):
    """Full per-peptide table for one contrast: every own-library peptide,
    untested (pre-filtered) ones with NA statistics."""
    trt_ids = ip_samples.loc[ip_samples["treatment"] == treatment, "tc_id"]
    ref_ids = ip_samples.loc[ip_samples["treatment"] == REFERENCE, "tc_id"]
    out = pd.DataFrame(index=pd.Index(all_refs, name="ref_id"))
    out["ip"] = ip
    out["treatment"] = treatment
    out["reference"] = REFERENCE
    out["n_treatment"] = len(trt_ids)
    out["n_reference"] = len(ref_ids)
    out = out.join(res.reindex(all_refs))
    out["filtered"] = ~out.index.isin(tested_refs)
    out["lfc_unbounded"] = out["lfc_unbounded"].fillna(False).astype(bool)
    out["mean_cpm_treatment"] = cpm.loc[all_refs, trt_ids].mean(axis=1)
    out["mean_cpm_reference"] = cpm.loc[all_refs, ref_ids].mean(axis=1)
    sig = out["padj"] < ALPHA
    out["enriched"] = sig & (out["log2FoldChange"] > LFC_THRESHOLD)
    out["depleted"] = sig & (out["log2FoldChange"] < -LFC_THRESHOLD)
    ann = meta.loc[all_refs]
    out["peptide_label"] = [peptide_label(a) for _, a in ann.iterrows()]
    for col in ["organism", "protein_name", "start", "end", "peptide_aa"]:
        out[col] = ann[col].values
    return out.reset_index().sort_values(["padj", "pvalue"], na_position="last")


def panel_grid(n):
    ncols = min(PANEL_NCOLS, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.3 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    for ax in axes_flat[n:]:
        ax.axis("off")
    return fig, axes_flat


def status_colors(tab, color):
    """Nx4 RGBA: treatment colour if enriched, grey if depleted, pale otherwise.
    (Built as an array, not np.where -- an RGBA tuple would broadcast wrongly.)"""
    out = np.tile(mcolors.to_rgba(NOT_SIG_COLOR), (len(tab), 1))
    out[tab["depleted"].to_numpy(dtype=bool)] = mcolors.to_rgba(DEPLETED_COLOR)
    out[tab["enriched"].to_numpy(dtype=bool)] = mcolors.to_rgba(color)
    return out


def status_legend(fig):
    fig.legend(handles=[Patch(color="gray", label=f"enriched (padj < {ALPHA}, log2FC > {LFC_THRESHOLD}); "
                                                  "treatment colour"),
                        Patch(color=DEPLETED_COLOR, label=f"depleted (padj < {ALPHA}, log2FC < -{LFC_THRESHOLD})"),
                        Patch(color=NOT_SIG_COLOR, label="not significant")],
               loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=9)


def plot_volcano(tables, ip, colors, out_dir):
    lfc = PLOT_LFC
    fig, axes = panel_grid(len(tables))
    for ax, (treatment, tab) in zip(axes, tables.items()):
        tab = tab.loc[tab["padj"].notna()]
        y = -np.log10(tab["padj"].clip(lower=1e-300))
        ax.scatter(tab[lfc], y, c=status_colors(tab, colors[treatment]), s=4, linewidths=0, rasterized=True)
        ax.axhline(-np.log10(ALPHA), color="black", linestyle=":", linewidth=0.8)
        for x in (-LFC_THRESHOLD, LFC_THRESHOLD):
            ax.axvline(x, color="black", linestyle=":", linewidth=0.8)
        # one label per organism (BH ties otherwise stack many labels on one row),
        # staggered downwards from the top-right so they do not overlap each other
        top = (tab.loc[tab["enriched"]].assign(_org=lambda d: d["organism"].fillna("unknown"))
               .drop_duplicates("_org").head(N_LABEL_TOP))
        for k, (_, r) in enumerate(top.iterrows()):
            ax.annotate(f"{r['ref_id'].split('|')[-1]} {short_name(r['organism'], 22)}",
                        (r[lfc], -np.log10(max(r["padj"], 1e-300))),
                        xytext=(6, -11 * k - 2), textcoords="offset points", fontsize=6,
                        ha="left", va="top", annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", linewidth=0.4, color="gray",
                                        shrinkA=0, shrinkB=1))
        ax.margins(x=0.30)
        n_up, n_down = int(tab["enriched"].sum()), int(tab["depleted"].sum())
        flags = []
        if tab["n_treatment"].iloc[0] < 2:
            flags.append("n=1")
        frac = (tab["pvalue"] < ALPHA).sum() / max(tab["pvalue"].notna().sum(), 1)
        if frac > CALIBRATION_WARN_FRAC:
            flags.append(f"{frac:.0%} of raw p < {ALPHA}: anti-conservative")
        suffix = ("\n[" + "; ".join(flags) + "]") if flags else ""
        ax.set_title(f"{treatment} vs {REFERENCE}\n{n_up} enriched, {n_down} depleted{suffix}",
                     fontsize=10, color="firebrick" if flags else "black")
        ax.set_xlabel("log2 fold change" + (" (shrunken)" if SHRINK_LFC else ""))
        ax.set_ylabel("-log10(adjusted p)")
    fig.suptitle(f"Differential enrichment vs {REFERENCE} ({ip})", fontsize=13)
    status_legend(fig)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(out_dir / "volcano.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_ma(tables, ip, colors, out_dir):
    lfc = PLOT_LFC
    fig, axes = panel_grid(len(tables))
    for ax, (treatment, tab) in zip(axes, tables.items()):
        tab = tab.loc[tab["baseMean"].notna()]
        ax.scatter(tab["baseMean"], tab[lfc], c=status_colors(tab, colors[treatment]), s=4, linewidths=0,
                   rasterized=True)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_title(f"{treatment} vs {REFERENCE}", fontsize=10)
        ax.set_xlabel("mean of normalised counts")
        ax.set_ylabel("log2 fold change" + (" (shrunken)" if SHRINK_LFC else ""))
    fig.suptitle(f"MA plots vs {REFERENCE} ({ip})", fontsize=13)
    status_legend(fig)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(out_dir / "ma_plots.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pvalue_histograms(tables, ip, colors, out_dir):
    fig, axes = panel_grid(len(tables))
    for ax, (treatment, tab) in zip(axes, tables.items()):
        ax.hist(tab["pvalue"].dropna(), bins=50, color=colors[treatment], edgecolor="black", linewidth=0.3)
        ax.set_title(f"{treatment} vs {REFERENCE}", fontsize=10)
        ax.set_xlabel("raw p-value")
        ax.set_ylabel("peptides")
    fig.suptitle(f"p-value distributions ({ip}): expect flat + a spike near 0", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir / "pvalue_histograms.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_dispersions(dds_by_treatment, ip, out_dir):
    """One dispersion panel per contrast (each is its own fit)."""
    fig, axes = panel_grid(len(dds_by_treatment))
    for ax, (treatment, dds) in zip(axes, dds_by_treatment.items()):
        var = dds.var
        means = var["_normed_means"]
        ax.scatter(means, var["genewise_dispersions"], s=2, color="black", alpha=0.3, linewidths=0,
                   rasterized=True, label="peptide-wise estimate")
        ax.scatter(means, var["dispersions"], s=2, color="tab:blue", alpha=0.3, linewidths=0, rasterized=True,
                   label="final (MAP)")
        order = np.argsort(means.values)
        ax.plot(means.values[order], var["fitted_dispersions"].values[order], color="tab:red", linewidth=1.5,
                label="fitted trend")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{treatment} vs {REFERENCE}", fontsize=10)
        ax.set_xlabel("mean of normalised counts")
        ax.set_ylabel("dispersion")
        ax.legend(markerscale=4, fontsize=7)
    fig.suptitle(f"Dispersion estimates ({ip})", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir / "dispersions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_enriched_heatmap(tables, ip, ip_samples, cpm, meta, colors, out_dir):
    """Union of enriched peptides across this IP's contrasts (capped at
    HEATMAP_MAX_PEPTIDES by best padj), log10(CPM + 1) across all of the
    IP's samples."""
    best = pd.concat([t.loc[t["enriched"], ["ref_id", "padj"]] for t in tables.values()])
    if best.empty:
        print(f"  no enriched peptides in {ip}: skipping heatmap", file=sys.stderr)
        return
    best = best.groupby("ref_id")["padj"].min().sort_values()
    n_total = len(best)
    ref_ids = best.head(HEATMAP_MAX_PEPTIDES).index
    group = {"name": ip, "samples": ip_samples, "fig_dir": out_dir, "title_suffix": f" ({ip})"}
    title = (f"Enriched peptides vs {REFERENCE} ({ip}): {len(ref_ids)} of {n_total} "
             f"(union over contrasts, best padj first; rows clustered)")
    draw_peptide_heatmap(np.log10(cpm.loc[ref_ids] + 1), meta, group, colors, title,
                         out_dir / "enriched_heatmap.png")


def plot_enrichment_counts(summary, ip, colors, out_dir):
    sub = summary.loc[summary["ip"] == ip]
    fig, ax = plt.subplots(figsize=(0.8 * len(sub) + 3, 4.5))
    x = np.arange(len(sub))
    bar_colors = [colors[t] for t in sub["treatment"]]
    ax.bar(x, sub["n_enriched"], color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.bar(x, -sub["n_depleted"], color=bar_colors, edgecolor="black", linewidth=0.5, alpha=0.4, hatch="//")
    for xi, (up, down) in enumerate(zip(sub["n_enriched"], sub["n_depleted"])):
        ax.text(xi, up, f"{up}", ha="center", va="bottom", fontsize=8)
        ax.text(xi, -down, f"{down}", ha="center", va="top", fontsize=8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sub["treatment"], rotation=45, ha="right")
    ax.set_ylabel("peptides  (enriched up, depleted down)")
    ax.set_title(f"Significant peptides vs {REFERENCE} ({ip})\n(padj < {ALPHA}, |log2FC| > {LFC_THRESHOLD})")
    pad = max(sub["n_enriched"].max(), sub["n_depleted"].max(), 1) * 0.15
    ax.set_ylim(-sub["n_depleted"].max() - pad, sub["n_enriched"].max() + pad)
    fig.savefig(out_dir / "enrichment_counts.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_organisms(org_summary, ip, colors, out_dir):
    sub = org_summary.loc[org_summary["ip"] == ip]
    if sub.empty:
        return
    top_orgs = sub.groupby("organism")["n_enriched"].sum().sort_values(ascending=False).head(ORGANISM_PLOT_MAX)
    treatments = list(dict.fromkeys(sub["treatment"]))
    sub = sub.loc[sub["organism"].isin(top_orgs.index)]
    y_of = {o: i for i, o in enumerate(top_orgs.index[::-1])}
    x_of = {t: i for i, t in enumerate(treatments)}
    fig, ax = plt.subplots(figsize=(0.8 * len(treatments) + 6, 0.3 * len(top_orgs) + 2))
    scale = 400 / max(sub["n_enriched"].max(), 1)
    ax.scatter([x_of[t] for t in sub["treatment"]], [y_of[o] for o in sub["organism"]],
               s=sub["n_enriched"] * scale, c=[colors[t] for t in sub["treatment"]], edgecolor="black",
               linewidth=0.5)
    ax.set_xticks(range(len(treatments)))
    ax.set_xticklabels(treatments, rotation=45, ha="right")
    ax.set_yticks(range(len(y_of)))
    # strains are often distinguished only by the parenthesised part, so keep it here
    ax.set_yticklabels([short_name(o, 52, strip_parens=False) for o in top_orgs.index[::-1]], fontsize=8)
    ax.set_xlim(-0.6, len(treatments) - 0.4)
    ax.grid(color="lightgray", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_title(f"Enriched peptides by organism vs {REFERENCE} ({ip})\n"
                 f"(top {len(top_orgs)} organisms by total enriched peptides)")
    sizes = sorted({int(v) for v in np.quantile(sub["n_enriched"], [0.5, 0.9, 1.0])})
    legend_outside(ax, [Line2D([], [], marker="o", linestyle="", color="gray", markeredgecolor="black",
                               markersize=np.sqrt(s * scale), label=f"{s} peptides") for s in sizes],
                   fontsize=8, labelspacing=1.5)
    fig.savefig(out_dir / "organism_enrichment.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_summary_md(summary, org_summary, ips):
    lines = [f"# Differential enrichment vs {REFERENCE}", "",
             f"One DESeq2-style negative-binomial fit per contrast (treatment vs {REFERENCE}, within each IP), "
             f"two-sided Wald test, Benjamini-Hochberg adjusted p-values.",
             f"Peptides are filtered within each contrast ({GROUP_FILTER_DESCRIPTION}); the rest are reported "
             f"with `filtered=True` and NA statistics.",
             f"**Enriched** = padj < {ALPHA} and log2FC > {LFC_THRESHOLD}; **depleted** = padj < {ALPHA} and "
             f"log2FC < -{LFC_THRESHOLD}.", "",
             "With 2 replicates per treatment, power is limited and DESeq2's Cook's-distance outlier filtering "
             "cannot run (it needs >= 3). `lfc_unbounded` marks peptides absent from one group entirely, whose "
             "fold change has no finite estimate.", ""]
    bad = summary.loc[summary["frac_p_below_alpha"] > CALIBRATION_WARN_FRAC]
    if not bad.empty:
        lines += [f"> **Calibration warning.** In {len(bad)} of {len(summary)} contrasts, "
                  f"{bad['frac_p_below_alpha'].min():.0%}-{bad['frac_p_below_alpha'].max():.0%} of raw p-values "
                  f"fall below {ALPHA} (expected ~{ALPHA:.0%} under the null), so the model understates "
                  f"within-group variability and these adjusted p-values are anti-conservative -- the hit counts "
                  f"below are upper bounds. Affected: "
                  + ", ".join(f"{r.ip} {r.treatment}" for r in bad.itertuples()) + ". "
                  f"Check the p-value histograms; note the poor replicate correlation in the QC figures.", ""]
    lines += [
             "Full table (all peptides x all contrasts): [differential_results.csv](differential_results.csv) | "
             "per-organism counts: [organism_summary.csv](organism_summary.csv)", "",
             "## Contrast summary", "",
             summary.to_markdown(index=False, floatfmt=".3g"), ""]
    for ip in ips:
        lines += [f"## {ip}", ""]
        sub = org_summary.loc[org_summary["ip"] == ip]
        if not sub.empty:
            top = (sub.groupby("organism")["n_enriched"].sum().sort_values(ascending=False).head(10)
                   .rename("enriched peptides (all contrasts)").reset_index())
            lines += ["Most-enriched organisms:", "", top.to_markdown(index=False), ""]
        for name in ["volcano", "ma_plots", "enrichment_counts", "organism_enrichment", "enriched_heatmap",
                     "pvalue_histograms", "dispersions"]:
            lines += [f"![{name}]({ip}/{name}.png)", ""]
    (OUT_DIR / "summary.md").write_text("\n".join(lines))
    print(f"Wrote {OUT_DIR / 'summary.md'}", file=sys.stderr)


def main():
    samples, colors = load_samples()
    meta = pd.read_csv(META_CSV, dtype=str).set_index("ref_id")
    mat = pd.read_csv(RESULTS_DIR / "count_matrix.csv", index_col=0)
    libraries = samples["library"].unique()
    if len(libraries) != 1:
        sys.exit(f"ERROR: differential_enrichment.py expects samples from one library, got {', '.join(libraries)}")
    own_refs = meta.index[meta["library"] == libraries[0]].intersection(mat.index)
    counts_all = mat.loc[own_refs, samples["tc_id"]]
    cpm = counts_all / counts_all.sum(axis=0) * 1e6
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_tables, summary_rows = [], []
    ip_tables = {}
    for ip, ip_samples in samples.groupby("ip", sort=True):
        ip_samples = ip_samples.reset_index(drop=True)
        if REFERENCE not in set(ip_samples["treatment"]):
            print(f"Skipping {ip}: no '{REFERENCE}' samples", file=sys.stderr)
            continue
        out_dir = OUT_DIR / ip
        out_dir.mkdir(parents=True, exist_ok=True)
        counts_ip = counts_all[ip_samples["tc_id"]]
        treatments = [t for t in dict.fromkeys(ip_samples["treatment"]) if t != REFERENCE]
        print(f"{ip}: {len(treatments)} contrasts vs {REFERENCE}, filtered per contrast "
              f"({GROUP_FILTER_DESCRIPTION}), size factors: {SIZE_FACTORS_FIT_TYPE}, {N_CPUS} CPUs",
              file=sys.stderr)

        tables, dds_by_treatment, size_factor_rows = {}, {}, []
        for treatment in treatments:
            dds, res, tested = fit_contrast(counts_ip, ip_samples, treatment)
            dds_by_treatment[treatment] = dds
            size_factor_rows.append(pd.DataFrame({"treatment": treatment, "tc_id": dds.obs_names,
                                                  "size_factor": dds.obs["size_factors"].values}))
            tab = contrast_table(res, treatment, ip, ip_samples, own_refs, tested, cpm, meta)
            tab.to_csv(out_dir / f"{safe_name(treatment)}_vs_{safe_name(REFERENCE)}.csv", index=False)
            tables[treatment] = tab
            enriched = tab.loc[tab["enriched"]]
            summary_rows.append({
                "ip": ip, "treatment": treatment, "reference": REFERENCE,
                "n_treatment": int(tab["n_treatment"].iloc[0]), "n_reference": int(tab["n_reference"].iloc[0]),
                "n_tested": int((~tab["filtered"]).sum()),
                "n_lfc_unbounded": int(tab["lfc_unbounded"].sum()),
                "frac_p_below_alpha": float((tab["pvalue"] < ALPHA).sum() / max(tab["pvalue"].notna().sum(), 1)),
                "n_with_padj": int(tab["padj"].notna().sum()),
                "n_enriched": int(tab["enriched"].sum()), "n_depleted": int(tab["depleted"].sum()),
                "top_enriched_peptide": enriched["peptide_label"].iloc[0] if len(enriched) else "",
                "top_enriched_log2FC": enriched["log2FoldChange"].iloc[0] if len(enriched) else np.nan,
                "top_enriched_padj": enriched["padj"].iloc[0] if len(enriched) else np.nan,
            })
        pd.concat(size_factor_rows).to_csv(out_dir / "size_factors.csv", index=False)
        all_tables.extend(tables.values())
        ip_tables[ip] = (tables, ip_samples, dds_by_treatment)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "enrichment_summary.csv", index=False)
    combined = pd.concat(all_tables, ignore_index=True)
    combined.to_csv(OUT_DIR / "differential_results.csv", index=False)
    org_summary = (combined.loc[combined["enriched"]]
                   .assign(organism=lambda d: d["organism"].fillna("unknown"))
                   .groupby(["ip", "treatment", "organism"]).size().rename("n_enriched").reset_index()
                   .sort_values(["ip", "treatment", "n_enriched"], ascending=[True, True, False]))
    org_summary.to_csv(OUT_DIR / "organism_summary.csv", index=False)
    print(summary[["ip", "treatment", "n_tested", "n_enriched", "n_depleted"]].to_string(index=False),
          file=sys.stderr)

    for ip, (tables, ip_samples, dds_by_treatment) in ip_tables.items():
        out_dir = OUT_DIR / ip
        plot_volcano(tables, ip, colors, out_dir)
        plot_ma(tables, ip, colors, out_dir)
        plot_pvalue_histograms(tables, ip, colors, out_dir)
        plot_dispersions(dds_by_treatment, ip, out_dir)
        plot_enriched_heatmap(tables, ip, ip_samples, cpm[ip_samples["tc_id"]], meta, colors, out_dir)
        plot_enrichment_counts(summary, ip, colors, out_dir)
        plot_organisms(org_summary, ip, colors, out_dir)
        print(f"Wrote figures to {out_dir}", file=sys.stderr)
    write_summary_md(summary, org_summary, list(ip_tables))
    print(f"Wrote {OUT_DIR / 'differential_results.csv'} ({len(combined)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
