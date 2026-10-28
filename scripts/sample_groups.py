"""Shared sample-sheet handling for the figure-generating scripts: merges
reference/samples.tsv with the per-sample metadata sheet
(reference/samples_meta.csv: sample, ip, treatment), fixes one
treatment -> colour mapping for every figure, and defines the figure
groups each script draws (all samples together, then one group per IP
condition, written to results/figures/<ip>/).

Imported by summarize.py, completeness_figures.py,
evenness_summary_figure.py and overview_figure.py (all run as
`python scripts/<name>.py`, so this sibling module is importable).
"""
# Note: top_peptides.py and differential_enrichment.py also import this module.
import json  # parses reference/figure_params.json
import sys  # sys.argv (optional CLI paths) and sys.exit (error out)
from pathlib import Path  # object-oriented file paths

import pandas as pd  # sample sheet tables
import matplotlib.pyplot as plt  # only used for the tab20 colour map
from matplotlib.lines import Line2D  # invisible "proxy" artists used as legend entries

# Repo root (this file is scripts/sample_groups.py, so up two levels).
PROJECT_DIR = Path(__file__).resolve().parent.parent
# Default sample sheet: tab-separated, columns tc_id, library, sample_prefix.
SAMPLES_TSV = PROJECT_DIR / "reference" / "samples.tsv"
# Default per-sample metadata: comma-separated, columns sample, ip, treatment.
SAMPLES_META_CSV = PROJECT_DIR / "reference" / "samples_meta.csv"

# Figure settings shared by all figure scripts.
_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
# Optional {"<treatment>": "<matplotlib colour>"} overrides; .get(..., {})
# means "no overrides" when the key is absent from the JSON.
TREATMENT_COLOR_OVERRIDES = _fig_params.get("treatment_colors", {})
# (x, y) in axes-fraction coordinates where legends are pinned, e.g. (1.02, 0.5)
# = just right of the plot, vertically centred. JSON lists become tuples here.
LEGEND_BBOX_ANCHOR = tuple(_fig_params["legend_bbox_anchor"])

# Line styles / marker shapes for the per-group "style" column: ip in the
# all-samples figures, replicate number in the per-IP figures.
LINESTYLE_CYCLE = ["-", "--", ":", "-."]  # solid, dashed, dotted, dash-dot
# Marker shapes: circle, square, triangle, diamond, down-triangle, plus, x,
# star, hexagon, left/right triangles. Reused (cycled) if there are more values.
MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]


# load_samples()
#   In:  optional CLI args: sys.argv[1] = sample sheet path, sys.argv[2] =
#        metadata path (defaults: SAMPLES_TSV / SAMPLES_META_CSV).
#   Out: (samples DataFrame: tc_id, library, sample_prefix, ip, treatment,
#         replicate - sorted in plotting order; dict treatment -> colour).
#   Called by: summarize.py, completeness_figures.py, evenness_summary_figure.py,
#              overview_figure.py, top_peptides.py, differential_enrichment.py.
def load_samples():
    """samples.tsv merged with samples_meta.csv (paths overridable as the
    first/second CLI arguments), in plotting order: by ip, then treatment
    (in order of first appearance in samples_meta.csv), then the original
    samples.tsv order. Adds a `replicate` column (1, 2, ... within each
    ip + treatment). Returns (samples, treatment_colors)."""
    # sys.argv[0] is the script name, so argv[1]/argv[2] are the first/second
    # user arguments. Fall back to the repo defaults when not given.
    samples_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    meta_path = Path(sys.argv[2]) if len(sys.argv) > 2 else SAMPLES_META_CSV
    samples = pd.read_csv(samples_path, sep="\t")  # sep="\t": tab-separated file
    sample_meta = pd.read_csv(meta_path, dtype=str)  # dtype=str: keep every value as text

    # Left join: keep every samples.tsv row and attach its ip/treatment by
    # matching tc_id to the metadata's `sample` column (samples with no match
    # get NaN). Metadata rows for other samples are dropped. The duplicate key
    # column `sample` is then removed.
    samples = samples.merge(sample_meta, left_on="tc_id", right_on="sample", how="left").drop(columns=["sample"])
    # tc_ids of samples lacking an ip or a treatment (| = element-wise OR).
    # .loc[row_mask, "tc_id"] selects those rows' tc_id column.
    missing = samples.loc[samples["ip"].isna() | samples["treatment"].isna(), "tc_id"].tolist()
    if missing:
        # sys.exit with a string prints it to stderr and exits with status 1.
        sys.exit(f"ERROR: no ip/treatment in {meta_path} for: {', '.join(missing)}")

    # Unique treatments in order of first appearance in samples_meta.csv.
    # dict.fromkeys() de-duplicates while preserving order (unlike set()).
    treatment_order = list(dict.fromkeys(sample_meta["treatment"]))
    # Temporary column recording the original samples.tsv row order (0, 1, 2...),
    # used as the final tie-breaker when sorting.
    samples["_order"] = range(len(samples))
    # Make treatment an *ordered categorical* so sorting follows treatment_order
    # (samples_meta.csv order) rather than alphabetical order.
    samples["treatment"] = pd.Categorical(samples["treatment"], categories=treatment_order, ordered=True)
    # Sort by ip (alphabetical), then treatment (category order), then original
    # order; drop the helper column; renumber the index 0..n-1.
    samples = samples.sort_values(["ip", "treatment", "_order"]).drop(columns=["_order"]).reset_index(drop=True)
    # Back to plain strings so downstream code can treat treatment as text.
    samples["treatment"] = samples["treatment"].astype(str)
    # Replicate number within each ip + treatment combination: cumcount()
    # numbers rows 0, 1, 2... within each group (in the sorted order), +1 makes
    # it 1-based, and astype(str) makes it a label ("1", "2", ...) for legends.
    samples["replicate"] = (samples.groupby(["ip", "treatment"]).cumcount() + 1).astype(str)

    # Colours are built from *all* treatments in samples_meta.csv (not just
    # those in this sample sheet), so they stay the same across subsets.
    return samples, treatment_colors(treatment_order)


# treatment_colors(treatment_order)
#   In:  list of treatment names, in the desired colour-assignment order.
#   Out: dict treatment -> matplotlib colour (RGBA tuple or override string).
#   Called by: load_samples() above.
def treatment_colors(treatment_order):
    """One colour per treatment, fixed across every figure (so IP4 and IP5
    figures are directly comparable), from tab20 in samples_meta.csv
    order; figure_params.json's optional `treatment_colors` overrides
    individual entries."""
    # tab20: matplotlib's qualitative palette of 20 distinct colours.
    cmap = plt.get_cmap("tab20")
    # Dict comprehension: i-th treatment -> i-th tab20 colour. cmap(int) picks
    # that palette entry; `% cmap.N` (N = 20) wraps around, so a 21st
    # treatment would reuse the first colour.
    colors = {t: cmap(i % cmap.N) for i, t in enumerate(treatment_order)}
    # Replace any entries given explicitly in figure_params.json.
    colors.update(TREATMENT_COLOR_OVERRIDES)
    return colors


# figure_groups(samples, fig_dir)
#   In:  samples DataFrame from load_samples(); fig_dir = base figures Path.
#   Out: list of group dicts: first "all", then one per ip (alphabetical).
#        Creates each group's output directory as a side effect.
#   Called by: summarize.py, completeness_figures.py, evenness_summary_figure.py,
#              overview_figure.py, top_peptides.py.
def figure_groups(samples, fig_dir):
    """The sample subsets each figure is drawn for: all samples (written to
    fig_dir, styled by ip) and each IP condition on its own (written to
    fig_dir/<ip>/, styled by replicate). Each group is a dict with keys
    name, samples, fig_dir, title_suffix, style_col."""
    # The all-samples group: line style / marker distinguishes the IPs.
    groups = [{"name": "all", "samples": samples, "fig_dir": fig_dir,
               "title_suffix": " (all samples)", "style_col": "ip"}]
    # groupby("ip") yields (ip value, sub-DataFrame) pairs, sorted by ip.
    for ip, sub in samples.groupby("ip", sort=True):
        # Per-IP group, written to e.g. figures/IP4/. Within one IP the style
        # distinguishes replicates. reset_index renumbers the subset's rows 0..n-1.
        groups.append({"name": ip, "samples": sub.reset_index(drop=True), "fig_dir": fig_dir / ip,
                       "title_suffix": f" ({ip})", "style_col": "replicate"})
    # Make sure every output folder exists (no error if it already does).
    for g in groups:
        g["fig_dir"].mkdir(parents=True, exist_ok=True)
    return groups


# style_map(group, cycle)
#   In:  a group dict from figure_groups(); cycle = LINESTYLE_CYCLE or MARKER_CYCLE.
#   Out: dict style-column value -> line style / marker.
#   Called by: summarize.py, completeness_figures.py.
def style_map(group, cycle):
    """Map each value of the group's style column to an entry of `cycle`."""
    # Distinct values of the style column (ip names or replicate labels),
    # sorted as strings so the assignment is deterministic.
    values = sorted(group["samples"][group["style_col"]].unique())
    # i-th value -> i-th style, wrapping round the cycle if it runs out.
    return {v: cycle[i % len(cycle)] for i, v in enumerate(values)}


# sample_label(srow, group)
#   In:  srow = one row (Series) of the samples DataFrame; group = group dict.
#   Out: label string, e.g. "IP4 PhIP_Seq1_5" (all-samples) or "PhIP_Seq1_5".
#   Called by: summarize.py, completeness_figures.py, evenness_summary_figure.py,
#              overview_figure.py, top_peptides.py.
def sample_label(srow, group):
    """Tick/title label for one sample: ip-prefixed in the all-samples
    figures, where IP4 and IP5 samples sit side by side."""
    # Conditional expression: prefix with the ip only in the "all" group.
    return f"{srow['ip']} {srow['tc_id']}" if group["name"] == "all" else srow["tc_id"]


# treatment_legend_handles(samples, colors, kind="patch")
#   In:  samples DataFrame (usually one group's), the treatment colour dict,
#        and kind: "patch" (squares), "line" (coloured lines) or anything
#        else (circles; callers pass "marker").
#   Out: list of Line2D legend handles, starting with a text-only header.
#   Called by: summarize.py, completeness_figures.py, evenness_summary_figure.py,
#              overview_figure.py, top_peptides.py.
def treatment_legend_handles(samples, colors, kind="patch"):
    """Legend handles (a header plus one entry per treatment present in
    `samples`) as filled squares (kind="patch"), lines or circles."""
    # An empty, invisible line whose label acts as a section header in the legend.
    handles = [Line2D([], [], linestyle="", label="— treatment —")]
    # Each treatment present in these samples, once, in their plotting order.
    for t in dict.fromkeys(samples["treatment"]):
        if kind == "line":
            # A short coloured line (for line plots such as rarefaction curves).
            handles.append(Line2D([], [], color=colors[t], linewidth=2, label=t))
        else:
            # No line, just a marker: square for "patch", circle otherwise.
            marker = "s" if kind == "patch" else "o"
            handles.append(Line2D([], [], marker=marker, linestyle="", color=colors[t],
                                  markeredgecolor="black", markersize=8, label=t))
    return handles


# style_legend_handles(group, styles, kind="line")
#   In:  group dict; styles = dict from style_map(); kind "line" or anything
#        else (markers).
#   Out: list of grey Line2D legend handles, starting with a header naming
#        the style column (e.g. "— replicate —").
#   Called by: summarize.py, completeness_figures.py.
def style_legend_handles(group, styles, kind="line"):
    """Legend handles explaining the group's style column (line styles or
    marker shapes)."""
    # Invisible header entry, e.g. "— ip —" or "— replicate —".
    handles = [Line2D([], [], linestyle="", label=f"— {group['style_col']} —")]
    for value, style in styles.items():
        if kind == "line":
            # Grey (colour-neutral) line in this value's line style.
            handles.append(Line2D([], [], color="gray", linestyle=style, linewidth=1.5, label=value))
        else:
            # Grey marker in this value's marker shape.
            handles.append(Line2D([], [], marker=style, linestyle="", color="gray",
                                  markeredgecolor="black", markersize=8, label=value))
    return handles


# legend_outside(target, handles, **kwargs)
#   In:  target = a matplotlib Axes or Figure; handles = legend handles;
#        **kwargs = any extra legend() options (e.g. fontsize=8), passed through.
#   Out: the created Legend object.
#   Called by: summarize.py, completeness_figures.py, evenness_summary_figure.py,
#              overview_figure.py, differential_enrichment.py.
def legend_outside(target, handles, **kwargs):
    """Legend to the right of an axes (or figure) rather than over the data."""
    # loc="center left" means the legend's left-middle point is placed at
    # LEGEND_BBOX_ANCHOR, i.e. just outside the right edge of the plot.
    return target.legend(handles=handles, loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR,
                         borderaxespad=0, **kwargs)
