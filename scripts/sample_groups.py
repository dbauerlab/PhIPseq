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
import json
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT_DIR = Path(__file__).resolve().parent.parent
SAMPLES_TSV = PROJECT_DIR / "reference" / "samples.tsv"
SAMPLES_META_CSV = PROJECT_DIR / "reference" / "samples_meta.csv"

_fig_params = json.loads((PROJECT_DIR / "reference" / "figure_params.json").read_text())
TREATMENT_COLOR_OVERRIDES = _fig_params.get("treatment_colors", {})
LEGEND_BBOX_ANCHOR = tuple(_fig_params["legend_bbox_anchor"])

# Line styles / marker shapes for the per-group "style" column: ip in the
# all-samples figures, replicate number in the per-IP figures.
LINESTYLE_CYCLE = ["-", "--", ":", "-."]
MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">"]


def load_samples():
    """samples.tsv merged with samples_meta.csv (paths overridable as the
    first/second CLI arguments), in plotting order: by ip, then treatment
    (in order of first appearance in samples_meta.csv), then the original
    samples.tsv order. Adds a `replicate` column (1, 2, ... within each
    ip + treatment). Returns (samples, treatment_colors)."""
    samples_path = Path(sys.argv[1]) if len(sys.argv) > 1 else SAMPLES_TSV
    meta_path = Path(sys.argv[2]) if len(sys.argv) > 2 else SAMPLES_META_CSV
    samples = pd.read_csv(samples_path, sep="\t")
    sample_meta = pd.read_csv(meta_path, dtype=str)

    samples = samples.merge(sample_meta, left_on="tc_id", right_on="sample", how="left").drop(columns=["sample"])
    missing = samples.loc[samples["ip"].isna() | samples["treatment"].isna(), "tc_id"].tolist()
    if missing:
        sys.exit(f"ERROR: no ip/treatment in {meta_path} for: {', '.join(missing)}")

    treatment_order = list(dict.fromkeys(sample_meta["treatment"]))
    samples["_order"] = range(len(samples))
    samples["treatment"] = pd.Categorical(samples["treatment"], categories=treatment_order, ordered=True)
    samples = samples.sort_values(["ip", "treatment", "_order"]).drop(columns=["_order"]).reset_index(drop=True)
    samples["treatment"] = samples["treatment"].astype(str)
    samples["replicate"] = (samples.groupby(["ip", "treatment"]).cumcount() + 1).astype(str)

    return samples, treatment_colors(treatment_order)


def treatment_colors(treatment_order):
    """One colour per treatment, fixed across every figure (so IP4 and IP5
    figures are directly comparable), from tab20 in samples_meta.csv
    order; figure_params.json's optional `treatment_colors` overrides
    individual entries."""
    cmap = plt.get_cmap("tab20")
    colors = {t: cmap(i % cmap.N) for i, t in enumerate(treatment_order)}
    colors.update(TREATMENT_COLOR_OVERRIDES)
    return colors


def figure_groups(samples, fig_dir):
    """The sample subsets each figure is drawn for: all samples (written to
    fig_dir, styled by ip) and each IP condition on its own (written to
    fig_dir/<ip>/, styled by replicate). Each group is a dict with keys
    name, samples, fig_dir, title_suffix, style_col."""
    groups = [{"name": "all", "samples": samples, "fig_dir": fig_dir,
               "title_suffix": " (all samples)", "style_col": "ip"}]
    for ip, sub in samples.groupby("ip", sort=True):
        groups.append({"name": ip, "samples": sub.reset_index(drop=True), "fig_dir": fig_dir / ip,
                       "title_suffix": f" ({ip})", "style_col": "replicate"})
    for g in groups:
        g["fig_dir"].mkdir(parents=True, exist_ok=True)
    return groups


def style_map(group, cycle):
    """Map each value of the group's style column to an entry of `cycle`."""
    values = sorted(group["samples"][group["style_col"]].unique())
    return {v: cycle[i % len(cycle)] for i, v in enumerate(values)}


def sample_label(srow, group):
    """Tick/title label for one sample: ip-prefixed in the all-samples
    figures, where IP4 and IP5 samples sit side by side."""
    return f"{srow['ip']} {srow['tc_id']}" if group["name"] == "all" else srow["tc_id"]


def treatment_legend_handles(samples, colors, kind="patch"):
    """Legend handles (a header plus one entry per treatment present in
    `samples`) as filled squares (kind="patch"), lines or circles."""
    handles = [Line2D([], [], linestyle="", label="— treatment —")]
    for t in dict.fromkeys(samples["treatment"]):
        if kind == "line":
            handles.append(Line2D([], [], color=colors[t], linewidth=2, label=t))
        else:
            marker = "s" if kind == "patch" else "o"
            handles.append(Line2D([], [], marker=marker, linestyle="", color=colors[t],
                                  markeredgecolor="black", markersize=8, label=t))
    return handles


def style_legend_handles(group, styles, kind="line"):
    """Legend handles explaining the group's style column (line styles or
    marker shapes)."""
    handles = [Line2D([], [], linestyle="", label=f"— {group['style_col']} —")]
    for value, style in styles.items():
        if kind == "line":
            handles.append(Line2D([], [], color="gray", linestyle=style, linewidth=1.5, label=value))
        else:
            handles.append(Line2D([], [], marker=style, linestyle="", color="gray",
                                  markeredgecolor="black", markersize=8, label=value))
    return handles


def legend_outside(target, handles, **kwargs):
    """Legend to the right of an axes (or figure) rather than over the data."""
    return target.legend(handles=handles, loc="center left", bbox_to_anchor=LEGEND_BBOX_ANCHOR,
                         borderaxespad=0, **kwargs)
