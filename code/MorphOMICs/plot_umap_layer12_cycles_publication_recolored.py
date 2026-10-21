#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Rebuild the layer-12 UMAP panels (per phase).

The shared display limits are estimated robustly from the archived reduced
MorphOMICs coordinates. Control is grey and CCI is raspberry.
"""
from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(os.path.abspath(__file__)).parent
BUNDLE_ROOT = SCRIPT_DIR.parent.parent
sys.path.append(str(BUNDLE_ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402
ns.set_style()
PUB = {"title": ns.FS_TITLE, "axis_title": ns.FS_AXIS,
       "axis_text": ns.FS_TICK + 2, "legend": ns.FS_LEGEND}

REDUCED = (BUNDLE_ROOT / "data" / "MorphOMICs" / "saved_instances_cam" /
           "Morphomics.PID_CAM_CCI_260224.cam_reduced.csv")
OUT_DIR = BUNDLE_ROOT / "generated_outputs" / "MorphOMICs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Panels ordered Male -> Proestrus -> Estrus -> Diestrus.
PHASE_ORDER = ["Male", "Proestrus", "Estrus", "Diestrus"]
CYCLE_EN = {"Proestrus": "Proestrus", "Estrus": "Estrus", "Diestrus": "Diestrus", "Male": "Male"}
COND_COLORS = {"Control": "#999999", "CCI": "#CC79A7"}  # unified condition palette
COND_ORDER = ["Control", "CCI"]

df = pd.read_csv(REDUCED, index_col=0)
df["Phase"] = df["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
df["Condition"] = df["Condition"].astype(str)
d12 = df[df["Layer"].astype(str) == "12"].copy()
ux, uy = "umap___1", "umap___2"

cycles = [c for c in PHASE_ORDER if c in set(d12["Phase"])]

# Robust shared display limits for the publication panels.
median_xy = d12[[ux, uy]].median()
mad_xy = (d12[[ux, uy]] - median_xy).abs().median().replace(0, 1)
robust_distance = np.sqrt((((d12[[ux, uy]] - median_xy) / mad_xy) ** 2).sum(axis=1))
ROBUST_DISTANCE_CUTOFF = 6.0
d12["_umap_in_view"] = robust_distance <= ROBUST_DISTANCE_CUTOFF
visible = d12[d12["_umap_in_view"]].copy()

xlim = (visible[ux].min() - 0.5, visible[ux].max() + 0.5)
ylim = (visible[uy].min() - 0.5, visible[uy].max() + 0.5)

crop_report = (
    d12.groupby(["Phase", "Condition"])["_umap_in_view"]
    .agg(total_points="size", shown_points="sum")
    .reset_index()
)
crop_report["cropped_points"] = crop_report["total_points"] - crop_report["shown_points"]
crop_report["robust_distance_cutoff"] = ROBUST_DISTANCE_CUTOFF
crop_report.to_csv(OUT_DIR / "cam_umap_layer12_visual_crop_report.csv", index=False)

fig, axes = plt.subplots(1, len(cycles),
                         figsize=(ns.mm2inch(40 * len(cycles)), ns.mm2inch(46)),
                         sharex=True, sharey=True)
if len(cycles) == 1:
    axes = [axes]
fig.patch.set_facecolor("white")

for k, (ax, cycle) in enumerate(zip(axes, cycles)):
    sub_all = d12[d12["Phase"] == cycle]
    sub = sub_all[sub_all["_umap_in_view"]]
    for cond in COND_ORDER:
        s = sub[sub["Condition"] == cond]
        for point in s.itertuples(index=False):
            ax.plot(
                [getattr(point, ux)],
                [getattr(point, uy)],
                marker="o",
                linestyle="none",
                markersize=np.sqrt(22),
                markerfacecolor=COND_COLORS[cond],
                markeredgecolor="none",
                alpha=0.80,
                zorder=3 if cond == "CCI" else 2,
            )
    # Sample counts and points outside the display limits are exported separately.
    ax.set_title(CYCLE_EN.get(cycle, cycle), fontsize=PUB["title"], fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=PUB["axis_title"], fontweight="bold")
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.tick_params(labelsize=ns.FS_TICK)
    if k == 0:
        ax.set_ylabel("UMAP 2", fontsize=PUB["axis_title"], fontweight="bold")

handles = [Line2D([0], [0], marker="o", linestyle="", markersize=11,
                  markerfacecolor=COND_COLORS[c], markeredgecolor="none", label=c)
           for c in COND_ORDER]
fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
           fontsize=PUB["legend"] + 1, bbox_to_anchor=(0.5, 1.06),
           handletextpad=0.3, columnspacing=1.2)
fig.subplots_adjust(left=0.055, right=0.995, top=0.82, bottom=0.16, wspace=0.07)

png = OUT_DIR / "cam_umap_layer12_cycles_publication.png"
svg = OUT_DIR / "cam_umap_layer12_cycles_publication.svg"
for cycle, target_ax in zip(cycles, axes):
    for candidate_ax in axes:
        candidate_ax.set_visible(candidate_ax is target_ax)
    for legend in fig.legends:
        legend.set_visible(False)
    part_base = OUT_DIR / f"cam_umap_layer12_{CYCLE_EN.get(cycle, cycle)}_publication"
    svg_part = part_base.with_suffix(".svg")
    pdf_part = part_base.with_suffix(".pdf")
    fig.savefig(
        svg_part,
        format="svg",
        bbox_inches="tight",
        pad_inches=0.12,
        facecolor="white",
        transparent=False,
    )
    fig.savefig(
        pdf_part,
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.12,
        facecolor="white",
        transparent=False,
    )
    print("[OK]", svg_part)
    print("[OK]", pdf_part)
for candidate_ax in axes:
    candidate_ax.set_visible(True)
for legend in fig.legends:
    legend.set_visible(True)
plt.savefig(png, dpi=300, bbox_inches="tight")
plt.savefig(svg, bbox_inches="tight")
plt.close("all")
print("[OK]", png)
print("[OK]", svg)
