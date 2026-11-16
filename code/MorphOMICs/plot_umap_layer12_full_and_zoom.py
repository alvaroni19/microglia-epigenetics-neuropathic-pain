#!/usr/bin/env python3
"""Render the ipsilateral layer-12 UMAP (full view and zoom).

The input contains one global UMAP embedding. Each point is a bootstrap bag of
12 cell barcodes sampled with replacement. Main axes use robust, data-derived
limits for detail; every ipsilateral layer-12 bag remains visible in the inset.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402


SOURCE = (
    ROOT
    / "data"
    / "MorphOMICs"
    / "saved_instances_cam"
    / "Morphomics.PID_CAM_CCI_260224.cam_reduced.csv"
)
OUTPUT = ROOT / "generated_outputs" / "MorphOMICs"

CYCLE_ORDER = ["Male", "Proestrus", "Estrus", "Diestrus"]
CYCLE_LABELS = {
    "Male": "Male",
    "Proestrus": "Proestrus",
    "Estrus": "Estrus",
    "Diestrus": "Diestrus",
}
CONDITION_ORDER = ["Control", "CCI"]
CONDITION_COLOURS = {"Control": "#999999", "CCI": "#CC79A7"}
ROBUST_DISTANCE_CUTOFF = 6.0


def io_path(path: Path) -> str:
    absolute = str(path.absolute())
    if os.name == "nt" and not absolute.startswith("\\\\?\\"):
        if absolute.startswith("\\\\"):
            return "\\\\?\\UNC\\" + absolute.lstrip("\\")
        return "\\\\?\\" + absolute
    return absolute


def padded_limits(values: pd.Series, fraction: float = 0.04) -> tuple[float, float]:
    low = float(values.min())
    high = float(values.max())
    padding = max((high - low) * fraction, 0.25)
    return low - padding, high + padding


def draw_points(ax: plt.Axes, frame: pd.DataFrame, x: str, y: str, size: float) -> None:
    for condition in CONDITION_ORDER:
        subset = frame[frame["Condition"] == condition]
        ax.scatter(
            subset[x],
            subset[y],
            s=size,
            color=CONDITION_COLOURS[condition],
            edgecolors="none",
            alpha=0.80,
            rasterized=False,
            zorder=3 if condition == "CCI" else 2,
        )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(io_path(SOURCE), index_col=0)
    data["Phase"] = data["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
    data["Condition"] = data["Condition"].astype(str)
    layer12_all = data[data["Layer"].astype(str) == "12"].copy()
    layer12 = layer12_all[layer12_all["Side"].astype(str).str.upper() == "IPSI"].copy()
    if len(data) != 1920 or len(layer12_all) != 960 or len(layer12) != 480:
        raise ValueError(
            "Expected 1,920 total bags, 960 layer-12 bags and 480 ipsilateral "
            f"layer-12 bags; found {len(data)}, {len(layer12_all)} and {len(layer12)}"
        )

    x, y = "umap___1", "umap___2"
    if not {x, y}.issubset(layer12.columns):
        raise KeyError("The archived global UMAP coordinates are missing.")

    median = layer12[[x, y]].median()
    mad = (layer12[[x, y]] - median).abs().median().replace(0, 1)
    distance = np.sqrt((((layer12[[x, y]] - median) / mad) ** 2).sum(axis=1))
    layer12["shown_in_zoom_view"] = distance <= ROBUST_DISTANCE_CUTOFF
    zoom = layer12[layer12["shown_in_zoom_view"]]
    zoom_xlim = padded_limits(zoom[x])
    zoom_ylim = padded_limits(zoom[y])
    full_xlim = padded_limits(layer12[x])
    full_ylim = padded_limits(layer12[y])

    cycles = [cycle for cycle in CYCLE_ORDER if cycle in set(layer12["Phase"])]
    ns.set_style()
    fig, axes = plt.subplots(
        1,
        len(cycles),
        figsize=(ns.mm2inch(42 * len(cycles)), ns.mm2inch(52)),
        sharex=True,
        sharey=True,
    )
    axes = np.atleast_1d(axes)

    for index, (axis, cycle) in enumerate(zip(axes, cycles)):
        cycle_all = layer12[layer12["Phase"] == cycle]
        cycle_zoom = cycle_all[cycle_all["shown_in_zoom_view"]]
        draw_points(axis, cycle_zoom, x, y, 12)
        axis.set_xlim(zoom_xlim)
        axis.set_ylim(zoom_ylim)
        axis.set_title(CYCLE_LABELS[cycle], fontsize=ns.FS_TITLE, fontweight="bold")
        axis.set_xlabel("UMAP 1", fontsize=ns.FS_AXIS, fontweight="bold")
        if index == 0:
            axis.set_ylabel("UMAP 2", fontsize=ns.FS_AXIS, fontweight="bold")

        inset = axis.inset_axes([0.58, 0.56, 0.38, 0.38])
        draw_points(inset, cycle_all, x, y, 4.5)
        inset.set_xlim(full_xlim)
        inset.set_ylim(full_ylim)
        inset.set_xticks([])
        inset.set_yticks([])
        inset.set_title("Full view", fontsize=max(ns.FS_SMALL - 0.5, 4.5), pad=1)
        inset.add_patch(
            Rectangle(
                (zoom_xlim[0], zoom_ylim[0]),
                zoom_xlim[1] - zoom_xlim[0],
                zoom_ylim[1] - zoom_ylim[0],
                fill=False,
                edgecolor="#333333",
                linewidth=0.55,
            )
        )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=5,
            markerfacecolor=CONDITION_COLOURS[condition],
            markeredgecolor="none",
            label=condition,
        )
        for condition in CONDITION_ORDER
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.01))
    fig.subplots_adjust(left=0.06, right=0.995, top=0.82, bottom=0.16, wspace=0.08)

    combined_base = OUTPUT / "cam_umap_layer12_cycles_full_and_zoom"
    fig.savefig(io_path(combined_base.with_suffix(".svg")), facecolor="white", transparent=False)
    fig.savefig(io_path(combined_base.with_suffix(".pdf")), facecolor="white", transparent=False)
    fig.savefig(io_path(combined_base.with_suffix(".png")), dpi=300, facecolor="white", transparent=False)

    for cycle, target_axis in zip(cycles, axes):
        for axis in axes:
            axis.set_visible(axis is target_axis)
        for legend in fig.legends:
            legend.set_visible(False)
        output_base = OUTPUT / f"cam_umap_layer12_{CYCLE_LABELS[cycle]}_publication"
        fig.savefig(io_path(output_base.with_suffix(".svg")), bbox_inches="tight", pad_inches=0.04, facecolor="white")
        fig.savefig(io_path(output_base.with_suffix(".pdf")), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    for axis in axes:
        axis.set_visible(True)
    for legend in fig.legends:
        legend.set_visible(True)
    plt.close(fig)

    report = (
        layer12.groupby(["Phase", "Condition"], as_index=False)
        .agg(total_bags=(x, "size"), bags_in_zoom=("shown_in_zoom_view", "sum"))
    )
    report["bags_in_full_view"] = report["total_bags"]
    report["bags_outside_zoom"] = report["total_bags"] - report["bags_in_zoom"]
    report["bootstrap_bags_per_stratum"] = 60
    report["cell_barcodes_per_bag"] = 12
    report["sampling"] = "with replacement"
    report["embedding"] = "single global UMAP"
    report["displayed_side"] = "IPSI"
    report["source"] = str(SOURCE.relative_to(ROOT)).replace("\\", "/")
    report.to_csv(io_path(OUTPUT / "cam_umap_layer12_full_and_zoom_report.csv"), index=False)

    print(f"[OK] {combined_base.with_suffix('.svg')}")
    print(f"[OK] all {len(layer12)} ipsilateral layer-12 bags are shown in full-view insets")


if __name__ == "__main__":
    main()
