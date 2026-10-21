#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the layer-12 ipsilateral persistence panels (persistence-difference panel)."""
from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(os.path.abspath(__file__)).parent
BUNDLE_ROOT = SCRIPT_DIR.parent.parent
sys.path.append(str(SCRIPT_DIR))
sys.path.append(str(BUNDLE_ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402
ns.set_style()
from morphomics_plot_utils import (  # noqa: E402
    compute_persistence_centroids,
    pick_representative_swc,
    read_swc_dataframe,
)

DATA = BUNDLE_ROOT / "data" / "MorphOMICs"
OUT_DIR = BUNDLE_ROOT / "generated_outputs" / "MorphOMICs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
prefix = "Morphomics.PID_CAM_CCI_260224"

# Unified condition palette (shared with every other figure).
COL_CCI = "#CC79A7"      # raspberry
COL_CONTROL = "#999999"  # grey

# Difference heatmap (CCI - Control): grey to white to raspberry.
from matplotlib.colors import LinearSegmentedColormap
DIFF_CMAP = LinearSegmentedColormap.from_list(
    "control_cci_diff", ["#5A5A5A", "#F4F4F4", COL_CCI])

TARGET_SIDE = "IPSI"
TARGET_LAYER = "12"

# Load the frozen persistence images from a non-executable NumPy archive.
vector_path = DATA / "saved_instances_cam" / "cam_persistence_vectors.npz"
with np.load(vector_path, allow_pickle=False) as vector_archive:
    group_columns = ["Condition", "Phase", "Sex", "Side", "Layer"]
    df_vectors = pd.DataFrame(
        {column: vector_archive[column] for column in group_columns}
    )
    df_vectors["pi"] = list(vector_archive["persistence_images"])
if "Phase" in df_vectors.columns:
    df_vectors["Phase"] = df_vectors["Phase"].replace({"Metestrus": "Diestrus"})

centroids = compute_persistence_centroids(
    df_vectors,
    group_columns=["Condition", "Phase", "Side", "Layer"],
    vector_column="pi",
)

# Manifest and representative SWC files for ipsilateral layer-12 groups.
manifest = pd.read_csv(DATA / "reports_cam" / "cam_manifest_selected.csv")
manifest["Phase"] = manifest["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
manifest["side"] = manifest["side"].astype(str)
manifest["layer"] = manifest["layer"].astype(str)
m12 = manifest[(manifest["side"] == TARGET_SIDE) & (manifest["layer"] == TARGET_LAYER)]


def to_bundled(path_value):
    """Resolve the archive-relative SWC path recorded in the manifest."""
    path = Path(str(path_value))
    return path if path.is_absolute() else BUNDLE_ROOT / path


rep_cache = {}
for (condition, cycle), sub in m12.groupby(["Condition", "Phase"]):
    paths = [to_bundled(p) for p in sub["target_path"].dropna().tolist()]
    paths = [p for p in paths if p.exists()]
    rep_cache[(str(condition), str(cycle))] = pick_representative_swc(paths) if paths else None


def plot_swc_trace(ax, swc_df, color, label):
    if swc_df is None or swc_df.empty:
        return
    if (swc_df["parent"] == -1).any():
        center = swc_df.loc[swc_df["parent"] == -1, ["x", "y"]].iloc[0].to_numpy(dtype=float)
    else:
        center = swc_df[["x", "y"]].to_numpy(dtype=float).mean(axis=0)
    node_xy = {int(r.n): np.array([r.x, r.y], dtype=float) - center
               for r in swc_df.itertuples(index=False)}
    for r in swc_df.itertuples(index=False):
        nid, pid = int(r.n), int(r.parent)
        if pid == -1 or pid not in node_xy or nid not in node_xy:
            continue
        x1, y1 = node_xy[nid]; x2, y2 = node_xy[pid]
        ax.plot([x1, x2], [y1, y2], color=color, lw=0.6, alpha=0.95)
    ax.plot([], [], color=color, lw=2.2, label=label)


def normalized_diff(a, b):
    return a / (a.max() + 1e-12) - b / (b.max() + 1e-12)


def build_pairs(cdict):
    by_ctx = {}
    for gk, data in cdict.items():
        combo = data["combo"]
        ctx = (combo.get("Phase"), combo.get("Side"), combo.get("Layer"))
        by_ctx.setdefault(ctx, {})[combo.get("Condition")] = gk
    out = []
    for ctx, cm in by_ctx.items():
        if "CCI" in cm and "Control" in cm:
            out.append((ctx, cm["CCI"], cm["Control"]))
    return sorted(out)


pairs = [p for p in build_pairs(centroids)
         if str(p[0][1]) == TARGET_SIDE and str(p[0][2]) == TARGET_LAYER]
print(f"Pairs IPSI/Layer12: {len(pairs)}")

# Font sizes (pt): titles/labels 7, ticks 5, legend 7.
TS_TITLE, TS_LABEL, TS_TICK, TS_LEG = ns.FS_TITLE, ns.FS_AXIS, ns.FS_TICK, ns.FS_LEGEND
for ctx, key_cci, key_ctrl in pairs:
    cycle, side, layer = ctx
    cci, ctrl = centroids[key_cci], centroids[key_ctrl]
    cci_img = gaussian_filter(cci["centroid"], sigma=1.2)
    ctrl_img = gaussian_filter(ctrl["centroid"], sigma=1.2)
    diff_img = gaussian_filter(normalized_diff(cci_img, ctrl_img), sigma=1.0)

    def _safe_read(p):
        return read_swc_dataframe(p) if p is not None else None
    rc, rt = rep_cache.get(("CCI", str(cycle))), rep_cache.get(("Control", str(cycle)))
    if rc is None or rt is None:
        print(f"  [warn] missing representative SWC for {cycle}: CCI={rc}, Control={rt}")
    swc_cci = _safe_read(rc)
    swc_ctrl = _safe_read(rt)

    # Physical size within a 183-mm two-column width. 4 panels in ~165 mm.
    fig, axes = plt.subplots(1, 4, figsize=(ns.mm2inch(165), ns.mm2inch(44)),
                             constrained_layout=True)
    im1 = axes[0].imshow(cci_img, origin="lower", cmap="magma")
    axes[0].set_title(f"CCI-IoN (n={cci['n_samples']})", fontsize=TS_TITLE, fontweight="bold")
    im2 = axes[1].imshow(ctrl_img, origin="lower", cmap="magma")
    axes[1].set_title(f"Control (n={ctrl['n_samples']})", fontsize=TS_TITLE, fontweight="bold")
    vmax = float(np.abs(diff_img).max()) or 1.0
    im3 = axes[2].imshow(diff_img, origin="lower", cmap=DIFF_CMAP, vmin=-vmax, vmax=vmax)
    axes[2].set_title("CCI-IoN − Control (normalized)", fontsize=TS_TITLE, fontweight="bold")
    for ax, im in zip(axes[:3], [im1, im2, im3]):
        ax.set_xlabel("Birth", fontsize=TS_LABEL, fontweight="bold")
        ax.set_ylabel("Death", fontsize=TS_LABEL, fontweight="bold")
        ax.tick_params(labelsize=TS_TICK)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=TS_TICK)

    axes[3].set_title("Representative SWC traces", fontsize=TS_TITLE, fontweight="bold")
    plot_swc_trace(axes[3], swc_cci, COL_CCI, "CCI-IoN")
    plot_swc_trace(axes[3], swc_ctrl, COL_CONTROL, "Control")
    axes[3].set_aspect("equal")
    axes[3].legend(loc="best", fontsize=TS_LEG)
    axes[3].axis("off")

    cycle_label = {"Male": "Male", "Proestrus": "Proestrus", "Estrus": "Estrus", "Diestrus": "Diestrus"}.get(str(cycle), str(cycle))
    fig.suptitle(f"Persistence comparison | Phase={cycle_label}, Side={side}, Layer={layer}",
                 fontsize=TS_TITLE, fontweight="bold")
    out_stem = f"cam_persistence_diff_Phase-{cycle_label}_Side-{side}_Layer-{layer}"
    out_name = f"{out_stem}.png"
    fig.savefig(OUT_DIR / out_name, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{out_stem}.svg", bbox_inches="tight")
    plt.close(fig)
    print("[OK]", out_name)
print("DONE")
