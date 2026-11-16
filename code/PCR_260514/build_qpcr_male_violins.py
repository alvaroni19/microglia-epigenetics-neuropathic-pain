# -*- coding: utf-8 -*-
"""Compact violin plots for qPCR fold-change data (male, Control vs CCI).

One ~30x36 mm panel per gene, reconstructed directly from the final Delta-Ct
sheets in `invivo_qpcr_delta_ct.xlsx`, with significance recomputed in Python.

Output: generated_outputs/PCR_violins/male/<Gene>_male.tiff (+ preview PNG), plus
one review grid per gene group (pain markers / main methylation regulators /
other methylation regulators / estrogen-related genes).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
from journal_style import (  # noqa: E402
    set_style, MM, FS_TICK, FS_AXIS, FS_TITLE, FS_ANNOT,
    gene_italic, save_panel,
)
from qpcr_violin_style import (  # noqa: E402
    GP_CONTROL, GP_CCI, draw_violin, add_bracket,
)
from pcr_python_stats import vivo_stars  # noqa: E402
from qpcr_raw import SCRIPT_TO_CANONICAL, load_invivo_long  # noqa: E402

set_style()

OUT = ROOT / "generated_outputs" / "PCR_violins" / "male"
VIVO_LONG = load_invivo_long()

PANEL_W_MM = 26.0
PANEL_H_MM = 36.0

# Gene groups used for the review grids.
PAIN_MARKER_GENES = ["CSF1R", "TLR4", "TNFa", "IL10", "KCC2", "BDNF"]
METHYLATION_MAIN_GENES = ["Dnmt1", "Dnmt3a", "Gadd45g"]
METHYLATION_OTHER_GENES = ["Dnmt3b", "Gadd45a", "Gadd45b"]
ESTROGEN_GENES = ["ERa", "ERb", "GPER", "ARO"]
ALL_GENES = PAIN_MARKER_GENES + METHYLATION_MAIN_GENES + METHYLATION_OTHER_GENES + ESTROGEN_GENES


def load_male_data(gene):
    canonical = SCRIPT_TO_CANONICAL[gene]
    subset = VIVO_LONG[
        (VIVO_LONG["gene"] == canonical) & (VIVO_LONG["phase"] == "Male")
    ]
    control = subset.loc[
        subset["condition"] == "Control", "relative_expression"
    ].dropna()
    cci = subset.loc[
        subset["condition"] == "CCI", "relative_expression"
    ].dropna()
    return control.to_numpy(dtype=float), cci.to_numpy(dtype=float)


def violin_panel(gene):
    control, cci = load_male_data(gene)
    stars = vivo_stars(gene, "Male CCI vs control")

    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM))
    data = [control, cci]
    violin_lo, violin_hi = [], []
    for i, (vals, col) in enumerate(zip(data, (GP_CONTROL, GP_CCI)), start=1):
        v_lo, v_hi = draw_violin(ax, vals, i, col, width=0.75)
        violin_lo.append(v_lo)
        violin_hi.append(v_hi)

    for i, vals in enumerate(data, start=1):
        med = np.median(vals)
        ax.plot([i - 0.18, i + 0.18], [med, med], color="black", lw=1.0, zorder=5)

    rng = np.random.default_rng(0)
    for i, vals in enumerate(data, start=1):
        x = i + rng.uniform(-0.07, 0.07, size=len(vals))
        ax.scatter(x, vals, s=3, color="black", zorder=6, linewidths=0)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Control", "CCI"], fontsize=FS_TICK, rotation=45, ha="right")
    ax.set_ylabel("Fold change", fontsize=FS_AXIS, fontweight="bold")
    ax.set_title(gene_italic(gene), fontsize=FS_TITLE, fontweight="bold", pad=2.0)

    lo = min(violin_lo)
    hi_data = max(control.max(), cci.max())
    hi_violin = max(violin_hi)
    span = hi_data - min(control.min(), cci.min())
    span = span if span > 0 else 1.0
    top = hi_violin
    if stars:
        y0 = hi_violin + 0.08 * span
        y1 = y0 + 0.06 * span
        add_bracket(ax, 1, 2, y0, y1, stars, FS_ANNOT)
        top = y1 + 0.06 * span

    ax.set_ylim(lo - 0.05 * span, top + 0.15 * span)
    ax.set_xlim(0.5, 2.5)
    return fig


def build_review_grid(genes, name):
    """Stack the saved panel TIFFs side by side for visual review."""
    imgs = [Image.open(OUT / f"{g}_male.tiff") for g in genes]
    h = max(im.height for im in imgs)
    w = sum(im.width for im in imgs)
    grid = Image.new("RGB", (w, h), "white")
    x = 0
    for im in imgs:
        grid.paste(im, (x, 0))
        x += im.width
    prev = OUT / "_previews"
    prev.mkdir(exist_ok=True)
    grid.save(prev / f"{name}_male_grid.png")
    print(f"[OK] {name}: review grid {grid.width}x{grid.height}px")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for gene in ALL_GENES:
        fig = violin_panel(gene)
        save_panel(fig, OUT / f"{gene}_male", vector=False, tiff=True, preview=False)
        print(f"[OK] {gene} male violin")

    build_review_grid(PAIN_MARKER_GENES, "pain_markers")
    build_review_grid(METHYLATION_MAIN_GENES, "methylation_main")
    build_review_grid(METHYLATION_OTHER_GENES, "methylation_other")
    build_review_grid(ESTROGEN_GENES, "estrogen_genes")


if __name__ == "__main__":
    main()
