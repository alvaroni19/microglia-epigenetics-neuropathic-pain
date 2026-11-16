# -*- coding: utf-8 -*-
"""Compact violin plots for qPCR fold-change data
(female, Control vs CCI x diestrus/proestrus/estrus).

6 violins per gene (Control/CCI paired within each estrous-cycle phase),
coloured by GraphPad's phase palette (same colour for the Control/CCI pair
within a phase, matching the original Prism plots). Values are reconstructed
directly from the final Delta-Ct sheets in `invivo_qpcr_delta_ct.xlsx`.

Significance brackets show the within-phase Control-vs-CCI contrasts,
recomputed in Python (two-way ANOVA + Bonferroni cell-mean comparisons; see
scripts/PCR_260514/recompute_all_python_stats.py).

Output: generated_outputs/PCR_violins/female/<Gene>_female.tiff (+ review
grids per destination figure group), plus a shared phase-colour legend.
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
    set_style, MM, FS_TICK, FS_AXIS, FS_TITLE, FS_LEGEND,
    gene_italic, save_panel,
)
from qpcr_violin_style import (  # noqa: E402
    draw_violin, add_bracket, PHASE_COLORS_GP, PHASE_EDGE_GP,
)
from pcr_python_stats import vivo_stars  # noqa: E402
from qpcr_raw import SCRIPT_TO_CANONICAL, load_invivo_long  # noqa: E402

set_style()

OUT = ROOT / "generated_outputs" / "PCR_violins" / "female"
VIVO_LONG = load_invivo_long()

PANEL_W_MM = 46.0
PANEL_H_MM = 36.0

PHASES = ["proestrus", "estrus", "diestrus"]  # order: proestrus -> estrus -> diestrus
GROUPS = ["control", "CCI"]

PAIN_MARKER_GENES = ["CSF1R", "TLR4", "TNFa", "IL10", "KCC2", "BDNF"]
METHYLATION_MAIN_GENES = ["Dnmt1", "Dnmt3a", "Gadd45g"]
METHYLATION_OTHER_GENES = ["Dnmt3b", "Gadd45a", "Gadd45b"]
ESTROGEN_GENES = ["ERa", "ERb", "GPER", "ARO"]
ALL_GENES = PAIN_MARKER_GENES + METHYLATION_MAIN_GENES + METHYLATION_OTHER_GENES + ESTROGEN_GENES


def load_female_data(gene):
    """Return {(phase, group): np.array of replicate values (NaNs dropped)}."""
    canonical = SCRIPT_TO_CANONICAL[gene]
    subset = VIVO_LONG[
        (VIVO_LONG["gene"] == canonical) & (VIVO_LONG["phase"] != "Male")
    ].copy()
    phase_names = {
        "Proestrus": "proestrus",
        "Estrus": "estrus",
        "Diestrus": "diestrus",
    }
    subset["plot_phase"] = subset["phase"].map(phase_names)
    data = {}
    for phase in PHASES:
        for condition in ("Control", "CCI"):
            values = subset.loc[
                (subset["plot_phase"] == phase)
                & (subset["condition"] == condition),
                "relative_expression",
            ].dropna()
            group_key = "control" if condition == "Control" else "CCI"
            data[(phase, group_key)] = values.to_numpy(dtype=float)
    return data


def load_brackets(gene):
    """Within-phase Control-vs-CCI brackets, recomputed in Python (two-way
    ANOVA + Bonferroni cell-mean comparisons)."""
    brackets = []
    for phase in PHASES:
        stars = vivo_stars(gene, f"{phase.capitalize()} CCI vs control")
        if stars:
            brackets.append(((phase, "control"), (phase, "CCI"), stars))
    return brackets


def violin_panel(gene):
    data = load_female_data(gene)
    brackets = load_brackets(gene)

    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM))

    positions = {}
    pos = 1.0
    xticks, xlabels = [], []
    group_centers = {}
    for phase in PHASES:
        phase_positions = []
        for group in GROUPS:
            positions[(phase, group)] = pos
            xticks.append(pos)
            xlabels.append("Control" if group == "control" else "CCI")
            phase_positions.append(pos)
            pos += 1.0
        group_centers[phase] = sum(phase_positions) / len(phase_positions)
        pos += 0.6  # gap between phases

    all_vals = []
    violin_lo, violin_hi = [], []
    for phase in PHASES:
        color = PHASE_COLORS_GP[phase]
        for group in GROUPS:
            vals = data[(phase, group)]
            all_vals.append(vals)
            x = positions[(phase, group)]
            v_lo, v_hi = draw_violin(ax, vals, x, color, width=0.75,
                                     edgecolor=PHASE_EDGE_GP.get(phase, "black"))
            violin_lo.append(v_lo)
            violin_hi.append(v_hi)
            med = np.median(vals)
            ax.plot([x - 0.18, x + 0.18], [med, med], color="black", lw=1.0, zorder=5)
            rng = np.random.default_rng(0)
            jitter = x + rng.uniform(-0.07, 0.07, size=len(vals))
            ax.scatter(jitter, vals, s=2.5, color="black", zorder=6, linewidths=0)

    lo = min(violin_lo)
    hi_violin = max(violin_hi)
    span = max(v.max() for v in all_vals) - min(v.min() for v in all_vals)
    span = span if span > 0 else 1.0

    # significance brackets, stacked in tiers so overlapping x-ranges don't
    # collide; non-overlapping brackets share a tier.
    top = hi_violin
    sig_brackets = [b for b in brackets if b[2]]
    sorted_brackets = sorted(
        sig_brackets, key=lambda b: abs(positions[b[1]] - positions[b[0]])
    )
    bracket_h = 0.05 * span
    text_gap = 0.03 * span
    tier_step = bracket_h + 0.07 * span + 0.06 * span
    tier_intervals = []  # list of lists of (x0, x1) occupied per tier
    for g1, g2, stars in sorted_brackets:
        x0, x1 = sorted((positions[g1], positions[g2]))
        margin = 0.15 * len(stars)
        cx0, cx1 = x0 - margin, x1 + margin
        tier = 0
        while tier < len(tier_intervals) and any(
            cx0 < ix1 and cx1 > ix0 for ix0, ix1 in tier_intervals[tier]
        ):
            tier += 1
        if tier == len(tier_intervals):
            tier_intervals.append([])
        tier_intervals[tier].append((cx0, cx1))
        y0 = hi_violin + 0.08 * span + tier * tier_step
        y1 = y0 + bracket_h
        add_bracket(ax, x0, x1, y0, y1, stars, FS_TICK, text_gap=text_gap)
        top = max(top, y1 + text_gap + 0.06 * span)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=FS_TICK, rotation=45, ha="right")
    ax.set_ylabel("Fold change", fontsize=FS_AXIS, fontweight="bold")
    ax.set_title(gene_italic(gene), fontsize=FS_TITLE, fontweight="bold", pad=2.0)

    ax.set_ylim(lo - 0.05 * span, top + 0.15 * span)
    ax.set_xlim(0.3, pos - 0.6 + 0.3)
    return fig


from matplotlib.legend_handler import HandlerTuple

def build_phase_legend():
    """One shared legend mapping phase -> GraphPad colour, reused across figures."""
    fig, ax = plt.subplots(figsize=(40 * MM, 8 * MM))
    ax.axis("off")
    
    handles = []
    labels = []
    
    # Add male split handle
    h_m1 = plt.Rectangle((0, 0), 0.5, 1, facecolor="#999999", edgecolor="black", linewidth=0.4)
    h_m2 = plt.Rectangle((0, 0), 0.5, 1, facecolor="#CC79A7", edgecolor="black", linewidth=0.4)
    handles.append((h_m1, h_m2))
    labels.append("Male (Control / CCI)")
    
    # Add phases
    for p in PHASES:
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=PHASE_COLORS_GP[p], edgecolor="black", linewidth=0.4))
        labels.append(p.capitalize())
        
    ax.legend(handles, labels, loc="center", ncol=4,
              frameon=False, fontsize=FS_LEGEND, handlelength=1.2, handletextpad=0.4,
              columnspacing=1.0, handler_map={tuple: HandlerTuple(ndivide=None, pad=0)})
    return fig


def build_review_grid(genes, name):
    imgs = [Image.open(OUT / f"{g}_female.tiff") for g in genes]
    h = max(im.height for im in imgs)
    w = sum(im.width for im in imgs)
    grid = Image.new("RGB", (w, h), "white")
    x = 0
    for im in imgs:
        grid.paste(im, (x, 0))
        x += im.width
    prev = OUT / "_previews"
    prev.mkdir(exist_ok=True)
    grid.save(prev / f"{name}_female_grid.png")
    print(f"[OK] {name}: review grid {grid.width}x{grid.height}px")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for gene in ALL_GENES:
        fig = violin_panel(gene)
        save_panel(fig, OUT / f"{gene}_female", vector=False, tiff=True, preview=False)
        print(f"[OK] {gene} female violin")

    legend = build_phase_legend()
    save_panel(legend, OUT / "phase_legend", vector=False, tiff=True, preview=False)
    print("[OK] phase legend")

    build_review_grid(PAIN_MARKER_GENES, "pain_markers")
    build_review_grid(METHYLATION_MAIN_GENES, "methylation_main")
    build_review_grid(METHYLATION_OTHER_GENES, "methylation_other")
    build_review_grid(ESTROGEN_GENES, "estrogen_genes")


if __name__ == "__main__":
    main()
