# -*- coding: utf-8 -*-
"""Compact violin plots for in-vitro qPCR fold-change data
(male/female primary microglia x control/CSF1/CSF1+E2).

6 violins per gene, coloured by GraphPad's male (salmon) / female (teal)
palette as in the original GraphPad Prism graphs. Built from
`invitro_qpcr_relative_expression.xlsx` ("<Gene> 1 way without E2" sheets).

Significance brackets show CSF1-vs-control and CSF1+E2-vs-CSF1 per sex,
recomputed in Python (two-way ANOVA + Bonferroni cell-mean comparisons; see
scripts/PCR_260514/recompute_all_python_stats.py).

Output: generated_outputs/PCR_violins/invitro/<Gene>_invitro.tiff (+ review
grids per destination figure group).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
from PIL import Image

ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
from journal_style import (  # noqa: E402
    set_style, MM, FS_TICK, FS_AXIS, FS_TITLE,
    gene_italic, save_panel,
)
from qpcr_violin_style import draw_violin, add_bracket  # noqa: E402
from invitro_violin_style import GP_MALE, GP_FEMALE  # noqa: E402
from pcr_python_stats import invitro_stars  # noqa: E402

set_style()

XLSX = (
    ROOT
    / "data"
    / "PCR_260514"
    / "invitro_qpcr_relative_expression.xlsx"
)
OUT = ROOT / "generated_outputs" / "PCR_violins" / "invitro_6group"

PANEL_W_MM = 46.0
PANEL_H_MM = 36.0

GROUPS = ["male_control", "male_CSF1", "male_CSF1E2",
          "female_control", "female_CSF1", "female_CSF1E2"]
GROUP_COLORS = {g: (GP_MALE if g.startswith("male") else GP_FEMALE) for g in GROUPS}
GROUP_XLSX_COLS = {
    "male_control": "male control", "male_CSF1": "male CSF1", "male_CSF1E2": "male CSF1+E2",
    "female_control": "female control", "female_CSF1": "female CSF1", "female_CSF1E2": "female CSF1+E2",
}
# Hyphen between sex and treatment.
GROUP_LABELS = {
    "male_control": "M-Ctrl", "male_CSF1": "M-CSF1", "male_CSF1E2": "M-CSF1+E2",
    "female_control": "F-Ctrl", "female_CSF1": "F-CSF1", "female_CSF1E2": "F-CSF1+E2",
}

METHYLATION_GENES = ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"]
ESTROGEN_GENES = ["ERa", "ERb", "GPER", "ARO"]
CSF1R_GENES = ["CSF1R"]
ALL_GENES = METHYLATION_GENES + ESTROGEN_GENES + CSF1R_GENES


SHEET_NAME_OVERRIDES_1WAY = {
    "Dnmt3a": "Dnmt 3a 1 way without E2",
    "Dnmt3b": "Dnmt 3b 1 way without E2",
}


def load_data(gene):
    sheet_name = SHEET_NAME_OVERRIDES_1WAY.get(gene, f"{gene} 1 way without E2")
    ws = openpyxl.load_workbook(XLSX, data_only=True)[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(h).strip() for h in rows[0]]
    cols = {h: i for i, h in enumerate(header)}
    data = {}
    for group in GROUPS:
        col_name = GROUP_XLSX_COLS[group]
        idx = cols[col_name]
        vals = [row[idx] for row in rows[1:] if row[idx] is not None]
        data[group] = np.array(vals, dtype=float)
    return data


def load_brackets(gene):
    """CSF1-vs-control and CSF1+E2-vs-CSF1 brackets per sex, recomputed in
    Python."""
    brackets = []
    for sex, group in (("Male", "male"), ("Female", "female")):
        star_csf1 = invitro_stars(gene, f"{sex} CSF1 vs control")
        if star_csf1:
            brackets.append((f"{group}_control", f"{group}_CSF1", star_csf1))
        star_e2 = invitro_stars(gene, f"{sex} CSF1+E2 vs CSF1")
        if star_e2:
            brackets.append((f"{group}_CSF1", f"{group}_CSF1E2", star_e2))
    return brackets


def violin_panel(gene):
    data = load_data(gene)
    brackets = load_brackets(gene)

    fig, ax = plt.subplots(figsize=(PANEL_W_MM * MM, PANEL_H_MM * MM))

    positions = {g: i + 1.0 for i, g in enumerate(GROUPS)}
    violin_lo, violin_hi = [], []
    rng = np.random.default_rng(0)
    for g in GROUPS:
        x = positions[g]
        vals = data[g]
        v_lo, v_hi = draw_violin(ax, vals, x, GROUP_COLORS[g], width=0.75)
        violin_lo.append(v_lo)
        violin_hi.append(v_hi)
        med = np.median(vals)
        ax.plot([x - 0.18, x + 0.18], [med, med], color="black", lw=1.0, zorder=5)
        jitter = x + rng.uniform(-0.07, 0.07, size=len(vals))
        ax.scatter(jitter, vals, s=2.5, color="black", zorder=6, linewidths=0)

    lo = min(violin_lo)
    hi_violin = max(violin_hi)
    span = max(v.max() for v in data.values()) - min(v.min() for v in data.values())
    span = span if span > 0 else 1.0

    # Stack brackets in tiers, ordered by the width of the span they cover
    # (narrowest first); overlapping x-ranges get separate tiers, while
    # non-overlapping brackets share a tier.
    top = hi_violin
    sig_brackets = [b for b in brackets if b[2]]
    sorted_brackets = sorted(sig_brackets, key=lambda b: abs(positions[b[1]] - positions[b[0]]))
    bracket_h = 0.05 * span
    text_gap = 0.03 * span
    tier_step = bracket_h + 0.07 * span + 0.06 * span
    tier_intervals = []
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

    ax.set_xticks(list(positions.values()))
    ax.set_xticklabels([GROUP_LABELS[g] for g in GROUPS], fontsize=FS_TICK, rotation=45, ha="right")
    ax.set_ylabel("Fold change", fontsize=FS_AXIS, fontweight="bold")
    ax.set_title(gene_italic(gene), fontsize=FS_TITLE, fontweight="bold", pad=2.0)

    ax.set_ylim(lo - 0.05 * span, top + 0.15 * span)
    ax.set_xlim(0.3, len(GROUPS) + 0.7)
    return fig


def build_review_grid(genes, name):
    imgs = [Image.open(OUT / f"{g}_invitro.tiff") for g in genes]
    h = max(im.height for im in imgs)
    w = sum(im.width for im in imgs)
    grid = Image.new("RGB", (w, h), "white")
    x = 0
    for im in imgs:
        grid.paste(im, (x, 0))
        x += im.width
    prev = OUT / "_previews"
    prev.mkdir(exist_ok=True)
    grid.save(prev / f"{name}_invitro_grid.png")
    print(f"[OK] {name}: review grid {grid.width}x{grid.height}px")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for gene in ALL_GENES:
        fig = violin_panel(gene)
        save_panel(fig, OUT / f"{gene}_invitro", vector=False, tiff=True, preview=False)
        print(f"[OK] {gene} in-vitro 6-group violin")

    build_review_grid(METHYLATION_GENES, "methylation_regulators")
    build_review_grid(ESTROGEN_GENES, "estrogen_genes")
    build_review_grid(CSF1R_GENES, "csf1r")


if __name__ == "__main__":
    main()
