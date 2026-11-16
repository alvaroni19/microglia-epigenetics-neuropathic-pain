# -*- coding: utf-8 -*-
"""Compact violin plots for in-vitro qPCR fold-change data
(2-way: Control/CSF1 x Male/Female primary microglia).

4 violins per gene, coloured by GraphPad's male (salmon) / female (teal)
palette as in the original GraphPad Prism graphs. Built from
`invitro_qpcr_relative_expression.xlsx` ("<Gene> 2 way" sheets).

Significance brackets show the CSF1-vs-control contrast per sex, recomputed
in Python (two-way ANOVA + Bonferroni cell-mean comparisons; see
code/PCR_260514/recompute_all_python_stats.py).

Output: generated_outputs/PCR_violins/invitro_2way/<Gene>_2way.tiff (+ a
review grid).
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
OUT = ROOT / "generated_outputs" / "PCR_violins" / "invitro_2way"

PANEL_W_MM = 34.0
PANEL_H_MM = 36.0

GROUPS = ["control_male", "control_female", "csf1_male", "csf1_female"]
GROUP_COLORS = {g: (GP_FEMALE if g.endswith("female") else GP_MALE) for g in GROUPS}
# x-axis groups as M-C, F-C, M-CSF1, F-CSF1.
GROUP_LABELS = {
    "control_male": "M-C", "control_female": "F-C",
    "csf1_male": "M-CSF1", "csf1_female": "F-CSF1",
}

SHEET_NAME_OVERRIDES_2WAY = {
    "Dnmt3a": " Dnmt3a 2 way",
}

METHYLATION_GENES = ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"]


def load_data(gene):
    sheet_name = SHEET_NAME_OVERRIDES_2WAY.get(gene, f"{gene} 2 way")
    ws = openpyxl.load_workbook(XLSX, data_only=True)[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(h).strip() if h is not None else "" for h in rows[0]]
    male_cols = [i for i, h in enumerate(header) if h.startswith("Male")]
    female_cols = [i for i, h in enumerate(header) if h.startswith("Female")]

    data = {}
    for row in rows[1:]:
        rowtitle = str(row[0]).strip().lower() if row[0] is not None else ""
        if rowtitle not in ("control", "csf1"):
            continue
        male_vals = [row[i] for i in male_cols if row[i] is not None]
        female_vals = [row[i] for i in female_cols if row[i] is not None]
        data[f"{rowtitle}_male"] = np.array(male_vals, dtype=float)
        data[f"{rowtitle}_female"] = np.array(female_vals, dtype=float)
    return data


def load_brackets(gene):
    """CSF1-vs-control brackets per sex, recomputed in Python."""
    brackets = []
    for sex, group in (("Male", "male"), ("Female", "female")):
        stars = invitro_stars(gene, f"{sex} CSF1 vs control")
        if stars:
            brackets.append((f"control_{group}", f"csf1_{group}", stars))
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
    imgs = [Image.open(OUT / f"{g}_2way.tiff") for g in genes]
    h = max(im.height for im in imgs)
    w = sum(im.width for im in imgs)
    grid = Image.new("RGB", (w, h), "white")
    x = 0
    for im in imgs:
        grid.paste(im, (x, 0))
        x += im.width
    prev = OUT / "_previews"
    prev.mkdir(exist_ok=True)
    grid.save(prev / f"{name}_2way_grid.png")
    print(f"[OK] {name}: review grid {grid.width}x{grid.height}px")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for gene in METHYLATION_GENES:
        fig = violin_panel(gene)
        save_panel(fig, OUT / f"{gene}_2way", vector=False, tiff=True, preview=False)
        print(f"[OK] {gene} in-vitro 2-way violin")

    build_review_grid(METHYLATION_GENES, "methylation_regulators")


if __name__ == "__main__":
    main()
