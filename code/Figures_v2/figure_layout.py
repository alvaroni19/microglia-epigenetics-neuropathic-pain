"""Common style and slide layout for the figure panels.

Every panel block is drawn on a figure whose size equals its slot on the
180 x 170 mm slides of the figure PowerPoints, so the block is placed at 100 %
scale and all text keeps its nominal size.  Coordinates are in millimetres
with the origin at the top-left corner of the slide.

Style (final size): Arial; tick labels 6 pt; axis and panel titles 7 pt bold;
gene titles 7 pt italic; legend entries 6 pt and legend titles 7 pt;
significance asterisks 8 pt (black; white with a dark outline inside dark
bubbles or cells).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

MM = 1 / 25.4

# ----------------------------------------------------------------------------- style
FS_TICK = 6.0
FS_LABEL = 7.0
FS_TITLE = 7.0
FS_LEGEND = 6.0
FS_LEGEND_TITLE = 7.0
FS_STAR = 8.0
FS_NS = 6.0
LW_AXIS = 0.6
LW_BRACKET = 0.5
DOT = 3.0            # data points inside violins (points^2)
SCATTER = 14.0       # data points in scatter plots (points^2)
LEG_MS = 4.0         # legend marker size (points)
TICK_LEN = 2.0


def rc():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": FS_TICK, "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.linewidth": LW_AXIS, "xtick.major.width": LW_AXIS, "ytick.major.width": LW_AXIS,
        "xtick.major.size": TICK_LEN, "ytick.major.size": TICK_LEN, "xtick.major.pad": 1,
        "ytick.major.pad": 1, "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "axes.labelsize": FS_LABEL, "axes.labelweight": "bold", "axes.labelpad": 1,
        "axes.titlesize": FS_TITLE, "axes.titlepad": 2, "legend.fontsize": FS_LEGEND,
        "legend.title_fontsize": FS_LEGEND_TITLE, "legend.frameon": False,
        "legend.handletextpad": 0.2, "legend.labelspacing": 0.15, "legend.borderaxespad": 0,
        "legend.borderpad": 0, "legend.columnspacing": 0.8, "legend.handlelength": 1.0,
        "lines.solid_capstyle": "butt",
    })


def star_text(ax, x, y, s, dark=False, **kw):
    """Asterisks in the common style; dark=True for dark bubble/cell backgrounds."""
    t = ax.text(x, y, s, fontsize=FS_STAR if s != "n.s." else FS_NS,
                color="white" if dark else "black", zorder=8, **kw)
    if dark:
        t.set_path_effects([pe.withStroke(linewidth=1.2, foreground="#202020")])
    return t


def bracket(ax, x0, x1, y, h, s):
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color="black", lw=LW_BRACKET, zorder=7)
    star_text(ax, (x0 + x1) / 2, y + (0.15 if s != "n.s." else 1.4) * h, s, ha="center", va="bottom")


def finish(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(LW_AXIS)
    ax.tick_params(labelsize=FS_TICK, width=LW_AXIS, length=TICK_LEN, pad=1)


# ----------------------------------------------------------------------------- block canvas
class Block:
    """Figure of exactly w x h mm; axes are placed with mm coordinates from the
    top-left corner of the block."""

    def __init__(self, w, h):
        self.w, self.h = float(w), float(h)
        self.fig = plt.figure(figsize=(self.w * MM, self.h * MM))

    def ax(self, x, y, w, h, **kw):
        return self.fig.add_axes([x / self.w, 1 - (y + h) / self.h, w / self.w, h / self.h], **kw)

    def text(self, x, y, s, **kw):
        return self.fig.text(x / self.w, 1 - y / self.h, s, **kw)

    def fx(self, x):
        return x / self.w

    def fy(self, y):
        return 1 - y / self.h

    def save(self, stem):
        stem = Path(stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        for ext in ("svg", "pdf", "png"):
            self.fig.savefig(stem.with_suffix(f".{ext}"), dpi=600, facecolor="white")
        plt.close(self.fig)


# ----------------------------------------------------------------------------- cell sizes (mm)
# in vivo gene cell: male axes + female axes (3 phases), rotated Control/CCI-IoN labels
VIVO_CELL = dict(w=58.5, h=33.5, left=8.0, male_w=11.5, gap=7.5, female_w=30.0, top=3.5, ax_h=21.5)
# in vitro 6-group gene cell
VITRO_CELL = dict(w=45.0, h=38.0, left=8.0, ax_w=35.5, top=3.5, ax_h=23.5)
LEGEND_ROW = 5.0
BUBBLE_CELL = 7.5     # bubble matrices, square cell
HEAT_CELL = (11.0, 3.8)

# ----------------------------------------------------------------------------- block sizes
# Final printed size (width, height) in mm of every data panel drawn by
# build_figure_blocks.py; block ids are the output file stems.
BLOCK_SIZES = {
    "medip_dmr_directionality": (118.0, 34.0),
    "medip_go_rank_distribution": (112.0, 78.0),
    "invitro_qpcr_methylation_regulators_csf1_bubble": (55.0, 57.0),
    "medip_candidate_gene_heatmap": (55.0, 82.0),
    "invitro_qpcr_methylation_regulators_bubble": (76.0, 66.0),
    "invitro_qpcr_csf1r_violin": (76.0, 64.0),
    "medip_csf1r_locus_dmr_ere_map": (175.0, 71.0),
    "behaviour_change_score_violins": (118.0, 64.0),
    "invivo_qpcr_pain_microglia_markers": (175.5, 72.0),
    "invivo_qpcr_methylation_regulators": (117.0, 67.0),
    "invivo_marker_pain_spearman": (57.0, 58.0),
    "invivo_estrogen_genes_bubble": (82.0, 62.0),
    "invitro_estrogen_genes_bubble": (82.0, 62.0),
    "estrogenic_score_vs_pain": (70.0, 44.0),
    "invitro_estrogenic_score": (82.0, 50.0),
    "invivo_mfa_phase_condition": (76.0, 56.0),
    "invitro_estrogen_genes_pca": (82.0, 47.0),
    "morphology_ramification_index": (102.0, 34.0),
    "morphology_umap": (102.0, 30.0),
    "morphology_persistence_difference": (104.0, 52.0),
    "morphology_ramification_pain_partial_spearman": (100.0, 48.0),
    "medip_robpca": (75.0, 57.0),
    "medip_global_methylation_heatmap": (62.0, 170.0),
    "medip_dmr_gene_set_upset": (77.0, 60.0),
    "medip_candidate_window_map": (77.0, 49.0),
    "invitro_qpcr_methylation_regulators_csf1_violins": (34.0, 167.0),
    "invitro_qpcr_methylation_regulators_violins": (135.0, 81.0),
    "invivo_qpcr_gadd45ab_violins": (117.0, 38.5),
    "invivo_qpcr_estrogen_genes_violins": (117.0, 72.0),
    "invitro_qpcr_estrogen_genes_violins": (180.0, 43.0),
}
