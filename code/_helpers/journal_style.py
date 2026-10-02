# -*- coding: utf-8 -*-
"""Shared dimensions, typography, colours and export functions for figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------- geometry
# Nature Communications physical sizes with explicit millimetre conversion.
MM_PER_INCH = 25.4
MM = 1.0 / MM_PER_INCH               # 1 mm in inches


def mm2inch(*values_mm):
    """Explicit millimetre -> inch conversion (returns a tuple, or scalar if 1)."""
    out = tuple(v / MM_PER_INCH for v in values_mm)
    return out[0] if len(out) == 1 else out


WIDTH_1COL_MM = 89.0                 # Nature single column
WIDTH_2COL_MM = 183.0                # Nature double column
WIDTH_1P5COL_MM = 120.0
MAX_HEIGHT_MM = 247.0                # Nature maximum figure height
DPI = 600                            # native render dpi (>= 300, no upscaling)
GUTTER_MM = 4.0                      # space between panels in the assembler
PAD_IN = 0.05                        # savefig tight pad (<= 0.1 in, per spec)

# ----------------------------------------------------------------------------- fonts (pt)
# Nature spec: panel letters 8 pt bold; axis titles & legends 7 pt; ticks 5 pt.
# In-plot text is floored at 7 pt.
FS_TICK = 5.0
FS_AXIS = 7.0
FS_TITLE = 7.0
FS_LEGEND = 7.0
FS_LEGEND_TITLE = 7.0
FS_ANNOT = 7.0      # gene labels, asterisks, in-plot text
FS_SMALL = 7.0      # significance captions / in-plot annotations
FS_PANEL_LETTER = 8.0   # bold panel letter (A, B, ...) drawn by the assembler

# ----------------------------------------------------------------------------- colours
GREY = "#BDBDBD"            # Control (condition axis)
CCI = "#D1495B"            # CCI  (the user's "magenta", matches approved figs)
CONTROL_DARK = "#7A7A7A"

# MeDIP in-vitro 4-group palette (sex x treatment).
GROUP_COLORS = {
    "FC": "#F4A582", "FCSF1": "#B2182B", "MC": "#92C5DE", "MCSF1": "#2166AC",
    "Female control": "#F4A582", "Female CSF": "#B2182B", "Female CSF1": "#B2182B",
    "Male control": "#92C5DE", "Male CSF": "#2166AC", "Male CSF1": "#2166AC",
}

PHASE_COLORS = {
    "Male": "#2166AC", "Proestrus": "#E7298A", "Estrus": "#F4C95D", "Diestrus": "#E8826B",
    "Macho": "#2166AC", "Proestro": "#E7298A", "Estro": "#F4C95D", "Diestro": "#E8826B",
}

# Methylation direction (separate semantic; blue/red is the field convention, kept).
METHYL_COLORS = {"Hyper": "#C0392B", "Hypo": "#2471A3", "NO": "#BDBDBD"}

# GO biological themes (categorical, ColorBrewer Dark2-ish).
THEME_COLORS = {
    "Morphology": "#1B9E77",
    "Inflammation / neuroinflammation": "#D95F02",
    "Sex hormones / estrogen": "#7570B3",
    "Mechanical stimulus / pain": "#E7298A",
    "Other": "#BDBDBD",
    "None": "#D9D9D9",
}

CONDITION_COLORS = {"Control": "#BDBDBD", "CCI": "#D1495B", "CCI-IoN": "#D1495B"}


def control_cci_cmap():
    """Sequential-diverging colormap grey(Control) -> white -> CCI, for Control/CCI heatmaps."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("control_cci", [GREY, "#FFFFFF", CCI])


def diverging_logfc_cmap():
    """Perceptual diverging map for log2FC (ColorBrewer RdBu reversed)."""
    return plt.get_cmap("RdBu_r")


# ----------------------------------------------------------------------------- rcParams
def set_style():
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": DPI,
        "pdf.fonttype": 42,          # embed TrueType (editable text in vector PDF)
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": FS_TICK,
        "axes.titlesize": FS_TITLE,
        "axes.titleweight": "bold",
        "axes.labelsize": FS_AXIS,
        "axes.labelweight": "bold",
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "legend.title_fontsize": FS_LEGEND_TITLE,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "xtick.major.pad": 1.5,
        "ytick.major.pad": 1.5,
        "axes.titlepad": 2.0,
        "axes.labelpad": 1.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 0.7,
        "patch.linewidth": 0.4,
        "legend.frameon": False,
        "legend.handlelength": 1.2,
        "legend.handletextpad": 0.4,
        "legend.columnspacing": 0.8,
        "legend.labelspacing": 0.3,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def fig_mm(width_mm: float, height_mm: float):
    """Create a figure whose physical size is exactly width_mm x height_mm."""
    return plt.subplots(figsize=(width_mm * MM, height_mm * MM))


def new_fig_mm(width_mm: float, height_mm: float):
    return plt.figure(figsize=(width_mm * MM, height_mm * MM))


# ----------------------------------------------------------------------------- gene italics
GENE_LABELS = {
    # Display alias -> official rat gene symbol (italicised)
    "KCC2": "Slc12a5", "CSF1R": "Csf1r", "Csf1r": "Csf1r", "TLR4": "Tlr4",
    "TNFa": "Tnf", "IL10": "Il10", "BDNF": "Bdnf",
    "ERa": "Esr1", "ERb": "Esr2", "GPER": "Gper1", "ARO": "Cyp19a1",
    "Dnmt1": "Dnmt1", "Dnmt3a": "Dnmt3a", "Dnmt3b": "Dnmt3b",
    "Gadd45a": "Gadd45a", "Gadd45b": "Gadd45b", "Gadd45g": "Gadd45g",
}


def gene_symbol(name: str) -> str:
    return GENE_LABELS.get(name, name)


def _title_case_symbol(s: str) -> str:
    """Proper rodent gene-symbol case: first letter upper, the rest lower
    (e.g. CSF1R -> Csf1r, KCC2 alias already mapped). Greek suffixes preserved."""
    if not s:
        return s
    return s[0].upper() + s[1:].lower()


def gene_italic(name: str) -> str:
    """Return a mathtext italic gene symbol in proper rodent Title case.

    Global rules (user spec): Title case + italic; Tnf is displayed as Tnfα.
    """
    s = _title_case_symbol(gene_symbol(name))
    if s == "Tnf":
        return "$\\mathit{Tnf}$α"
    safe = s.replace(" ", r"\ ")
    return rf"$\mathit{{{safe}}}$"


# Convenience alias matching the new theme naming.
gene_label = gene_italic


def italicize_gene_ticklabels(ax, axis="y"):
    """Set the given axis tick labels to italic font style (for gene names)."""
    labels = ax.get_yticklabels() if axis == "y" else ax.get_xticklabels()
    for t in labels:
        t.set_fontstyle("italic")


# ----------------------------------------------------------------------------- repulsion
def repel(ax, texts, **kwargs):
    """Repel a list of matplotlib Text objects (anti-overlap). No-op if empty."""
    if not texts:
        return
    try:
        from adjustText import adjust_text
    except Exception:
        return
    params = dict(
        arrowprops=dict(arrowstyle="-", color="#555555", lw=0.4),
        expand=(1.05, 1.2),
        force_text=(0.3, 0.5),
        min_arrow_len=2,
    )
    params.update(kwargs)
    adjust_text(texts, ax=ax, **params)


# ----------------------------------------------------------------------------- saving
def save_panel(fig, out_base: Path, vector: bool, tiff: bool = True,
               preview: bool = True, pad_mm: float = 0.4, svg: bool = True):
    """
    Save a panel rendered at its true physical size.
      vector=True  -> also write a .pdf (for scatter/line/scheme panels)
      tiff=True    -> write a .tiff at DPI (LZW)  (deliverable / assembly source)
      svg=True     -> write an editable .svg (svg.fonttype='none' -> text stays
                      text; transparent so panels overlay cleanly when the user
                      ungroups them in PowerPoint / Illustrator)
      preview=True -> write a small .png for visual inspection (cleaned up later)
    `out_base` is a path WITHOUT extension.
    """
    out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    # Preserve the declared physical canvas. A tight bounding box changes the
    # output dimensions and prevents 1:1 panel assembly.
    saved = []
    if tiff:
        p = out_base.with_suffix(".tiff")
        fig.savefig(p, dpi=DPI, bbox_inches=None,
                    facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
        saved.append(p)
    if vector:
        p = out_base.with_suffix(".pdf")
        fig.savefig(p, bbox_inches=None, facecolor="white")
        saved.append(p)
    if svg:
        p = out_base.with_suffix(".svg")
        fig.savefig(p, format="svg", bbox_inches=None,
                    facecolor="white", transparent=False)
        saved.append(p)
    if preview:
        prev = out_base.parent / "_previews"
        prev.mkdir(parents=True, exist_ok=True)
        p = prev / (out_base.name + ".png")
        fig.savefig(p, dpi=200, bbox_inches=None)
        saved.append(p)
    plt.close(fig)
    return saved
