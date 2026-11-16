# -*- coding: utf-8 -*-
"""
Unified colour-blind-safe palette + categorical ordering for the F1-F6 refactor.

Single source of truth so every figure is internally consistent and orthogonal
biological variables never collide. Based on the Okabe-Ito set (max contrast for
deuteranopia / protanopia).

Variable families (each kept visually separable from the others):
  * Sex            -> Male = vermillion, Female = blue.   (NOTE: this REVERSES the
                      old style mapping where Male was blue.)
  * Estrous phase  -> nested within Female (sky / green / yellow).
  * Condition      -> Control vs CCI = grey vs reddish-purple.
  * Baseline/Test  -> a SEPARATE neutral pair (light-grey vs teal) so the F3
                      baseline-vs-test contrast never reads as the F6
                      control-vs-CCI contrast.
  * Methylation    -> hypo (blue, down) / hyper (red, down) for the genomic
                      track only (field convention; a different semantic context
                      from the sex blue/red, never co-plotted with sex points).

Ordering rule (global): MALES FIRST, then females; within females
proestrus -> estrus -> diestrus.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

# Reuse the shared geometry / rcParams / saving so panels stay 1:1 mm and
# SVG export (svg.fonttype='none', transparent) is inherited.
try:  # works both as a package import and a plain-path import
    from . import journal_style as ns
except Exception:  # pragma: no cover - direct script execution
    import journal_style as ns  # type: ignore

# ----------------------------------------------------------------------------- sex
SEX = {
    "Male": "#D55E00", "Female": "#0072B2",
    "Male": "#D55E00", "Female": "#0072B2", "M": "#D55E00", "F": "#0072B2",
}

# ----------------------------------------------------------------------------- estrous phase (nested in Female)
PHASE = {
    "Proestrus": "#56B4E9", "Estrus": "#009E73", "Diestrus": "#F0E442",
    "Proestrus": "#56B4E9", "Estrus": "#009E73", "Diestrus": "#F0E442",
}
# Yellow needs a dark edge to read on white; everything else edges black.
PHASE_EDGE = {
    "Diestrus": "#7A6E00", "Diestrus": "#7A6E00",
    "Proestrus": "#000000", "Estrus": "#000000",
    "Proestrus": "#000000", "Estrus": "#000000",
}

# ----------------------------------------------------------------------------- condition
CONDITION = {"Control": "#999999", "CCI": "#CC79A7", "CCI-IoN": "#CC79A7"}

# ----------------------------------------------------------------------------- baseline vs test (F3) -- intentionally distinct from CONDITION
BASELINE_TEST = {"Baseline": "#CFCFCF", "Test": "#44AA99"}

# ----------------------------------------------------------------------------- sex x treatment, in-vitro MeDIP (two-way violins, CSF1 bubble, Csf1r violin)
# Hue encodes sex; the CSF1 level is a lighter tint of the same hue (50% toward
# white) so treatment never borrows a phase or condition colour.
GROUP_SEX_TREAT = {
    "M-C": "#D55E00", "M-CSF1": "#EAAE7F",
    "F-C": "#0072B2", "F-CSF1": "#7FB8D8",
    # English / legacy aliases mapped onto the same four colours
    "MC": "#D55E00", "MCSF1": "#EAAE7F", "FC": "#0072B2", "FCSF1": "#7FB8D8",
    "Male control": "#D55E00", "Male CSF1": "#EAAE7F",
    "Female control": "#0072B2", "Female CSF1": "#7FB8D8",
}

# ----------------------------------------------------------------------------- methylation direction (genomic track F1C / F2C)
METHYL = {
    "Hypo": "#2166AC", "Hipo": "#2166AC",   # blue, density drawn downward
    "Hyper": "#B2182B", "Hiper": "#B2182B",  # red, density drawn upward
    "NO": "#BDBDBD",
}
ERE_MARKER = "#222222"   # neutral dark square for ERE coordinates (recolour in SVG if desired)

# ----------------------------------------------------------------------------- GO biological themes (F1E) -- Okabe-Ito categorical
THEME = {
    "Morphology": "#009E73",
    "Inflammation / neuroinflammation": "#D55E00",
    "Sex hormones / estrogen": "#CC79A7",
    "Mechanical stimulus / pain": "#E69F00",
    "Other": "#999999",
    "None": "#D9D9D9",
}

# ----------------------------------------------------------------------------- canonical ordering
SEX_ORDER = ["Male", "Female"]
PHASE_ORDER = ["Proestrus", "Estrus", "Diestrus"]
# Combined: males first, then female-pool, then the three phases.
GROUP_ORDER = ["Male", "Female", "Proestrus", "Estrus", "Diestrus"]
GROUP_SEX_TREAT_ORDER = ["M-C", "M-CSF1", "F-C", "F-CSF1"]

_PHASE_CANON = {"proestrus": "Proestrus", "estrus": "Estrus", "diestrus": "Diestrus",
                "proestrus": "Proestrus", "estrus": "Estrus", "diestrus": "Diestrus"}
_SEX_CANON = {"male": "Male", "male": "Male", "m": "Male",
              "female": "Female", "female": "Female", "f": "Female"}


def _rank(label: str) -> int:
    """Sort key implementing 'males first, then females (pro->est->die)'."""
    key = str(label).strip().lower()
    canon = _SEX_CANON.get(key) or _PHASE_CANON.get(key) or str(label)
    try:
        return GROUP_ORDER.index(canon)
    except ValueError:
        return len(GROUP_ORDER)  # unknown labels go last, stable


def order_categories(values):
    """Return `values` sorted by the canonical order (males first)."""
    return sorted(values, key=_rank)


def color_for(label: str, family: str = "auto") -> str:
    """Look up a colour for a category. family in
    {'auto','sex','phase','condition','baseline_test','group','methyl','theme'}."""
    tables = {
        "sex": SEX, "phase": PHASE, "condition": CONDITION,
        "baseline_test": BASELINE_TEST, "group": GROUP_SEX_TREAT,
        "methyl": METHYL, "theme": THEME,
    }
    if family != "auto":
        return tables[family].get(label, "#000000")
    for t in (SEX, PHASE, CONDITION, BASELINE_TEST, GROUP_SEX_TREAT, METHYL, THEME):
        if label in t:
            return t[label]
    return "#000000"


def edge_for(label: str) -> str:
    """Edge colour for a fill (dark edge on yellow diestrus, else black)."""
    return PHASE_EDGE.get(label, "#000000")


# ----------------------------------------------------------------------------- styling entry point
def set_style():
    """Apply the shared rcParams, then nothing else to override (palette is
    accessed explicitly via the dicts above, not via the prop_cycle)."""
    ns.set_style()


# Re-export the shared tapered-KDE violin so F4B and F6C are drawn by the SAME
# kernel/bandwidth/width (the spec requires F6C violins identical to F4B).
try:
    from .qpcr_violin_style import draw_violin
except Exception:  # pragma: no cover
    from qpcr_violin_style import draw_violin  # type: ignore


def draw_violin_unified(ax, data, x, label, *, family="auto", width=0.7,
                        lw=0.4, zorder=2):
    """Tapered violin filled with the unified palette colour for `label`."""
    return draw_violin(ax, data, x, color_for(label, family), width=width,
                       edgecolor=edge_for(label), lw=lw, zorder=zorder)


# Convenience re-exports for geometry / saving.
fig_mm = ns.fig_mm
new_fig_mm = ns.new_fig_mm
save_panel = ns.save_panel
gene_italic = ns.gene_italic
MM = ns.MM
WIDTH_1COL_MM = ns.WIDTH_1COL_MM
WIDTH_2COL_MM = ns.WIDTH_2COL_MM


def tnf_alpha(label: str) -> str:
    """Normalise any TNF spelling to the required 'TNF-α' display label."""
    s = str(label)
    for variant in ("TNFa", "TNFA", "TNFα", "TNF-a", "Tnf", "TNF"):
        if s.strip() == variant:
            return "TNF-α"
    return s
