# -*- coding: utf-8 -*-
"""Shared violin, significance and colour helpers for manuscript figures."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde, t as t_dist

# Unified condition and estrous-phase palette.
GP_CONTROL = "#999999"  # Control  (grey)
PANEL_W_MM = 34.0
GP_CCI = "#CC79A7"      # CCI      (reddish purple)
GP_DIESTRUS = "#F0E442"   # diestrus (yellow; dark edge applied when drawn)
GP_PROESTRUS = "#56B4E9"  # proestrus (sky blue)
GP_ESTRUS = "#009E73"     # estrus (bluish green)

# Canonical phase order used throughout the manuscript figures.
PHASE_COLORS_GP = {
    "proestrus": GP_PROESTRUS,
    "estrus": GP_ESTRUS,
    "diestrus": GP_DIESTRUS,
}

# Dark edge on the yellow diestrus fill so it reads on white.
PHASE_EDGE_GP = {"diestrus": "#7A6E00", "proestrus": "#000000", "estrus": "#000000"}


def draw_violin(ax, data, x, color, width=0.7, edgecolor="black", lw=0.4, zorder=2):
    """Draw a single tapered (pointed-end) violin at position x.

    Returns (lo, hi): the full vertical extent of the drawn shape (including
    the tapered tails, which extend beyond the raw data range), so callers
    can size ylim / bracket placement to avoid clipping the points.
    """
    data = np.asarray(data, dtype=float)
    if len(data) < 2 or np.allclose(data, data[0]):
        # degenerate case: draw a thin sliver so the shape doesn't error out
        y = data.mean() if len(data) else 0.0
        ax.plot([x - width / 4, x + width / 4], [y, y], color=color, lw=2, zorder=zorder)
        return y, y
    kde = gaussian_kde(data)
    bw = kde.factor * data.std(ddof=1)
    lo, hi = data.min() - 3 * bw, data.max() + 3 * bw
    ys = np.linspace(lo, hi, 200)
    dens = kde(ys)
    dens = dens / dens.max() * (width / 2)
    ax.fill_betweenx(ys, x - dens, x + dens, facecolor=color, edgecolor=edgecolor,
                      linewidth=lw, zorder=zorder)
    return lo, hi


def stars_or_none(summary: str | None):
    """Return the significance string, or None if not significant / not available."""
    if summary is None:
        return None
    summary = summary.strip()
    if summary in ("", "ns", "ns."):
        return None
    return summary


def _to_float(s):
    return float(s.strip().replace(",", "."))


def read_cell_means_pairs(path: Path):
    """Parse a Prism 'Compare cell means regardless of rows and columns'
    (Bonferroni multiple comparisons) export into
    {(group1, group2): summary_string}, keys lower-cased
    e.g. ("control:diestrus", "cci:diestrus") -> "*"."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter=";"))
    out = {}
    in_table = False
    for row in rows:
        if not row:
            continue
        label = row[0].strip()
        if "multiple comparisons test" in label.lower():
            in_table = True
            continue
        if label.startswith("Test details"):
            break
        if in_table and " vs. " in label:
            g1, g2 = [g.strip().lower() for g in label.split(" vs. ")]
            summary = row[4].strip()
            out[(g1, g2)] = summary
            out[(g2, g1)] = summary
    return out


def is_cell_means_format(path: Path) -> bool:
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if "Compare cell means" in line:
                return True
            if line.strip().startswith("Residual;"):
                return False
    return False


def read_residual_ms_df(path: Path):
    """Read 'Residual' MS and DF from a Prism 'Two-way ANOVA' CSV export.

    The ANOVA table section looks like:
        ANOVA table;SS (Type III);DF;MS;F (DFn, DFd);P value
            Residual;0,5929;19;0,03121;;
    Returns (MS, DF) as floats.
    """
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter=";"))
    for row in rows:
        if row and row[0].strip() == "Residual":
            return _to_float(row[3]), _to_float(row[2])
    raise ValueError(f"Residual row not found in {path}")


def p_to_stars(p: float):
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def sidak_posthoc_control_vs_cci(data: dict, ms_resid: float, df_resid: float, phases):
    """Sidak-corrected Control-vs-CCI t-tests within each phase, sharing the
    pooled residual variance/df from the omnibus 2-way ANOVA (matches Prism's
    "Sidak's multiple comparisons test" using the model's error term).

    `data` maps (phase, "control"/"CCI") -> np.array of replicate values.
    Returns {phase: stars_or_None}.
    """
    m = len(phases)
    out = {}
    for phase in phases:
        a = data[(phase, "control")]
        b = data[(phase, "CCI")]
        na, nb = len(a), len(b)
        se = np.sqrt(ms_resid * (1.0 / na + 1.0 / nb))
        t_val = (b.mean() - a.mean()) / se
        p = 2 * (1 - t_dist.cdf(abs(t_val), df_resid))
        p_adj = 1 - (1 - p) ** m
        out[phase] = stars_or_none(p_to_stars(p_adj))
    return out


def add_bracket(ax, x0, x1, y0, y1, text, fontsize, lw=0.6, text_gap=0.0):
    """Draw a significance bracket with centred text above it.

    The label is offset by a small fixed amount in points (not data units),
    so it sits just above the bracket line regardless of the y-axis scale.
    """
    ax.plot([x0, x0, x1, x1], [y0, y1, y1, y0], lw=lw, color="black", zorder=10)
    
    # Asterisks render high in their text box, so they need a negative offset
    # to sit close to the bracket without floating too high.
    is_stars = text and all(c == '*' for c in text.strip())
    y_offset = -3.0 if is_stars else 1.5
    
    ax.annotate(text, xy=((x0 + x1) / 2, y1), xytext=(0, y_offset),
                 textcoords="offset points",
                 ha="center", va="bottom", fontsize=fontsize, fontweight="bold", zorder=10)
