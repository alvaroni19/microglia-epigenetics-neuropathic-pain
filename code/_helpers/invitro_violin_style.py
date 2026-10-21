# -*- coding: utf-8 -*-
"""Shared helpers for the in-vitro (microglia) violin replicas.

Colours sampled from the original GraphPad Prism graphs (in vitro):
male = light red/pink, female = teal.
"""
from __future__ import annotations

import sys
from pathlib import Path

import fitz

# Violin fill colours: unified Okabe-Ito sex palette (Male = vermillion,
# Female = blue) as in the six-group Csf1r violin.
GP_MALE = "#D55E00"
GP_FEMALE = "#0072B2"

# RGB of the ORIGINAL GraphPad PDF fills (used only to locate violins/brackets
# when parsing the source PDF -- must stay as the PDF's colours, not the new ones).
MALE_RGB = (1.0, 0.502, 0.502)
FEMALE_RGB = (0.059, 0.6, 0.698)
PANEL_W_MM = 34.0


def _close(c1, c2, tol=0.01):
    return all(abs(a - b) < tol for a, b in zip(c1, c2))


def violin_centers(pdf_path, page_num, color_rgb):
    """Return sorted list of x-centers of the violin-fill rects of the given
    colour on a page, excluding small colour-key legend swatches (height<8)."""
    doc = fitz.open(pdf_path)
    try:
        page = doc[page_num - 1]
        centers = []
        for dr in page.get_drawings():
            f = dr.get("fill")
            if f and _close(f, color_rgb):
                r = dr["rect"]
                if r.x0 > 200:
                    continue  # legend swatch (drawn outside the plot area)
                centers.append((r.x0 + r.x1) / 2)
        return sorted(centers)
    finally:
        doc.close()


def match_brackets_to_groups(brackets, group_centers, tol=3.5):
    """Given a list of (x0, x1, stars) brackets and an ordered list of
    (label, x_center) for each violin group, return a list of
    (label_a, label_b, stars) by matching bracket edges to nearest centers."""
    out = []
    for x0, x1, stars in brackets:
        a = min(group_centers, key=lambda gc: abs(gc[1] - x0))
        b = min(group_centers, key=lambda gc: abs(gc[1] - x1))
        if abs(a[1] - x0) > tol or abs(b[1] - x1) > tol:
            continue
        if a[0] == b[0]:
            continue
        out.append((a[0], b[0], stars))
    return out
