# -*- coding: utf-8 -*-
"""Lookup helpers for the qPCR violin significance annotations, reading the
fully Python-recomputed contrast tables in generated_outputs/PCR_260514/
(produced by code/PCR_260514/recompute_all_python_stats.py) instead of
GraphPad PDF/CSV exports.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).absolute().parents[2]
OUT_NATIVE = ROOT / "generated_outputs" / "PCR_260514"
OUT = (
    Path("\\\\?\\" + str(OUT_NATIVE.absolute()))
    if sys.platform == "win32" and not str(OUT_NATIVE).startswith("\\\\?\\")
    else OUT_NATIVE
)

_VIVO = None
_INVITRO = None

def _load_vivo():
    global _VIVO
    if _VIVO is None:
        files = [
            "invivo_pain_microglia_markers_stats.csv",
            "invivo_methylation_regulators_stats.csv",
            "invivo_estrogen_genes_stats.csv",
        ]
        _VIVO = pd.concat([pd.read_csv(OUT / f) for f in files], ignore_index=True)
    return _VIVO


def _load_invitro():
    global _INVITRO
    if _INVITRO is None:
        files = [
            "invitro_csf1_by_sex_methylation_regulators_stats.csv",
            "invitro_methylation_regulators_stats.csv",
            "invitro_estrogen_genes_stats.csv",
            "invitro_csf1r_stats.csv",
        ]
        _INVITRO = pd.concat([pd.read_csv(OUT / f) for f in files], ignore_index=True)
    return _INVITRO


def vivo_stars(gene: str, comparison: str):
    """Significance stars (e.g. '**') for an in-vivo CCI-vs-control contrast,
    or None if not significant."""
    df = _load_vivo()
    rows = df[(df.gene == gene) & (df.comparison == comparison)]
    if rows.empty:
        return None
    s = rows.iloc[0]["stars"]
    return s if isinstance(s, str) else None


def invitro_stars(gene: str, comparison: str):
    """Significance stars for an in-vitro CSF1-vs-control or
    CSF1+E2-vs-CSF1 contrast, or None if not significant."""
    df = _load_invitro()
    rows = df[(df.gene == gene) & (df.comparison == comparison)]
    if rows.empty:
        return None
    s = rows.iloc[0]["stars"]
    return s if isinstance(s, str) else None
