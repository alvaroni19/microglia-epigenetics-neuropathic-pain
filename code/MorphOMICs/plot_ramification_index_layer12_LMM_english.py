#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot layer-12 ramification indices using animal-level mixed-model tests.

Control is grey, CCI is raspberry, cells are points, and group means are
joined by a neutral line.
Significance is obtained from a linear mixed model with animal as a random
effect, avoiding cell-level pseudoreplication.
"""
from pathlib import Path
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR  = Path(__file__).resolve().parent
BUNDLE_ROOT = SCRIPT_DIR.parent.parent
sys.path.append(str(BUNDLE_ROOT / "code" / "_helpers"))
import nature_style as ns  # noqa: E402
from violin_style import draw_violin  # noqa: E402
ns.set_style()
# Fixed publication font sizes used throughout this panel.
PUB = {"star": ns.FS_ANNOT, "title": ns.FS_TITLE,
       "axis_text": ns.FS_TICK, "axis_title": ns.FS_AXIS}

PERCELL_GENERATED = BUNDLE_ROOT / "generated_outputs" / "MorphOMICs" / "cam_sholl_per_cell.csv"
PERCELL_ARCHIVED = BUNDLE_ROOT / "data" / "MorphOMICs" / "reports_cam" / "cam_sholl_per_cell.csv"
PERCELL = PERCELL_GENERATED if PERCELL_GENERATED.exists() else PERCELL_ARCHIVED
OUT_DIR = BUNDLE_ROOT / "generated_outputs" / "MorphOMICs"
LMM_CSV = OUT_DIR / "cam_ramification_index_layer12_animal_and_mixedmodel.csv"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RI = "ShollRamificationIndex"
LAYER = 12
CICLO_ORDER = ["Macho", "Proestro", "Estro", "Diestro"]
CYCLE_EN = {"Proestro": "Proestrus", "Estro": "Estrus", "Diestro": "Diestrus", "Macho": "Male"}
COL_MEAN = "#222222"   # Neutral connector between group means.
COND_COLORS = {"Control": "#999999", "CCI": "#CC79A7"}  # unified condition palette

# Linear mixed-model estimates.
def fit_mixed_models():
    import statsmodels.formula.api as smf
    dd = pd.read_csv(PERCELL)
    dd = dd[dd["Capa"] == LAYER].dropna(subset=[RI])
    p_values = {}
    rows = []
    for cycle, sub in dd.groupby("Ciclo"):
        sub = sub.copy()
        sub["Condicion"] = pd.Categorical(sub["Condicion"], ["Control", "CCI"])
        formula = "ShollRamificationIndex ~ C(Condicion)"
        if sub["Lado"].nunique() > 1:
            formula += " + C(Lado)"
        model = None
        optimizer = None
        for candidate in ("bfgs", "powell"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fitted = smf.mixedlm(formula, sub, groups=sub["AnimalID"]).fit(
                    reml=True, method=candidate, disp=False
                )
            term = next(name for name in fitted.params.index if "Condicion" in name)
            estimates = np.array([
                fitted.params[term], fitted.bse[term], fitted.tvalues[term], fitted.pvalues[term]
            ], dtype=float)
            if fitted.converged and np.isfinite(estimates).all():
                model = fitted
                optimizer = candidate
                break
        if model is None:
            raise RuntimeError(f"Mixed model failed for cycle {cycle!r}.")
        p_values[cycle] = float(model.pvalues[term])
        rows.append({
            "cycle": CYCLE_EN.get(cycle, cycle),
            "n_cells": len(sub),
            "n_animals": sub["AnimalID"].nunique(),
            "comparison": "CCI vs Control",
            "beta_CCI": float(model.params[term]),
            "standard_error": float(model.bse[term]),
            "z_value": float(model.tvalues[term]),
            "p_value": float(model.pvalues[term]),
            "random_intercept_variance": float(model.cov_re.iloc[0, 0]),
            "boundary_random_effect": bool(model.cov_re.iloc[0, 0] < 1e-6),
            "converged": bool(model.converged),
            "optimizer": optimizer,
            "estimation": "REML",
        })
    pd.DataFrame(rows).to_csv(LMM_CSV, index=False)
    return p_values

def stars(p):
    if p is None or np.isnan(p): return "ns"
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

df = pd.read_csv(PERCELL)
df = df[df["Capa"] == LAYER].dropna(subset=[RI]).copy()
df["Condicion"] = df["Condicion"].astype(str)
df["Ciclo"] = df["Ciclo"].astype(str)
cycles = [c for c in CICLO_ORDER if c in set(df["Ciclo"])]
pvals = fit_mixed_models()

# Physical size matched to the PCR violin panels (each cycle ~26 mm wide, ~42 mm
# tall) so F6 violins read identically to Fig 4B in shape, height and width.
fig, axes = plt.subplots(1, len(cycles),
                         figsize=(ns.mm2inch(26 * len(cycles) + 12), ns.mm2inch(42)),
                         sharey=True)
if len(cycles) == 1:
    axes = [axes]
fig.patch.set_facecolor("white")
fig.subplots_adjust(left=0.085, right=0.995, top=0.86, bottom=0.11, wspace=0.06)
rng = np.random.default_rng(7)
groups = ["Control", "CCI"]

for k, (ax, cycle) in enumerate(zip(axes, cycles)):
    sub = df[df["Ciclo"] == cycle]
    means, data_by_g = [], []
    for g in groups:
        vals = sub.loc[sub["Condicion"] == g, RI].to_numpy(dtype=float)
        data_by_g.append(vals)
        means.append(np.mean(vals) if len(vals) else np.nan)
    # F4B-style tapered violins (same draw_violin kernel/width as Fig 4B)
    for x, g, vals in zip([1, 2], groups, data_by_g):
        if len(vals):
            draw_violin(ax, vals, x, COND_COLORS[g], width=0.75, edgecolor="black", lw=0.5)
    # Cell-level observations.
    for j, vals in enumerate(data_by_g, start=1):
        xj = j + rng.normal(0, 0.055, size=len(vals))
        ax.scatter(xj, vals, s=8, color="black", alpha=0.30, zorder=3, linewidths=0)
    # Connector between group means.
    ax.plot([1, 2], means, "-", color=COL_MEAN, linewidth=2.6, zorder=4)
    ax.scatter([1, 2], means, s=95, color=COL_MEAN, edgecolor="white",
               linewidth=1.0, zorder=5)
    # Significance bracket based only on the mixed model.
    st = stars(pvals.get(cycle, np.nan))
    if st != "ns":
        ymax = max(np.nanmax(data_by_g[0]), np.nanmax(data_by_g[1]))
        ybar = ymax + 0.35
        ax.plot([1, 1, 2, 2], [ybar, ybar + 0.2, ybar + 0.2, ybar],
                color="black", linewidth=1.6)
        ax.text(1.5, ybar + 0.18, st, ha="center", va="bottom", fontsize=PUB["star"] + 3)
    ax.set_title(CYCLE_EN.get(cycle, cycle), fontsize=PUB["title"], fontweight="bold")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Control", "CCI"], fontsize=PUB["axis_text"], rotation=30, ha="right")
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(0.4, 2.6)
    if k > 0:
        ax.tick_params(axis="y", left=False)

axes[0].set_ylabel("Ramification index", fontsize=PUB["axis_title"], fontweight="bold")
axes[0].set_ylim(2, 9.2)

png = OUT_DIR / "cam_ramification_index_layer12_LMM_by_cycle_english.png"
tif = OUT_DIR / "cam_ramification_index_layer12_LMM_by_cycle_english.tiff"
svg = OUT_DIR / "cam_ramification_index_layer12_LMM_by_cycle_english.svg"
plt.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(tif, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(svg, bbox_inches="tight", facecolor="white")  # editable SVG/PPTX, white bg
plt.close("all")
print("[OK] Mixed-model figure saved:")
print("  ", png)
print("  ", tif)
print("Mixed-model p values:", {CYCLE_EN.get(k, k): round(v, 4) for k, v in pvals.items()})
