# -*- coding: utf-8 -*-
"""
In vivo multivariate panels (delta estrogenic score vs delta pain; marker-pain
Spearman bars) at true mm size.
Reads sources regenerated from final Delta-Ct and behavioural inputs.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).absolute().parents[2]
sys.path.append(str(ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402
import unified_palette as up  # noqa: E402

OUT = ROOT / "generated_outputs" / "PCR_MFA"
SRC = OUT


def build_delta_estrogenic_score_vs_delta_pain():
    """Phase-matched CCI-minus-control changes, female animals."""
    ns.set_style()
    df = pd.read_csv(SRC / "invivo_estrogenic_score_delta_pain_per_animal.csv")
    fig, ax = ns.fig_mm(80, 62)
    x = df["delta_estrogenic_score"].to_numpy(dtype=float)
    y = df["delta_pain_016"].to_numpy(dtype=float)
    if x.size >= 3 and np.unique(x).size >= 2:
        slope, intercept = np.polyfit(x, y, deg=1)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, color="#202020", lw=1.0, zorder=1)
    # Unified colour-blind palette (Male=vermillion, phases sky/green/yellow).
    for phase, sub in df.groupby("phase", sort=False):
        ax.scatter(sub["delta_estrogenic_score"], sub["delta_pain_016"],
                   s=14, color=up.color_for(phase), edgecolor="white",
                   linewidth=0.4, label=phase, alpha=0.95, zorder=3)
    raw = stats.spearmanr(df["delta_estrogenic_score"], df["delta_pain_016"])
    ax.set_xlabel("Δ estrogenic expression score", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.set_ylabel("Δ mechanical response (0.16 g)", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.tick_params(labelsize=ns.FS_TICK)
    ax.set_title("Animal-level association", fontsize=ns.FS_TITLE, fontweight="bold")
    ax.annotate(
        f"Spearman rho={raw.statistic:.2f}, p={raw.pvalue:.3f}, n={len(df)}",
        xy=(0.02, 0.98), xycoords="axes fraction", ha="left", va="top",
        fontsize=ns.FS_SMALL
    )
    fig.legend(title="Phase", frameon=False, loc="center left",
               bbox_to_anchor=(1.0, 0.5), bbox_transform=ax.transAxes,
               fontsize=ns.FS_LEGEND, title_fontsize=ns.FS_LEGEND_TITLE)
    ns.save_panel(fig, OUT / "invivo_delta_estrogenic_score_vs_delta_pain", vector=True, tiff=True)
    print("[OK] delta estrogenic score vs delta pain (CCI minus phase-matched control medians)")


def build_marker_pain_spearman_bars():
    """Spearman correlation of molecular markers with pain. Larger,
    italic gene labels on the y-axis."""
    ns.set_style()
    df = pd.read_csv(SRC / "invivo_marker_pain_spearman.csv").copy()
    forbidden = {"ar", "hsd17b3", "pkcg", "prkcg", "pgr", "p4r", "p2rx4", "p2x4"}
    present = forbidden.intersection(
        set(df["gene_label"].astype(str).str.lower())
    )
    if present:
        raise ValueError(f"Forbidden legacy markers in final correlation source: {present}")

    def _disp(g):
        s = ns.gene_symbol(str(g))
        s = (s[0].upper() + s[1:].lower()) if s else s
        return {"Tnf": "Tnfα"}.get(s, s)

    df["lab"] = df["gene_label"].map(_disp)   # Title-case italic + Tnfα + aliases
    df = df.sort_values("rho")
    y = np.arange(len(df))
    sig = df["q_BH"] < 0.05
    colors = np.where(df["rho"] > 0, ns.METHYL_COLORS["Hyper"], ns.METHYL_COLORS["Hypo"])
    fig, ax = ns.fig_mm(84, 90)
    ax.hlines(y, 0, df["rho"], color="#BBBBBB", lw=0.6, zorder=1)
    ax.scatter(df["rho"], y, s=np.where(sig, 22, 12), c=colors,
               edgecolor="#333333", linewidth=0.4, zorder=3)
    for yi, (_, r) in zip(y, df.iterrows()):
        if r["q_BH"] < 0.05:
            ax.text(r["rho"] + (0.03 if r["rho"] >= 0 else -0.03), yi, "*",
                    ha="left" if r["rho"] >= 0 else "right", va="center", fontsize=ns.FS_AXIS)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(df["lab"], fontsize=7.0, fontstyle="italic")  # larger gene labels
    ax.set_xlabel(r"Spearman correlation $\rho$", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.set_xlim(-1, 1)
    ax.tick_params(axis="x", labelsize=ns.FS_TICK)
    ax.tick_params(axis="y", length=0)
    # Internal title removed (goes in the figure caption per global rule).
    ns.save_panel(fig, OUT / "invivo_marker_pain_spearman_bars", vector=True, tiff=True)
    print("[OK] marker-pain Spearman bars (larger italic gene labels)")


_P = {"delta_estrogenic_vs_pain": build_delta_estrogenic_score_vs_delta_pain, "marker_pain_spearman": build_marker_pain_spearman_bars}
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--only", default="all")
    a = ap.parse_args()
    for fn in (_P.values() if a.only == "all" else [_P[a.only]]):
        fn()
