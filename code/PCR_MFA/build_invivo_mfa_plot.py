# -*- coding: utf-8 -*-
"""
MFA individuals plot, re-plotted in Python (editable SVG/PPTX) from the
coordinates exported by build_invivo_mfa.R.

The proestrus scatter points AND its 95% confidence ellipse are drawn
ON TOP of every other phase (highest z-order). Unified Okabe-Ito phase palette;
condition encoded by marker shape (Control = circle, CCI-IoN = triangle).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
import unified_palette as up  # noqa: E402
ns = up.ns

up.set_style()
OUT = ROOT / "generated_outputs" / "PCR_MFA"

# draw order: everything else first, PROESTRUS LAST (on top)
DRAW_ORDER = ["Male", "Estrus", "Diestrus", "Proestrus"]
COLOR = {"Male": up.SEX["Male"], "Proestrus": up.PHASE["Proestrus"],
         "Estrus": up.PHASE["Estrus"], "Diestrus": up.PHASE["Diestrus"]}
EDGE = {"Diestrus": "#7A6E00"}
SHAPE = {"Control": "o", "CCI": "^"}


def cov_ellipse(ax, x, y, color, zorder):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    s = np.sqrt(chi2.ppf(0.95, 2))
    w, h = 2 * s * np.sqrt(np.maximum(vals, 1e-9))
    ax.add_patch(Ellipse((x.mean(), y.mean()), w, h, angle=theta, facecolor=color,
                         edgecolor=color, lw=0.5, alpha=0.18, zorder=zorder))
    ax.add_patch(Ellipse((x.mean(), y.mean()), w, h, angle=theta, facecolor="none",
                         edgecolor=color, lw=0.6, zorder=zorder + 0.1))


def main():
    df = pd.read_csv(OUT / "invivo_mfa_coords.csv")
    eig = pd.read_csv(OUT / "invivo_mfa_eig.csv")
    var1, var2 = float(eig["var1"][0]), float(eig["var2"][0])

    fig, ax = up.fig_mm(86, 64)
    ax.axvline(0, ls="--", color="gray", lw=0.3, zorder=1)
    ax.axhline(0, ls="--", color="gray", lw=0.3, zorder=1)

    for i, phase in enumerate(DRAW_ORDER):
        sub = df[df["phase"] == phase]
        if sub.empty:
            continue
        z = 3 + i * 2  # proestrus (last) gets the highest z-order
        col = COLOR[phase]
        cov_ellipse(ax, sub["Dim1"].to_numpy(), sub["Dim2"].to_numpy(), col, zorder=z)
        for cond, mk in SHAPE.items():
            sub_c = sub[sub["condition"] == cond]
            ax.scatter(sub_c["Dim1"], sub_c["Dim2"], marker=mk, s=18, facecolor=col,
                       edgecolor=EDGE.get(phase, "black"), linewidth=0.4, zorder=z + 1,
                       label=None)

    # legends: phase (colour) + condition (shape)
    from matplotlib.lines import Line2D
    phase_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR[p],
                      markeredgecolor=EDGE.get(p, "black"), markersize=5,
                      label={"Male": "Male", "Proestrus": "Proestrus", "Estrus": "Estrus",
                             "Diestrus": "Diestrus"}[p])
               for p in ["Male", "Proestrus", "Estrus", "Diestrus"]]
    cond_h = [Line2D([0], [0], marker=SHAPE[c], color="w", markerfacecolor="#999999",
                     markeredgecolor="black", markersize=5, label=lab)
              for c, lab in (("Control", "Control"), ("CCI", "CCI-IoN"))]
    leg1 = ax.legend(handles=phase_h, title="Biological phase", loc="upper left",
                     bbox_to_anchor=(1.01, 1.0), fontsize=ns.FS_LEGEND, title_fontsize=ns.FS_LEGEND_TITLE)
    ax.add_artist(leg1)
    ax.legend(handles=cond_h, title="Condition", loc="lower left",
              bbox_to_anchor=(1.01, 0.0), fontsize=ns.FS_LEGEND, title_fontsize=ns.FS_LEGEND_TITLE)

    ax.set_xlabel(f"Dimension 1 ({var1:.1f}%)", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.set_ylabel(f"Dimension 2 ({var2:.1f}%)", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.tick_params(labelsize=ns.FS_TICK)

    ns.save_panel(fig, OUT / "invivo_mfa_proestrus_top", vector=True, tiff=True)
    print("[OK] invivo_mfa_proestrus_top (proestrus drawn on top)")


if __name__ == "__main__":
    main()
