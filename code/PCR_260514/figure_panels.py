# -*- coding: utf-8 -*-
"""
PCR_260514 publication panels rebuilt at true mm size (
5-7 pt fonts, italic gene symbols, RGB perceptual palettes).

Reads from the source stats dataframes already produced by the qPCR pipeline
(generated_outputs/PCR_260514/*_bubble_stats.csv), which are the exact data
layer behind the published bubble matrices.
"""
from __future__ import annotations

import sys
from pathlib import Path
from textwrap import wrap

import numpy as np
import pandas as pd
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse, Patch, Rectangle

ROOT = Path(__file__).absolute().parents[2]
sys.path.append(str(ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402

OUT = ROOT / "generated_outputs" / "PCR_260514"
STATS = OUT  # bubble stats csv live alongside the tiffs


def _bubble(df, comps, genes, width_mm, height_mm, cbar=True, sig_caption=True,
            box_w=1.0, col_labels=None):
    """Fixed-size bubble matrix: rows=genes (italic), cols=comparisons,
    colour = log2FC (diverging), asterisks = significance.

    box_w: width (in column units, <=1) of the grey background cell. Narrower
           cells hug the bubbles without changing bubble size or
           spacing, so bubbles never overlap.
    col_labels: optional {comparison: display label} for the x tick labels.
    """
    ns.set_style()
    df = df[df["comparison"].isin(comps)].copy()
    genes = [g for g in genes if g in set(df["gene_label"])]
    y_map = {g: i for i, g in enumerate(reversed(genes))}
    x_map = {c: i for i, c in enumerate(comps)}
    lim = float(np.nanmax(np.abs(df["log2fc"].to_numpy())))
    lim = (lim if np.isfinite(lim) and lim > 0 else 1.0) * 1.05

    fig, ax = ns.fig_mm(width_mm, height_mm)
    # background cells (narrowed to box_w, centred on the bubble)
    for x in range(len(comps)):
        for y in range(len(genes)):
            ax.add_patch(Rectangle((x - box_w / 2, y - 0.5), box_w, 1, facecolor="#F7F7F7",
                                   edgecolor="#DDDDDD", lw=0.4, zorder=0))
    cmap = ns.diverging_logfc_cmap()
    norm = Normalize(vmin=-lim, vmax=lim)
    # marker diameter sized to the tighter (row) pitch so bubbles never overlap
    row_pitch_in = (height_mm * ns.MM * 0.80) / max(len(genes), 1)
    s = (0.78 * row_pitch_in * 72) ** 2   # points^2
    for _, r in df.iterrows():
        x, y = x_map[r["comparison"]], y_map[r["gene_label"]]
        col = cmap(norm(r["log2fc"])) if np.isfinite(r["log2fc"]) else (1, 1, 1, 1)
        ax.scatter(x, y, s=s, facecolor=col, edgecolor="#333333", linewidth=0.5, zorder=3)
        if isinstance(r["stars"], str) and r["stars"]:
            text = ax.text(x, y, r["stars"], ha="center", va="center", fontsize=9.0,
                           fontweight="bold", color="black", zorder=4)
            text.set_path_effects([pe.withStroke(linewidth=1.0, foreground="white")])

    ax.set_xticks(range(len(comps)))
    xlabels = [(col_labels or {}).get(c, c) for c in comps]
    ax.set_xticklabels(xlabels, fontsize=ns.FS_AXIS, rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(list(reversed(genes)), fontsize=ns.FS_AXIS, fontstyle="italic")
    ax.tick_params(length=0)
    ax.set_xlim(-0.55, len(comps) - 0.45)
    ax.set_ylim(-0.55, len(genes) - 0.45)
    for sp in ax.spines.values():
        sp.set_visible(False)

    if cbar:
        import matplotlib.pyplot as plt
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.04, aspect=12)
        cb.solids.set_rasterized(False)
        cb.set_label("log$_2$ fold change", fontsize=ns.FS_AXIS, fontweight="bold")
        cb.ax.tick_params(labelsize=ns.FS_TICK, width=0.4, length=2)
        cb.outline.set_linewidth(0.4)
    # significance caption removed — rotated x-labels would overlap it
    return fig


def _disp_gene(g):
    # Map qPCR aliases to official rodent symbols (KCC2->Slc12a5, ERa->Esr1,
    # GPER->Gper1, ARO->Cyp19a1, ...), enforce Title case, and the global
    # Tnf -> Tnfα display rule. Italic is applied at draw time.
    s = ns.gene_symbol(str(g))
    s = (s[0].upper() + s[1:].lower()) if s else s
    return {"Tnf": "Tnfα"}.get(s, s)


def build_bubble(csv, out, comps=None, genes=None, w=None, h=None, title="",
                 col_labels=None, box_w=1.0):
    df = pd.read_csv(STATS / csv)
    df["gene_label"] = df["gene_label"].map(_disp_gene)
    if comps is None:
        comps = list(dict.fromkeys(df["comparison"]))
    if genes is None:
        genes = list(dict.fromkeys(df["gene_label"]))
    else:
        genes = [_disp_gene(g) for g in genes]
    if w is None:
        w = 30 + 16 * len(comps)
    if h is None:
        h = 20 + 9 * len(genes)
    fig = _bubble(df, comps, genes, width_mm=w, height_mm=h,
                  box_w=box_w, col_labels=col_labels)
    ns.save_panel(fig, OUT / out, vector=False, tiff=True)
    print(f"[OK] {out}.tiff  ({len(genes)}x{len(comps)})")


VIVO = ["Male CCI vs control", "Diestrus CCI vs control",
        "Proestrus CCI vs control", "Estrus CCI vs control"]
VITRO = ["Male CSF1 vs control", "Male CSF1+E2 vs CSF1",
         "Female CSF1 vs control", "Female CSF1+E2 vs CSF1"]


def build_invitro_csf1_by_sex_bubble():
    """In vitro methylation regulators, male/female CSF1 vs control only.
    Uses 2-way ANOVA + Bonferroni (6 comparisons) to match the two-way violins."""
    build_bubble("invitro_csf1_by_sex_methylation_regulators_stats.csv",
                 "invitro_csf1_by_sex_methylation_regulators_bubble",
                 comps=["Male CSF1 vs control", "Female CSF1 vs control"],
                 genes=["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"],
                 col_labels={"Male CSF1 vs control": "M-CSF1 vs C",
                             "Female CSF1 vs control": "F-CSF1 vs C"},
                 box_w=0.62, w=54, h=58)


VITRO_LABELS = {
    "Male CSF1 vs control": "M-CSF1 vs C",
    "Male CSF1+E2 vs CSF1": "M-CSF1+E2 vs M-CSF1",
    "Female CSF1 vs control": "F-CSF1 vs C",
    "Female CSF1+E2 vs CSF1": "F-CSF1+E2 vs F-CSF1",
}
VIVO_LABELS = {
    "Male CCI vs control": "M-CCI vs C",
    "Diestrus CCI vs control": "Diestrus-CCI vs C",
    "Proestrus CCI vs control": "Proestrus-CCI vs C",
    "Estrus CCI vs control": "Estrus-CCI vs C",
}


def build_invitro_methylation_regulators_bubble():
    """Full methylation-regulator pairwise matrix (6 genes x 4 contrasts).
    Same changes as F1E (Title-case italic genes, narrowed boxes,
    hyphenated comparison labels)."""
    build_bubble("invitro_methylation_regulators_stats.csv",
                 "invitro_methylation_regulators_bubble", comps=VITRO,
                 genes=["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"],
                 col_labels=VITRO_LABELS, box_w=0.62)


def build_invivo_pain_markers_bubble():
    build_bubble("invivo_pain_microglia_markers_stats.csv",
                 "invivo_pain_microglia_markers_bubble", comps=VIVO)


def build_invivo_methylation_regulators_bubble():
    build_bubble("invivo_methylation_regulators_stats.csv",
                 "invivo_methylation_regulators_bubble", comps=VIVO)


ER_ARO_GENES = ["ERa", "ERb", "GPER", "ARO"]   # -> Esr1, Esr2, Gper1, Cyp19a1


def build_invivo_estrogen_genes_bubble():
    """Presentation F5A: in-vivo ER/aromatase bubble. Same treatment as F1E
    (Title-case italic genes, narrowed boxes, hyphenated comparison labels)."""
    build_bubble("invivo_estrogen_genes_stats.csv",
                 "invivo_estrogen_genes_bubble", comps=VIVO, genes=ER_ARO_GENES,
                 col_labels=VIVO_LABELS, box_w=0.62)


def build_invitro_estrogen_genes_bubble():
    """Presentation F5D: in-vitro ER/aromatase bubble. Same treatment as F1E."""
    build_bubble("invitro_estrogen_genes_stats.csv",
                 "invitro_estrogen_genes_bubble", comps=VITRO, genes=ER_ARO_GENES,
                 col_labels=VITRO_LABELS, box_w=0.62)


SCR = ROOT / "generated_outputs" / "PCR_260514"
VITRO_GROUP_COLORS = {
    "Male Control": "#BDD7E7", "Male CSF1": "#6BAED6", "Male CSF1+E2": "#2171B5",
    "Female Control": "#FCAE91", "Female CSF1": "#FB6A4A", "Female CSF1+E2": "#CB181D",
}
VIVO_GROUP_COLORS = {
    "Male Control": "#9ECAE1", "Male CCI-IoN": "#2171B5",
    "Proestrus Control": "#B2DF8A", "Proestrus CCI-IoN": "#33A02C",
    "Estrus Control": "#FDBF6F", "Estrus CCI-IoN": "#FF7F00",
    "Diestrus Control": "#CAB2D6", "Diestrus CCI-IoN": "#6A3D9A",
}
PHASE_SHAPES = {"Male": "s", "Proestrus": "o", "Estrus": "^", "Diestrus": "D",
                "Male": "s", "Proestrus": "o", "Estrus": "^", "Diestrus": "D"}
SEX_COLORS = {"Male": "#D55E00", "Female": "#0072B2"}  # unified Okabe-Ito sex
VIVO_CONDITION_SHAPES = {"Control": "o", "CCI-IoN": "^"}
VITRO_CONDITION_SHAPES = {"Control": "o", "CSF1": "s", "CSF1+E2": "^"}


def _variance(path):
    v = pd.read_csv(path)
    return {r["PC"]: r["variance_explained_pct"] for _, r in v.iterrows()}


def _add_covariance_ellipse(ax, xy, color, level=0.95):
    """Draw a descriptive normal-theory data ellipse, not a confidence interval."""
    xy = np.asarray(xy, dtype=float)
    xy = xy[np.isfinite(xy).all(axis=1)]
    if len(xy) < 3:
        return
    covariance = np.cov(xy, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if np.any(eigenvalues <= 0):
        return
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    # sqrt(chi-square quantile with 2 df); fixed value avoids an extra dependency.
    radius = 2.44774683068 if level == 0.95 else np.sqrt(-2.0 * np.log1p(-level))
    width, height = 2.0 * radius * np.sqrt(eigenvalues)
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    center = xy.mean(axis=0)
    ax.add_patch(Ellipse(
        center, width, height, angle=angle, facecolor=color, edgecolor=color,
        alpha=0.18, linewidth=0.0, zorder=1,
    ))
    ax.add_patch(Ellipse(
        center, width, height, angle=angle, facecolor="none", edgecolor=color,
        linewidth=0.55, zorder=2,
    ))


def build_invitro_estrogenic_score_and_pca():
    """Directed estrogenic score and PCA from the same in-vitro data."""
    import matplotlib.pyplot as plt
    from scipy import stats

    ns.set_style()
    sc = pd.read_csv(SCR / "invitro_estrogen_genes_pca_scores.csv")
    var = _variance(SCR / "invitro_estrogen_genes_pca_variance.csv")
    score = pd.read_csv(SCR / "invitro_estrogenic_score_per_replicate.csv")

    sex_order = ["Male", "Female"]
    treatment_order = ["Control", "CSF1", "CSF1+E2"]
    # Hyphenated comparison labels.
    contrasts = [
        ("Male", "CSF1", "Control", "M-CSF1\nvs M-C"),
        ("Male", "CSF1+E2", "CSF1", "M-CSF1+E2\nvs M-CSF1"),
        ("Female", "CSF1", "Control", "F-CSF1\nvs F-C"),
        ("Female", "CSF1+E2", "CSF1", "F-CSF1+E2\nvs F-CSF1"),
    ]

    plot_rows = []
    stat_rows = []
    for contrast_idx, (sex, numerator, denominator, label) in enumerate(contrasts):
        numerator_values = score.loc[
            (score["sex"] == sex) & (score["treatment"] == numerator), "ESI_full"
        ].dropna().to_numpy(float)
        denominator_values = score.loc[
            (score["sex"] == sex) & (score["treatment"] == denominator), "ESI_full"
        ].dropna().to_numpy(float)
        baseline = float(np.median(denominator_values))
        delta_values = numerator_values - baseline
        test = stats.mannwhitneyu(
            numerator_values, denominator_values, alternative="two-sided", method="asymptotic"
        )
        for replicate, value in enumerate(delta_values, start=1):
            plot_rows.append({
                "contrast_idx": contrast_idx,
                "contrast": label.replace("\n", " "),
                "sex": sex,
                "numerator": numerator,
                "denominator": denominator,
                "replicate": replicate,
                "directed_delta_estrogenic_score": value,
            })
        stat_rows.append({
            "contrast_idx": contrast_idx,
            "contrast": label.replace("\n", " "),
            "sex": sex,
            "numerator": numerator,
            "denominator": denominator,
            "n_numerator": len(numerator_values),
            "n_denominator": len(denominator_values),
            "median_delta": float(np.median(delta_values)),
            "mannwhitney_u": float(test.statistic),
            "p_value": float(test.pvalue),
        })

    plot_df = pd.DataFrame(plot_rows)
    stat_df = pd.DataFrame(stat_rows)
    p_values = stat_df["p_value"].to_numpy(float)
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q_values = np.empty_like(ranked)
    q_values[order] = np.minimum(ranked, 1.0)
    stat_df["q_BH"] = q_values
    stat_df.to_csv(OUT / "invitro_estrogenic_score_stats.csv", index=False)

    fig, (ax_score, ax_pca) = plt.subplots(
        1, 2,
        figsize=(176 * ns.MM, 58 * ns.MM),
        gridspec_kw={"width_ratios": [1.02, 1.18]},
    )
    fig.subplots_adjust(left=0.065, right=0.83, bottom=0.25, top=0.88, wspace=0.34)

    rng = np.random.default_rng(20260608)
    x_positions = np.arange(len(contrasts))
    for x, (sex, _, _, _) in enumerate(contrasts):
        values = plot_df.loc[
            plot_df["contrast_idx"] == x, "directed_delta_estrogenic_score"
        ].to_numpy(float)
        median = float(stat_df.loc[stat_df["contrast_idx"] == x, "median_delta"].iloc[0])
        color = SEX_COLORS[sex]
        ax_score.bar(
            x, median, width=0.62, facecolor=color, edgecolor=color,
            linewidth=0.55, alpha=0.30, zorder=2,
        )
        jitter = rng.normal(0.0, 0.045, size=len(values))
        ax_score.scatter(
            np.full(len(values), x) + jitter, values, s=15,
            marker="o", facecolor=color,
            edgecolor="black", linewidth=0.4, alpha=0.95, zorder=3,
        )

    score_min = float(plot_df["directed_delta_estrogenic_score"].min())
    score_max = float(plot_df["directed_delta_estrogenic_score"].max())
    score_span = max(score_max - score_min, 1.0)
    ax_score.set_ylim(score_min - 0.16 * score_span, score_max + 0.12 * score_span)
    ax_score.axhline(0, color="#707070", ls="--", lw=0.45, zorder=1)
    ax_score.axvline(1.5, color="#A0A0A0", ls=":", lw=0.4, zorder=1)
    ax_score.set_xticks(x_positions)
    ax_score.set_xticklabels([c[3] for c in contrasts], fontsize=ns.FS_SMALL)
    ax_score.set_ylabel("Directed delta\nestrogenic score", fontsize=ns.FS_AXIS,
                        fontweight="bold")
    ax_score.set_title("Directed estrogenic score", fontsize=ns.FS_TITLE,
                       fontweight="bold")
    ax_score.tick_params(labelsize=ns.FS_TICK)
    ax_score.grid(axis="y", color="#E2E2E2", lw=0.35)
    ax_score.set_axisbelow(True)

    for sex in sex_order:
        sub = sc[sc["sex"] == sex]
        _add_covariance_ellipse(ax_pca, sub[["PC1", "PC2"]], SEX_COLORS[sex])
    for treatment in treatment_order:
        sub = sc[sc["treatment"] == treatment]
        for sex in sex_order:
            group = sub[sub["sex"] == sex]
            if len(group):
                ax_pca.scatter(
                    group["PC1"], group["PC2"], s=18, facecolor=SEX_COLORS[sex],
                    marker=VITRO_CONDITION_SHAPES[treatment], edgecolor="black",
                    linewidth=0.4, alpha=0.95, zorder=3,
                )
    ax_pca.axhline(0, color="#808080", ls="--", lw=0.4)
    ax_pca.axvline(0, color="#808080", ls="--", lw=0.4)
    ax_pca.set_xlabel(f"PC1 ({var.get('PC1',0):.0f}%)", fontsize=ns.FS_AXIS,
                      fontweight="bold")
    ax_pca.set_ylabel(f"PC2 ({var.get('PC2',0):.0f}%)", fontsize=ns.FS_AXIS,
                      fontweight="bold")
    ax_pca.set_title("In vitro ER / aromatase PCA", fontsize=ns.FS_TITLE,
                     fontweight="bold")
    ax_pca.tick_params(labelsize=ns.FS_TICK)
    sex_handles = [
        Patch(facecolor=SEX_COLORS[sex], edgecolor=SEX_COLORS[sex], alpha=0.7, label=sex)
        for sex in sex_order
    ]
    treatment_handles = [
        plt.Line2D([0], [0], marker=VITRO_CONDITION_SHAPES[treatment], linestyle="none",
                   markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.4,
                   markersize=4.2, label=treatment)
        for treatment in treatment_order
    ]
    legend_sex = fig.legend(
        handles=sex_handles, title="Biological group", loc="upper left",
        bbox_to_anchor=(0.845, 0.86), bbox_transform=fig.transFigure,
        borderaxespad=0.0,
    )
    fig.add_artist(legend_sex)
    fig.legend(
        handles=treatment_handles, title="Condition", loc="lower left",
        bbox_to_anchor=(0.845, 0.20), bbox_transform=fig.transFigure,
        borderaxespad=0.0,
    )
    ns.save_panel(fig, OUT / "invitro_estrogenic_score_and_pca", vector=True, tiff=True)
    print("[OK] in vitro directed estrogenic score + PCA")


_PCR = {
    "invitro_csf1_by_sex": build_invitro_csf1_by_sex_bubble, "invitro_methylation_regulators": build_invitro_methylation_regulators_bubble, "invivo_pain_markers": build_invivo_pain_markers_bubble,
    "invivo_methylation_regulators": build_invivo_methylation_regulators_bubble, "invivo_estrogen_genes": build_invivo_estrogen_genes_bubble, "invitro_estrogen_genes": build_invitro_estrogen_genes_bubble,
    "invitro_estrogenic_score": build_invitro_estrogenic_score_and_pca,
}

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="all")
    a = ap.parse_args()
    if a.only == "all":
        for fn in _PCR.values():
            fn()
    elif a.only in _PCR:
        _PCR[a.only]()
