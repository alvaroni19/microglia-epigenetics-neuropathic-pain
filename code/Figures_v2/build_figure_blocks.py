"""Draw every figure panel block at the exact size of its slot (figure_layout.py).

Data, statistics and display conventions are those of build_manuscript_panels_v2.py
and build_panels_v2_mfa_morphology.py (no test is recomputed); this module only
changes geometry and typography so that all blocks share one style and can be
placed at 100 % scale by assemble_figure_pptx.py.

Output: generated_outputs/Figures_v2/blocks/<block id>.{svg,pdf,png}
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path
from textwrap import wrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyArrowPatch, Patch, Rectangle
from matplotlib.textpath import TextPath

HERE = Path(__file__).absolute().parent
ROOT = HERE.parents[1]
CODE, GEN = ROOT / "code", ROOT / "generated_outputs"
OUT = GEN / "Figures_v2" / "blocks"

QPCR_VITRO = ROOT / "data" / "PCR_260514" / "invitro_qpcr_relative_expression.xlsx"
QPCR_VIVO = ROOT / "data" / "PCR_260514" / "invivo_qpcr_delta_ct.xlsx"
BEHAVIOUR = ROOT / "data" / "Allodynia_0_16" / "allodynia_change_scores.xlsx"
VON_FREY = ROOT / "data" / "Allodynia_0_16" / "allodynia_von_frey_scores.xlsx"
LINKAGE = ROOT / "data" / "MorphOMICs" / "morphology_behaviour_linkage.xlsx"
RESTRICTED = ("behaviour and RT-qPCR inputs are not distributed "
              "(available from the corresponding author on request)")


def _gate(label, *paths):
    """True when every input exists; otherwise report the panel as skipped."""
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        print(f"[SKIP] {label}: {missing[0].name} not found; {RESTRICTED}")
    return not missing

for sub in ("Figures_v2", "_helpers", "MeDIP", "PCR_260514", "MorphOMICs"):
    sys.path.insert(0, str(CODE / sub))

import figure_layout as L  # noqa: E402
import original_style as st  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


V2 = _load("build_manuscript_panels_v2", HERE / "build_manuscript_panels_v2.py")
MM2 = _load("build_panels_v2_mfa_morphology", HERE / "build_panels_v2_mfa_morphology.py")
GREEK, PHASES, TITLE = V2.GREEK, V2.PHASES, V2.TITLE
COND = MM2.COND
SLOT = dict(L.BLOCK_SIZES)


def stars(p):
    if p is None or not np.isfinite(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""


def text_mm(s, fs, italic=False, bold=False):
    """Width and height (mm) of a single-line string at font size fs."""
    fp = FontProperties(family="Arial", style="italic" if italic else "normal",
                        weight="bold" if bold else "normal")
    ext = TextPath((0, 0), s, size=fs, prop=fp).get_extents()
    return ext.width * 25.4 / 72, fs * 1.18 * 25.4 / 72


def block(bid):
    return L.Block(*SLOT[bid])


def legend_row(b, kind, y=0.3, x_center=None):
    """Shared one-line legend: 'Biological phase' (in vivo) or 'Sex' (in vitro)."""
    if kind == "phase":
        title = "Biological phase"
        hs = [Line2D([], [], marker="o", ls="", mfc=st.SAT[k], mec="black", mew=0.5, ms=L.LEG_MS + 1, label=k)
              for k in ["Male"] + PHASES]
    else:
        title = "Sex"
        hs = [Line2D([], [], marker="o", ls="", mfc="white", mec=st.SEX_EDGE[k], mew=1.3, ms=L.LEG_MS + 1,
                     label=k.capitalize()) for k in ("male", "female")]
    hs = [Line2D([], [], ls="", marker="", label=title)] + hs
    lg = b.fig.legend(handles=hs, ncol=len(hs), loc="upper center",
                      bbox_to_anchor=(b.fx(x_center if x_center is not None else b.w / 2), b.fy(y)),
                      fontsize=L.FS_LEGEND, handletextpad=0.15, columnspacing=0.7, handlelength=0.8)
    lg.get_texts()[0].set_fontsize(L.FS_LEGEND_TITLE)
    return lg


# ============================================================================ violin cells
def _vivo_stats():
    return pd.concat([pd.read_csv(GEN / "PCR_260514" / f) for f in (
        "invivo_pain_microglia_markers_stats.csv", "invivo_methylation_regulators_stats.csv",
        "invivo_estrogen_genes_stats.csv")], ignore_index=True)


def vivo_cell(b, x, y, gene, fc, stats, first_col=True):
    c = L.VIVO_CELL
    axm = b.ax(x + c["left"], y + c["top"], c["male_w"], c["ax_h"])
    axf = b.ax(x + c["left"] + c["male_w"] + c["gap"], y + c["top"], c["female_w"], c["ax_h"])
    alias = V2.VIVO_ALIAS.get(gene, gene)
    stx = stats[stats.gene == alias].set_index("comparison")["stars"]
    for ax, groups in ((axm, ["Male"]), (axf, PHASES)):
        his, ticks = [], []
        for i, ph in enumerate(groups):
            xc, xi = i * 2.3, i * 2.3 + 1
            sub = fc[fc.phase == ph]
            _, h1 = st.violin(ax, sub[sub.condition == "Control"][gene].dropna().values, xc, st.LIGHT[ph],
                              width=0.85, lw=0.5, point_size=L.DOT, median_w=0.5)
            _, h2 = st.violin(ax, sub[sub.condition == "CCI"][gene].dropna().values, xi, st.SAT[ph],
                              width=0.85, lw=0.5, point_size=L.DOT, median_w=0.5)
            his.append(np.nanmax([h1, h2]))
            ticks += [xc, xi]
        ymin, ymax = ax.get_ylim()
        span = ymax - ymin
        top = max(his)
        for i, ph in enumerate(groups):
            s = stx.get(f"{ph} CCI vs control")
            if isinstance(s, str) and s:
                L.bracket(ax, i * 2.3, i * 2.3 + 1, his[i] + 0.03 * span, 0.03 * span, s)
                top = max(top, his[i] + 0.20 * span)
        ax.set_ylim(ymin, max(ymax, top))
        ax.set_xlim(-0.7, (len(groups) - 1) * 2.3 + 1.7)
        ax.set_xticks(ticks)
        ax.set_xticklabels(["Control", "CCI-IoN"] * len(groups), rotation=45, ha="right", rotation_mode="anchor")
        ax.set_title(TITLE.get(gene, gene), fontstyle="italic", fontsize=L.FS_TITLE, pad=1.5)
        L.finish(ax)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    axm.set_ylabel("Fold change", labelpad=1)


def vivo_grid(bid, genes, ncols, legend=True):
    b = block(bid)
    fc = V2.invivo_fold_changes(pd.read_csv(GEN / "PCR_260514" / "invivo_qpcr_animal_level.csv"))
    stats = _vivo_stats()
    y0 = L.LEGEND_ROW if legend else 0.0
    if legend:
        legend_row(b, "phase")
    for k, g in enumerate(genes):
        r, cidx = divmod(k, ncols)
        vivo_cell(b, cidx * L.VIVO_CELL["w"], y0 + r * L.VIVO_CELL["h"], g, fc, stats)
    b.save(OUT / bid)


def vitro_axes_draw(ax, gene, data, brackets, order, labels, lw=1.0):
    pos = {g: i + 1.0 for i, g in enumerate(order)}
    his = {}
    for g in order:
        sex = "male" if g.startswith("male") or g.endswith("_male") else "female"
        treat = ("CSF1E2" if "E2" in g else "CSF1" if "csf1" in g.lower() else "control")
        _, his[g] = st.violin(ax, data[g], pos[g], st.TREAT_FILL[treat], edge=st.SEX_EDGE[sex], lw=lw,
                              width=0.82, point_size=L.DOT, median_w=0.45)
    ymin, ymax = ax.get_ylim()
    span = max(his.values()) - ymin
    tiers = []
    for g1, g2, s in sorted(brackets, key=lambda bb: abs(pos[bb[1]] - pos[bb[0]])):
        a, c = sorted((pos[g1], pos[g2]))
        level = max(his[g] for g in order if a <= pos[g] <= c) + 0.03 * span
        for (ta, tb, ty) in tiers:
            if not (c < ta or a > tb):
                level = max(level, ty + 0.13 * span)
        tiers.append((a, c, level))
        L.bracket(ax, a, c, level, 0.03 * span, s)
    t_top = max([t[2] for t in tiers], default=max(his.values()))
    # headroom for the asterisks above the highest bracket: 3.4 mm, whatever the axes height
    ax_h_mm = ax.get_position().height * ax.figure.get_figheight() * 25.4
    f = min(0.45, 3.4 / max(ax_h_mm, 1.0))
    top = max(t_top + 0.16 * span, (t_top - ymin * f) / (1 - f))
    ax.set_ylim(ymin, max(ymax, top))
    ax.set_xlim(0.4, len(order) + 0.6)
    ax.set_xticks(list(pos.values()))
    ax.set_xticklabels([labels[g] for g in order], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_title(TITLE.get(gene, gene), fontstyle="italic", fontsize=L.FS_TITLE, pad=1.5)
    ax.set_ylabel("Fold change", labelpad=1)
    L.finish(ax)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))


def _iv6():
    iv6 = importlib.import_module("build_invitro_6group_violins")
    lab6 = {g: iv6.GROUP_LABELS[g].replace("Ctrl", "C") for g in iv6.GROUPS}
    return iv6, lab6


def vitro_grid(bid, genes, ncols, legend=True):
    iv6, lab6 = _iv6()
    b = block(bid)
    c = L.VITRO_CELL
    y0 = L.LEGEND_ROW if legend else 0.0
    if legend:
        legend_row(b, "sex")
    for k, g in enumerate(genes):
        r, ci = divmod(k, ncols)
        ax = b.ax(ci * c["w"] + c["left"], y0 + r * c["h"] + c["top"], c["ax_w"], c["ax_h"])
        vitro_axes_draw(ax, g, iv6.load_data(g), iv6.load_brackets(g), iv6.GROUPS, lab6)
    b.save(OUT / bid)


def panel_invitro_qpcr_csf1r_violin():
    iv6, lab6 = _iv6()
    b = block("invitro_qpcr_csf1r_violin")
    ax = b.ax(14, 6, 58, 40)
    vitro_axes_draw(ax, "CSF1R", iv6.load_data("CSF1R"), iv6.load_brackets("CSF1R"), iv6.GROUPS, lab6, lw=1.2)
    legend_row(b, "sex", y=b.h - 4.2, x_center=14 + 29)
    b.save(OUT / "invitro_qpcr_csf1r_violin")


def panel_invitro_qpcr_methylation_regulators_csf1_violins():
    """CSF1 vs control in vitro panels, one column; groups ordered by sex (M-C, M-CSF1, F-C, F-CSF1)."""
    iv2 = importlib.import_module("build_invitro_2way_violins")
    b = block("invitro_qpcr_methylation_regulators_csf1_violins")
    order = ["control_male", "csf1_male", "control_female", "csf1_female"]
    lab = {"control_male": "M-C", "csf1_male": "M-CSF1", "control_female": "F-C", "csf1_female": "F-CSF1"}
    legend_row(b, "sex")
    cell_h = (b.h - L.LEGEND_ROW) / len(iv2.METHYLATION_GENES)
    for k, g in enumerate(iv2.METHYLATION_GENES):
        ax = b.ax(8.0, L.LEGEND_ROW + k * cell_h + 3.5, b.w - 9.0, cell_h - 3.5 - 8.5)
        br = [(x[0], x[1], x[2]) for x in iv2.load_brackets(g)]
        vitro_axes_draw(ax, g, iv2.load_data(g), br, order, lab)
    b.save(OUT / "invitro_qpcr_methylation_regulators_csf1_violins")


# ============================================================================ bubble matrices
def bubble(bid, df, comps, genes, col_labels, title_rows=None):
    b = block(bid)
    df = df[df["comparison"].isin(comps)].copy()
    genes = [g for g in genes if g in set(df["gene_label"])]
    nr, nc = len(genes), len(comps)
    cell = L.BUBBLE_CELL
    lab_w = max(text_mm(g, L.FS_LABEL, italic=True)[0] for g in genes) + 1.5
    xl = [col_labels.get(c, c) for c in comps]
    xl_w = max(text_mm(line, L.FS_TICK)[0] for s in xl for line in s.split("\n"))
    n_lines = max(s.count("\n") + 1 for s in xl)
    line_h = text_mm("Ag", L.FS_TICK)[1]
    rot_deg = 35 if n_lines == 1 else 55
    rot = np.deg2rad(rot_deg)
    bottom = xl_w * np.sin(rot) + n_lines * line_h * np.cos(rot) + 2.0
    left_need = max(lab_w, xl_w * np.cos(rot) + n_lines * line_h * np.sin(rot) - 0.5 * cell + 2.5)
    grid_w, grid_h = nc * cell, nr * cell
    cb_block = 3.0 + 2.6 + 9.0
    x0 = max(left_need, (b.w - grid_w - cb_block - left_need) / 2 + left_need)
    y0 = max(1.0, (b.h - grid_h - bottom) / 2)
    ax = b.ax(x0, y0, grid_w, grid_h)
    ymap = {g: nr - 1 - i for i, g in enumerate(genes)}
    xmap = {c: i for i, c in enumerate(comps)}
    for x in range(nc):
        for y in range(nr):
            ax.add_patch(Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor="#F7F7F7", edgecolor="#E2E2E2",
                                   lw=0.5, zorder=0))
    lim = float(np.nanmax(np.abs(df["log2fc"].to_numpy())))
    lim = (lim if np.isfinite(lim) and lim > 0 else 1.0) * 1.05
    norm = Normalize(-lim, lim)
    diam_pt = 0.84 * cell * L.MM * 72
    for _, r in df.iterrows():
        if r["gene_label"] not in ymap:
            continue
        x, y = xmap[r["comparison"]], ymap[r["gene_label"]]
        col = st.LOGFC_CMAP(norm(r["log2fc"])) if np.isfinite(r["log2fc"]) else (1, 1, 1, 1)
        ax.scatter(x, y, s=diam_pt ** 2, facecolor=col, edgecolor="#333333", lw=0.5, zorder=3)
        if isinstance(r["stars"], str) and r["stars"]:
            dark = abs(norm(r["log2fc"]) - 0.5) > 0.30
            L.star_text(ax, x, y, r["stars"], dark=dark, ha="center", va="center_baseline")
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(-0.5, nr - 0.5)
    ax.set_xticks(range(nc))
    ax.set_xticklabels(xl, rotation=rot_deg, ha="right", rotation_mode="anchor", fontsize=L.FS_TICK,
                       multialignment="center", linespacing=1.0)
    ax.set_yticks([ymap[g] for g in genes])
    ax.set_yticklabels(genes, fontsize=L.FS_LABEL, fontstyle="italic")
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb_h = min(grid_h, max(22.0, 0.8 * grid_h))
    cax = b.ax(x0 + grid_w + 3.0, y0 + (grid_h - cb_h) / 2, 2.6, cb_h)
    cb = b.fig.colorbar(plt.cm.ScalarMappable(cmap=st.LOGFC_CMAP, norm=norm), cax=cax)
    cb.solids.set_rasterized(False)
    cb.set_label("log$_2$ fold change", fontsize=L.FS_LABEL, fontweight="bold", labelpad=2)
    cb.ax.tick_params(labelsize=L.FS_TICK, width=0.4, length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)
    cb.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    b.save(OUT / bid)


def _bubble_df(fp, csv):
    df = pd.read_csv(GEN / "PCR_260514" / csv)
    df["gene_label"] = df["gene_label"].map(fp._disp_gene).map(lambda g: GREEK.get(g, g)).replace({"Tnfα": "Tnf-α"})
    return df


def bubbles():
    fp = _load("pcr_figure_panels", CODE / "PCR_260514" / "figure_panels.py")
    vitro = ["Male CSF1 vs control", "Male CSF1+E2 vs CSF1", "Female CSF1 vs control", "Female CSF1+E2 vs CSF1"]
    vlab = {"Male CSF1 vs control": "M-CSF1 vs M-C", "Male CSF1+E2 vs CSF1": "M-CSF1+E2 vs M-CSF1",
            "Female CSF1 vs control": "F-CSF1 vs F-C", "Female CSF1+E2 vs CSF1": "F-CSF1+E2 vs F-CSF1"}
    vivo = ["Male CCI vs control", "Proestrus CCI vs control", "Estrus CCI vs control", "Diestrus CCI vs control"]
    wlab = {"Male CCI vs control": "M-CCI-IoN\nvs M-C", "Proestrus CCI vs control": "Proestrus-CCI-IoN\nvs Proestrus-C",
            "Estrus CCI vs control": "Estrus-CCI-IoN\nvs Estrus-C", "Diestrus CCI vs control": "Diestrus-CCI-IoN\nvs Diestrus-C"}
    meth = [GREEK.get(g, g) for g in ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"]]
    er = [GREEK.get(fp._disp_gene(g), fp._disp_gene(g)) for g in fp.ER_ARO_GENES]
    d2 = _bubble_df(fp, "invitro_methylation_regulators_stats.csv")
    bubble("invitro_qpcr_methylation_regulators_csf1_bubble", d2, ["Male CSF1 vs control", "Female CSF1 vs control"], meth, vlab)
    bubble("invitro_qpcr_methylation_regulators_bubble", d2, vitro, meth, vlab)
    bubble("invivo_estrogen_genes_bubble", _bubble_df(fp, "invivo_estrogen_genes_stats.csv"), vivo, er, wlab)
    bubble("invitro_estrogen_genes_bubble", _bubble_df(fp, "invitro_estrogen_genes_stats.csv"), vitro, er, vlab)


# ============================================================================ MeDIP (DMR directionality, GO, candidate genes)
def panel_medip_dmr_directionality():
    mp = _load("medip_figure_panels", CODE / "MeDIP" / "figure_panels.py")
    stt = mp._directionality_stats().set_index("comparison")
    b = block("medip_dmr_directionality")
    specs = [("CSF1 effect within sex", [("F_CSF_vs_F_C", "F-CSF1 vs\nF-C"), ("M_CSF_vs_M_C", "M-CSF1 vs\nM-C")],
              lambda r: (r["prop_hypo"], r["prop_hyper"]),
              (st.CSF1_HYPO, "Hypomethylated after CSF1"), (st.CSF1_HYPER, "Hypermethylated after CSF1")),
             ("Sex effect (male vs female)", [("F_C_vs_M_C", "M-C vs\nF-C"), ("F_CSF_vs_M_CSF", "M-CSF1 vs\nF-CSF1")],
              lambda r: (r["prop_hyper"], r["prop_hypo"]),
              (st.SEX_FEMALE_HIGH, "Higher methylation in female"), (st.SEX_MALE_HIGH, "Higher methylation in male"))]
    half = b.w / 2
    for k, (title, rows, split, lsp, rsp) in enumerate(specs):
        ax = b.ax(k * half + 13.0, 4.0, half - 17.0, 13.0)
        for i, (key, lab) in enumerate(rows):
            lft, rgt = split(stt.loc[key])
            y = len(rows) - 1 - i
            ax.barh(y, -lft, color=lsp[0], height=0.66)
            ax.barh(y, rgt, color=rsp[0], height=0.66)
            for val, xc in ((lft, -lft / 2), (rgt, rgt / 2)):
                ax.text(xc, y, f"{100 * val:.0f}%", ha="center", va="center", color="white",
                        fontsize=L.FS_LABEL, fontweight="bold")
        ax.axvline(0, color="black", lw=L.LW_AXIS)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xticklabels(["100%", "50%", "0", "50%", "100%"])
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([lab for _, lab in rows][::-1], fontsize=L.FS_LABEL, fontweight="bold")
        ax.set_xlabel("Fraction of DMR regions")
        ax.set_title(title, fontweight="bold", pad=2)
        L.finish(ax)
        ax.tick_params(axis="y", length=0)
        b.fig.legend(handles=[Patch(fc=lsp[0], label=lsp[1]), Patch(fc=rsp[0], label=rsp[1])],
                     loc="upper left", bbox_to_anchor=(b.fx(k * half + 13.0), b.fy(25.5)), ncol=1,
                     handlelength=0.9, handleheight=0.7)
    b.save(OUT / "medip_dmr_directionality")


def _place_go_labels(b, ax, picks, pts_xy, colors):
    """Greedy placement: boxes go top-down; each takes the left-most position that
    overlaps neither the plotted points nor the boxes already placed."""
    b.fig.canvas.draw()
    to_mm = lambda xy: np.c_[xy[:, 0] / b.fig.dpi * 25.4, xy[:, 1] / b.fig.dpi * 25.4]
    disp = ax.transData.transform(pts_xy)
    P = to_mm(disp)
    bb = ax.get_window_extent()
    ax_x0, ax_y0 = bb.x0 / b.fig.dpi * 25.4, bb.y0 / b.fig.dpi * 25.4
    ax_x1, ax_y1 = bb.x1 / b.fig.dpi * 25.4, bb.y1 / b.fig.dpi * 25.4
    pad, gap = 0.35, 0.6
    placed = []
    y_cursor = ax_y1 - 0.6
    for row, txt, col in picks:
        lines = txt.split("\n")
        tw = max(text_mm(s, L.FS_LEGEND)[0] for s in lines) + 2 * pad
        th = len(lines) * L.FS_LEGEND * 1.15 * 25.4 / 72 + 2 * pad
        best = None
        y_top = y_cursor
        while best is None and y_top - th > ax_y0 + 0.5:
            for x_left in np.arange(ax_x0 + 0.12 * (ax_x1 - ax_x0), ax_x1 - tw - 0.3, 0.5):
                x_right, y_bot = x_left + tw, y_top - th
                hit = ((P[:, 0] > x_left - gap) & (P[:, 0] < x_right + gap) &
                       (P[:, 1] > y_bot - gap) & (P[:, 1] < y_top + gap)).any()
                over = any(not (x_right + gap < a or x_left - gap > c or y_top + gap < d or y_bot - gap > e)
                           for a, d, c, e in placed)
                if not hit and not over:
                    best = (x_left, y_bot, x_right, y_top)
                    break
            y_top -= 0.5
        if best is None:
            best = (ax_x1 - tw - 0.3, ax_y0 + 1, ax_x1 - 0.3, ax_y0 + 1 + th)
        placed.append(best)
        y_cursor = best[1] - gap
        # position in data coordinates (left, vertical centre)
        inv = ax.transData.inverted()
        lx, ly = inv.transform(((best[0] + pad) / 25.4 * b.fig.dpi, (best[1] + best[3]) / 2 / 25.4 * b.fig.dpi))
        ax.annotate(txt, xy=(row["rank_within_group"], row["neglog10_pvalue"]), xytext=(lx, ly),
                    fontsize=L.FS_LEGEND, color=col, ha="left", va="center", linespacing=1.0,
                    bbox=dict(boxstyle=f"square,pad={pad / (L.FS_LEGEND * 25.4 / 72):.2f}", fc="white", ec=col, lw=0.5),
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=0.35, ls=(0, (2, 1.5)),
                                    shrinkA=0, shrinkB=1.0), zorder=6)


def panel_medip_go_rank_distribution():
    from collections import Counter
    mp = _load("medip_figure_panels", CODE / "MeDIP" / "figure_panels.py")
    df = pd.read_excel(V2.DATA_MEDIP / "GO_relevance_rank_distribution_ws400.xlsx")
    df["biological_theme"] = df["biological_theme"].replace({"Steroid / sex hormones": "Sex hormones / estrogen"})
    df = df[np.isfinite(df["rank_within_group"]) & np.isfinite(df["neglog10_pvalue"])].copy()
    tokens = df["Description"].map(mp._label_tokens)
    doc_freq = Counter(t for ts in tokens for t in ts)
    panels = [("M", "GO M CSF1 vs M C", "M-CSF1 vs M-C"), ("F", "GO F CSF1 vs F C", "F-CSF1 vs F-C")]
    go_key = {"M": "M_CSF_vs_M_C", "F": "F_CSF_vs_F_C"}
    b = block("medip_go_rank_distribution")
    left, top, colgap, rowgap, legend_h = 11.0, 4.0, 6.0, 8.5, 6.0
    aw = (b.w - left - colgap - 1.0) / 2
    ah = (b.h - top - legend_h - rowgap - 7.0) / 2
    used = []
    for r, disp_dir in enumerate(["Hyper", "Hypo"]):
        for c, (sex, rel_comp, title) in enumerate(panels):
            rel_dir = st.medip_direction(go_key[sex], disp_dir)
            ax = b.ax(left + c * (aw + colgap), top + r * (ah + rowgap), aw, ah)
            sub = df[(df["comparison_label"] == rel_comp) & (df["direction"] == rel_dir)].copy()
            ax.scatter(sub["rank_within_group"], sub["neglog10_pvalue"], s=1.0, color="#C9C9C9", lw=0, zorder=2)
            for theme, col in mp.GO_THEME_COLORS.items():
                th = sub[sub["biological_theme"] == theme]
                ax.scatter(th["rank_within_group"], th["neglog10_pvalue"], s=3.0, color=col, lw=0, zorder=3)
            picks = []
            for theme, term in V2.OLD_GO_TERMS[(sex, disp_dir)].items():
                hit = sub[sub["Description"].str.lower() == term.lower()]
                if len(hit):
                    picks.append((hit.sort_values("neglog10_pvalue", ascending=False).iloc[0], "preferred term"))
                else:
                    rep = mp._representative_theme_term(sub[sub["biological_theme"] == theme], theme, doc_freq, len(df))
                    if len(rep):
                        picks.append((rep.iloc[0], f"representative term ('{term}' not enriched)"))
            picks = sorted(picks, key=lambda p: -p[0]["neglog10_pvalue"])
            x_max = max(1.0, sub["rank_within_group"].max())
            y_top = max(1.0, sub["neglog10_pvalue"].max())
            ax.set_xlim(-0.02 * x_max, x_max * 1.03)
            ax.set_ylim(0, y_top * 1.08)
            if r == 0:
                ax.set_title(title, fontweight="bold", pad=2)
            if c == 0:
                ax.set_ylabel(f"{disp_dir}\n-log$_{{10}}$(P value)", labelpad=1)
            if r == 1:
                ax.set_xlabel("Rank", labelpad=1)
            ax.grid(True, color="#EEEEEE", lw=0.3, zorder=0)
            L.finish(ax)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
            lab = [(row, "\n".join(wrap(str(row["Description"]), 32)),
                    mp.GO_THEME_COLORS.get(row["biological_theme"], "#666666")) for row, _ in picks]
            _place_go_labels(b, ax, lab, sub[["rank_within_group", "neglog10_pvalue"]].to_numpy(float), None)
            for row, src in picks:
                used.append(dict(panel=f"{title} {disp_dir}", release_contrast=rel_comp, release_direction=rel_dir,
                                 term=row["Description"], theme=row["biological_theme"], region=row["region"],
                                 rank=int(row["rank_within_group"]), neglog10_p=round(float(row["neglog10_pvalue"]), 3),
                                 source=src))
    handles = [Line2D([0], [0], marker="o", ls="", mfc=col, mec="none", ms=L.LEG_MS, label=mp.GO_THEME_LEGEND[t])
               for t, col in mp.GO_THEME_COLORS.items()]
    hs = [Line2D([], [], ls="", marker="", label="Biological theme")] + handles
    lg = b.fig.legend(handles=hs, loc="lower center", ncol=len(hs), bbox_to_anchor=(0.5, 0.0),
                      handletextpad=0.1, columnspacing=0.45, borderaxespad=0.1)
    lg.get_texts()[0].set_fontsize(L.FS_LEGEND_TITLE)
    b.save(OUT / "medip_go_rank_distribution")
    pd.DataFrame(used).to_csv(OUT / "medip_go_rank_distribution_labels_used.csv", index=False)


def panel_medip_candidate_gene_heatmap():
    hits = pd.read_excel(V2.DATA_MEDIP / "methylation_map_candidate_genes_hits_ws400.xlsx")
    keep = ["F_CSF_vs_F_C", "M_CSF_vs_M_C"]
    hits = hits[hits["comparison"].isin(keep)].copy()
    hits["v"] = hits["edgeR.logFC"] * hits["comparison"].map(lambda k: st.MEDIP_CONTRAST[k][1])
    gc = hits.groupby(["symbol", "comparison"]).agg(v=("v", "mean"), p=("edgeR.p.value", "min")).reset_index()
    order = [g for g in V2.CANDIDATE_GENE_ORDER if g in set(gc["symbol"])]
    z = gc.pivot(index="symbol", columns="comparison", values="v").reindex(index=order, columns=keep)
    pv = gc.pivot(index="symbol", columns="comparison", values="p").reindex(index=order, columns=keep)
    M = z.to_numpy(float)
    vmax = float(np.nanmax(np.abs(M)))
    nr, nc = M.shape
    cw, ch = L.HEAT_CELL
    b = block("medip_candidate_gene_heatmap")
    labs = [GREEK.get(g, g) for g in order]
    left = max(text_mm(s, L.FS_LABEL, italic=True)[0] for s in labs) + 2.0
    xl = [st.medip_label(k) for k in keep]
    bottom = max(text_mm(s, L.FS_TICK)[0] for s in xl) * np.sin(np.deg2rad(35)) + 3
    grid_w, grid_h = nc * cw, nr * ch
    x0 = max(left, (b.w - grid_w - 14.6) / 2)
    y0 = max(1.0, (b.h - grid_h - bottom) / 2)
    ax = b.ax(x0, y0, grid_w, grid_h)
    ax.add_patch(Rectangle((-0.5, -0.5), nc, nr, facecolor="#F7F7F7", edgecolor="none", zorder=0))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    ax.pcolormesh(np.arange(nc + 1) - 0.5, np.arange(nr + 1) - 0.5, np.ma.masked_invalid(M), cmap=st.LOGFC_CMAP,
                  norm=norm, shading="flat", edgecolors="none", zorder=1)
    for i in range(nr):
        for j in range(nc):
            s = stars(float(pv.iloc[i, j])) if np.isfinite(pv.iloc[i, j]) else ""
            if s:
                L.star_text(ax, j, i, s, dark=abs(M[i, j]) > 0.55 * vmax, ha="center", va="center_baseline")
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.set_yticks(range(nr))
    ax.set_yticklabels(labs, fontsize=L.FS_LABEL, fontstyle="italic")
    ax.set_xticks(range(nc))
    ax.set_xticklabels(xl, rotation=35, ha="right", rotation_mode="anchor", fontsize=L.FS_TICK,
                       multialignment="center", linespacing=1.0)
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = b.ax(x0 + grid_w + 3.0, y0 + 0.1 * grid_h, 2.6, 0.8 * grid_h)
    cb = b.fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=st.LOGFC_CMAP), cax=cax)
    cb.solids.set_rasterized(False)
    cb.set_label("log$_2$ fold change", fontsize=L.FS_LABEL, fontweight="bold", labelpad=2)
    cb.ax.tick_params(labelsize=L.FS_TICK, width=0.4, length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)
    b.save(OUT / "medip_candidate_gene_heatmap")


# ============================================================================ Csf1r locus map
def panel_medip_csf1r_locus_dmr_ere_map():
    import logomaker
    import map_csf1r_esr_motifs_integrated as mc
    gene_start, gene_stop, strand, exons = mc.read_gene_model()
    rec = mc.read_fasta(V2.DATA_MEDIP / "csf1r_context_rn6.fasta")[0]
    x0, x1 = int(rec.start), int(rec.start) + len(rec.sequence) - 1
    motifs = [mc.read_meme(ROOT / "data" / "MeDIP" / "motifs" / f) for f in ("MA0112.3.meme", "MA0258.2.meme")]
    pwm_hits = pd.concat([mc.scan_record(rec, m, 0.92) for m in motifs], ignore_index=True)
    win = pd.read_csv(GEN / "Figures_v2" / "panels" / "medip_csf1r_locus_windows_display.csv")
    win = win[win["plotted_colour_from_this_contrast"]].drop_duplicates("window_id").sort_values("start")
    b = block("medip_csf1r_locus_dmr_ere_map")
    L0, W = 12.0, b.w - 14.0
    gene_ax = b.ax(L0, 0.5, W, 7.5)
    cpg_ax = b.ax(L0, 9.0, W, 7.0)
    trk_ax = b.ax(L0, 16.5, W, 7.5)
    for a in (gene_ax, cpg_ax, trk_ax):
        a.set_xlim(x0, x1)
    gene_ax.add_patch(Rectangle((gene_start, -0.08), gene_stop - gene_start, 0.16, facecolor="#2CA02C", edgecolor="none"))
    for _, e in exons.iterrows():
        gene_ax.add_patch(Rectangle((e["start"], -0.4), max(int(e["stop"] - e["start"]), 120), 0.8,
                                    facecolor="#9B59B6", edgecolor="#6C3483", lw=0.3, zorder=3))
    tss = gene_start if strand == 1 else gene_stop
    gene_ax.add_patch(Rectangle((tss - 420, -0.48), 380, 0.96, facecolor="#F39C12", edgecolor="none", zorder=4))
    for xx, s, col in ((tss - 230, "TSS", "#F39C12"), (gene_start + 0.06 * (gene_stop - gene_start), "GENE BODY", "#2CA02C"),
                       (gene_start + 0.26 * (gene_stop - gene_start), "EXONS", "#9B59B6")):
        gene_ax.text(xx, -0.62, s, ha="center" if s == "TSS" else "left", va="top", fontsize=L.FS_TITLE,
                     color=col, fontweight="bold")
    gene_ax.text(x1, -0.62, "Csf1r (rn6)", ha="right", va="top", fontsize=L.FS_TITLE, fontstyle="italic")
    gene_ax.set_ylim(-1.6, 1.0)
    gene_ax.axis("off")
    seq = rec.sequence.upper()
    step = 400
    starts = np.arange(0, len(seq), step)
    dens = [1000 * sum(seq[i + k:i + k + 2] == "CG" for k in range(min(step, len(seq) - i) - 1)) / max(1, min(step, len(seq) - i))
            for i in starts]
    cpg_ax.bar(x0 + starts, dens, width=step * 0.92, align="edge", color="#C8C8C8", edgecolor="none")
    top = int(np.ceil(max(dens) / 20.0) * 20)
    cpg_ax.set_ylim(0, top)
    cpg_ax.set_yticks([0, top // 2, top])
    cpg_ax.text(x0 + 0.5 * (x1 - x0), top, "CpG / kb", ha="center", va="top", fontsize=L.FS_TICK)
    cpg_ax.spines[["top", "right", "bottom"]].set_visible(False)
    cpg_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    cpg_ax.tick_params(axis="y", labelsize=L.FS_TICK, length=L.TICK_LEN, width=L.LW_AXIS, pad=1)
    trk_ax.axhline(0.62, color="#BBBBBB", lw=1.0, zorder=1)
    for _, h in pwm_hits.iterrows():
        inside = ((win["start"] <= h["stop"]) & (win["stop"] >= h["start"])).any()
        trk_ax.add_patch(Rectangle((h["start"], 0.45), max(int(h["stop"] - h["start"]), 110), 0.34,
                                   facecolor="black" if inside else "white", edgecolor="black", lw=0.4,
                                   zorder=5 if inside else 3))
    for _, w in win.iterrows():
        col = st.HYPER if str(w["direction_display"]).upper().startswith("HYPER") else st.HYPO
        trk_ax.add_patch(Rectangle((w["start"], 0.36), w["stop"] - w["start"], 0.52, facecolor="white", edgecolor=col,
                                   lw=0.9, zorder=4))
        nb = win[(win["start"] - w["start"]).abs().between(1, 1300)]
        dx = 0 if nb.empty else (-450 if (nb["start"] > w["start"]).any() else 450)
        trk_ax.text((w["start"] + w["stop"]) / 2 + dx, 0.22, w["window_id"], ha="center", va="top",
                    fontsize=L.FS_TITLE, color=col, fontweight="bold")
    trk_ax.set_ylim(-0.75, 1.0)
    trk_ax.set_yticks([])
    trk_ax.spines[["top", "right", "left"]].set_visible(False)
    trk_ax.ticklabel_format(style="plain", axis="x", useOffset=False)
    trk_ax.tick_params(axis="x", labelsize=L.FS_TICK, length=L.TICK_LEN, width=L.LW_AXIS, pad=1)
    cx = L0 + W / 2
    b.text(cx - 1, 31.5, "MeDIP 400 bp windows", ha="right", va="center", fontsize=L.FS_LEGEND_TITLE, fontweight="bold")
    b.fig.legend(handles=[Patch(fc="white", ec=st.HYPO, lw=0.9, label="Hypomethylated DMR"),
                          Patch(fc="white", ec=st.HYPER, lw=0.9, label="Hypermethylated DMR")],
                 loc="center left", bbox_to_anchor=(b.fx(cx), b.fy(31.5)), ncol=2, handlelength=0.9,
                 handletextpad=0.3, columnspacing=0.8)
    b.text(cx - 1, 35.5, "ERE", ha="right", va="center", fontsize=L.FS_LEGEND_TITLE, fontweight="bold")
    b.text(cx, 35.5, "PWM relative score ≥ 0.92 (filled: inside a DMR window)", ha="left", va="center",
           fontsize=L.FS_LEGEND)
    lw_mm, lh = 46.0, 25.0
    for i, (m, name) in enumerate(zip(motifs, ["ERα", "ERβ"])):
        lx = cx - 8 - lw_mm if i == 0 else cx + 14
        lax = b.ax(lx, 44.0, lw_mm, lh)
        prob = pd.DataFrame(m.pwm, columns=list("ACGT"))
        info = logomaker.transform_matrix(prob, from_type="probability", to_type="information")
        logomaker.Logo(info, ax=lax, color_scheme="classic", font_name="Arial", vpad=0.02)
        lax.set_ylim(0, 2)
        lax.set_ylabel("Bits", labelpad=1)
        lax.set_xticks([])
        L.finish(lax)
        lax.spines["bottom"].set_visible(False)
        lax.text(0.5, 1.03, name, transform=lax.transAxes, fontsize=L.FS_TITLE, fontweight="bold", ha="center", va="bottom")
        b.fig.add_artist(FancyArrowPatch((b.fx(cx), b.fy(37.5)), (b.fx(lx + lw_mm / 2), b.fy(40.2)),
                                         transform=b.fig.transFigure, arrowstyle="-|>", mutation_scale=6, lw=0.6,
                                         color="black"))
    b.save(OUT / "medip_csf1r_locus_dmr_ere_map")


# ============================================================================ Behaviour
def _behaviour():
    B = GEN / "Behaviour_reanalysis"
    d = pd.read_csv(B / "plot_data.csv")
    emm = pd.read_csv(B / "emm_HC3.csv")
    sx = pd.read_excel(B / "behaviour_reanalysis_stats.xlsx", sheet_name=None)
    return d, emm, sx


def panel_behaviour_change_score_violins():
    d, emm, sx = _behaviour()
    inj = sx["injury_effect_by_group"].set_index("Group")
    cci = sx["CCI_phase_comparisons"].set_index("contrast")
    pint = sx["ANOVA_females_HC3"].set_index("term").loc["Condition:Group", "Pr(>F)"]
    star = lambda p: stars(p) or "n.s."
    b = block("behaviour_change_score_violins")
    legend_row(b, "phase", x_center=14 + (b.w - 15) / 2)
    ax = b.ax(14.0, 7.0, b.w - 15.0, b.h - 7.0 - 9.5)
    G = ["Male"] + PHASES
    rng = np.random.default_rng(1)
    pos, ext = {}, {}
    for i, g in enumerate(G):
        for j, (cond, col) in enumerate([("Control", st.LIGHT[g]), ("CCI", st.SAT[g])]):
            x = i * 2.6 + j
            pos[(g, cond)] = x
            y = d.query("Group==@g and Condition==@cond")["Delta"].values
            ext[(g, cond)] = st.violin(ax, y, x, col, edge="black", lw=0.5, width=0.8, point_size=0, median_w=0)
            ax.scatter(x + rng.uniform(-.12, .12, len(y)), y, s=L.DOT, c="black", lw=0, zorder=3)
            e = emm.query("Group==@g and Condition==@cond").iloc[0]
            ax.plot([x - .22, x + .22], [e.emmean] * 2, c="black", lw=1.2, zorder=4)
            ax.plot([x, x], [e["lower.CL"], e["upper.CL"]], c="black", lw=0.7, zorder=4)
    for g in G:
        top = max(ext[(g, "Control")][1], ext[(g, "CCI")][1])
        L.bracket(ax, pos[(g, "Control")], pos[(g, "CCI")], top + .18, .12, star(inj.loc[g, "p_holm"]))
    ytop = max(v[1] for v in ext.values()) + .75
    L.bracket(ax, pos[("Proestrus", "CCI")], pos[("Diestrus", "CCI")], ytop, .12,
              star(cci.loc["Proestrus - Diestrus", "p_holm_females"]))
    ax.axhline(0, ls=":", c="grey", lw=0.6, zorder=0)
    ax.set_xticks([pos[k] for k in pos])
    ax.set_xticklabels(["Control", "CCI-IoN"] * 4, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Δ mechanical response (0.16 g)\n(test − baseline)")
    ax.set_ylim(min(v[0] for v in ext.values()) - .15, ytop + 0.9)
    ax.set_xlim(-.7, pos[("Diestrus", "CCI")] + .7)
    ax.text(.01, .995, f"Condition × phase (females): P = {pint:.3f}", transform=ax.transAxes, ha="left", va="top",
            fontsize=L.FS_LEGEND)
    L.finish(ax)
    b.save(OUT / "behaviour_change_score_violins")


# ============================================================================ In vivo correlations, MFA and estrogenic score
def panel_invivo_marker_pain_spearman():
    df = pd.read_csv(GEN / "PCR_MFA" / "invivo_marker_pain_spearman.csv").sort_values("rho")
    labs = df["gene_label"].map(lambda g: GREEK.get(g, g)).tolist()
    b = block("invivo_marker_pain_spearman")
    left = max(text_mm(s, L.FS_LABEL, italic=True)[0] for s in labs) + 2.5
    ax = b.ax(left, 2.0, b.w - left - 2.5, b.h - 2.0 - 9.0)
    y = np.arange(len(df))
    ax.hlines(y, 0, df["rho"], color="#BBBBBB", lw=0.6, zorder=1)
    for yi, (_, r) in zip(y, df.iterrows()):
        col = "#D7301F" if r["rho"] > 0 else "#2C7FB8"
        sig = r["q_BH"] < 0.05
        ax.scatter(r["rho"], yi, s=60 if sig else L.SCATTER, color=col, edgecolor="#555555" if sig else "none",
                   lw=0.6, zorder=3)
        if sig:
            L.star_text(ax, r["rho"], yi, "*", dark=True, ha="center", va="center_baseline")
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=L.FS_LABEL, fontstyle="italic")
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xticklabels(["−1", "−0.5", "0", "0.5", "1"])
    ax.set_xlabel("Spearman correlation ρ")
    ax.set_ylim(-0.6, len(df) - 0.4)
    L.finish(ax)
    ax.tick_params(axis="y", length=0)
    b.save(OUT / "invivo_marker_pain_spearman")


def panel_estrogenic_score_vs_pain():
    d5 = pd.read_csv(GEN / "PCR_MFA" / "invivo_estrogenic_score_delta_pain_per_animal.csv")
    b = block("estrogenic_score_vs_pain")
    ax = b.ax(10.0, 2.0, b.w - 10.0 - 21.0, b.h - 2.0 - 8.0)
    x, y = d5["delta_estrogenic_score"].to_numpy(float), d5["delta_pain_016"].to_numpy(float)
    bb, a = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, a + bb * xs, color="black", lw=0.8, zorder=2)
    for ph in PHASES:
        s = d5[d5["phase"] == ph]
        ax.scatter(s["delta_estrogenic_score"], s["delta_pain_016"], s=L.SCATTER, color=st.SAT[ph],
                   edgecolor="black", lw=0.4, zorder=3)
    ax.set_xlabel("Δ estrogenic expression score")
    ax.set_ylabel("Δ mechanical response (0.16 g)")
    L.finish(ax)
    lg = ax.legend(handles=[Line2D([0], [0], marker="o", ls="", mfc=st.SAT[p], mec="black", mew=0.4, ms=L.LEG_MS,
                                   label=p) for p in PHASES], title="Biological phase", loc="upper left",
                   bbox_to_anchor=(1.02, 1.0))
    lg._legend_box.align = "left"
    b.save(OUT / "estrogenic_score_vs_pain")


def _two_legends(ax, h1, t1, h2, t2, y2=0.48):
    l1 = ax.legend(handles=h1, title=t1, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    l1._legend_box.align = "left"
    ax.add_artist(l1)
    l2 = ax.legend(handles=h2, title=t2, loc="upper left", bbox_to_anchor=(1.02, y2))
    l2._legend_box.align = "left"


def panel_invivo_mfa_phase_condition():
    df = pd.read_csv(GEN / "PCR_MFA" / "invivo_mfa_coords.csv")
    eig = pd.read_csv(GEN / "PCR_MFA" / "invivo_mfa_eig.csv")
    b = block("invivo_mfa_phase_condition")
    ax = b.ax(10.0, 2.0, b.w - 10.0 - 21.0, b.h - 2.0 - 8.0)
    ax.axvline(0, color="#BBBBBB", lw=0.5, zorder=1)
    ax.axhline(0, color="#BBBBBB", lw=0.5, zorder=1)
    for i, ph in enumerate(["Male", "Diestrus", "Estrus", "Proestrus"]):
        sub = df[df["phase"] == ph]
        MM2._cov_ellipse(ax, sub["Dim1"].to_numpy(), sub["Dim2"].to_numpy(), MM2.PHASE_FILL[ph], MM2.PHASE_EDGE[ph],
                         2 + 2 * i, level=0.68)
        for cond, mk in (("Control", "o"), ("CCI", "^")):
            sub_c = sub[sub["condition"] == cond]
            ax.scatter(sub_c["Dim1"], sub_c["Dim2"], marker=mk, s=L.SCATTER, facecolor=MM2.PHASE_FILL[ph],
                       edgecolor="#333333", lw=0.4, zorder=3 + 2 * i)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=MM2.PHASE_FILL[p], mec="#333333", mew=0.4, ms=L.LEG_MS, label=p)
          for p in ["Male", "Proestrus", "Estrus", "Diestrus"]]
    h2 = [Line2D([0], [0], marker=m, ls="", mfc="#999999", mec="#333333", mew=0.4, ms=L.LEG_MS, label=l)
          for m, l in (("o", "Control"), ("^", "CCI-IoN"))]
    _two_legends(ax, h1, "Biological phase", h2, "Condition", 0.42)
    ax.set_xlabel(f"Dimension 1 ({float(eig['var1'][0]):.1f}%)")
    ax.set_ylabel(f"Dimension 2 ({float(eig['var2'][0]):.1f}%)")
    L.finish(ax)
    b.save(OUT / "invivo_mfa_phase_condition")


def panel_invitro_estrogenic_score_and_pca():
    scr = GEN / "PCR_260514"
    score = pd.read_csv(scr / "invitro_estrogenic_score_per_replicate.csv")
    sc = pd.read_csv(scr / "invitro_estrogen_genes_pca_scores.csv")
    var = pd.read_csv(scr / "invitro_estrogen_genes_pca_variance.csv").set_index("PC")["variance_explained_pct"]
    SF = MM2.SEX_FILL
    contrasts = [("Male", "CSF1", "Control", "M-CSF1\nvs M-C"), ("Male", "CSF1+E2", "CSF1", "M-CSF1+E2\nvs M-CSF1"),
                 ("Female", "CSF1", "Control", "F-CSF1\nvs F-C"), ("Female", "CSF1+E2", "CSF1", "F-CSF1+E2\nvs F-CSF1")]
    b = block("invitro_estrogenic_score")
    ax = b.ax(10.0, 7.0, b.w - 12.0, b.h - 7.0 - 8.0)
    rng = np.random.default_rng(20260608)
    for x, (sex, num, den, lab) in enumerate(contrasts):
        a = score.loc[(score.sex == sex) & (score.treatment == num), "ESI_full"].dropna().to_numpy(float)
        d = a - np.median(score.loc[(score.sex == sex) & (score.treatment == den), "ESI_full"].dropna().to_numpy(float))
        ax.bar(x, np.median(d), width=0.62, color=SF[sex], alpha=0.35, edgecolor=SF[sex], lw=0.5, zorder=2)
        ax.scatter(x + rng.normal(0, 0.05, len(d)), d, s=L.SCATTER, facecolor=SF[sex], edgecolor="#333333", lw=0.4, zorder=3)
    ax.axhline(0, color="#555555", lw=0.5, ls=(0, (1, 1)))
    ax.axvline(1.5, color="#BBBBBB", lw=0.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([c[3] for c in contrasts])
    ax.set_ylabel("Δ estrogenic score")
    ax.grid(axis="y", color="#E6E6E6", lw=0.4)
    ax.set_axisbelow(True)
    L.finish(ax)
    hs = [Line2D([], [], ls="", marker="", label="Biological group")] + \
         [Patch(fc=SF[s], alpha=0.5, ec=SF[s], label=s) for s in ("Male", "Female")]
    lg = b.fig.legend(handles=hs, loc="upper center", ncol=3, bbox_to_anchor=(b.fx(10 + (b.w - 12) / 2), b.fy(0.3)))
    lg.get_texts()[0].set_fontsize(L.FS_LEGEND_TITLE)
    b.save(OUT / "invitro_estrogenic_score")

    b = block("invitro_estrogen_genes_pca")
    ax = b.ax(10.0, 2.0, b.w - 10.0 - 18.0, b.h - 2.0 - 8.0)
    for sex in ("Male", "Female"):
        s = sc[sc.sex == sex]
        MM2._cov_ellipse(ax, s["PC1"].to_numpy(), s["PC2"].to_numpy(), SF[sex], SF[sex], 1)
    shapes = {"Control": "o", "CSF1": "s", "CSF1+E2": "^"}
    for tr, mk in shapes.items():
        for sex in ("Male", "Female"):
            g = sc[(sc.sex == sex) & (sc.treatment == tr)]
            ax.scatter(g["PC1"], g["PC2"], marker=mk, s=L.SCATTER, facecolor=SF[sex], edgecolor="#333333", lw=0.4, zorder=3)
    ax.axhline(0, color="#777777", lw=0.4, ls=(0, (1, 1)))
    ax.axvline(0, color="#777777", lw=0.4, ls=(0, (1, 1)))
    ax.set_xlabel(f"Dimension 1 ({var['PC1']:.0f}%)")
    ax.set_ylabel(f"Dimension 2 ({var['PC2']:.0f}%)")
    L.finish(ax)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=SF[s], mec="#333333", mew=0.4, ms=L.LEG_MS, label=s) for s in ("Male", "Female")]
    h2 = [Line2D([0], [0], marker=m, ls="", mfc="white", mec="black", mew=0.5, ms=L.LEG_MS, label=t) for t, m in shapes.items()]
    _two_legends(ax, h1, "Sex", h2, "Condition", 0.55)
    b.save(OUT / "invitro_estrogen_genes_pca")


# ============================================================================ Morphology
def panel_morphology_ramification_index():
    pc = GEN / "MorphOMICs" / "cam_sholl_per_cell.csv"
    if not pc.exists():
        pc = ROOT / "data" / "MorphOMICs" / "reports_cam" / "cam_sholl_per_cell.csv"
    d = pd.read_csv(pc)
    d = d[(d["Layer"] == 12) & (d["Side"] == "IPSI")].dropna(subset=["ShollRamificationIndex"])
    lmm = pd.read_csv(GEN / "MorphOMICs" / "cam_ramification_index_layer12_animal_and_mixedmodel.csv").set_index("cycle")
    b = block("morphology_ramification_index")
    left, gap = 9.0, 1.5
    aw = (b.w - left - 3 * gap - 0.5) / 4
    rng = np.random.default_rng(7)
    ri = "ShollRamificationIndex"
    axes = []
    for k, cy in enumerate(["Male", "Proestrus", "Estrus", "Diestrus"]):
        ax = b.ax(left + k * (aw + gap), 4.0, aw, b.h - 4.0 - 4.5, sharey=axes[0] if axes else None)
        axes.append(ax)
        sub = d[d["Phase"].replace({"Metestrus": "Diestrus"}) == cy]
        tops = []
        for x, g in ((1, "Control"), (2, "CCI")):
            v = sub.loc[sub["Condition"] == g, ri].to_numpy(float)
            _, hi = st.violin(ax, v, x, COND[g], width=0.8, lw=0.5, point_size=0, median_w=0)
            ax.scatter(x + rng.normal(0, 0.06, len(v)), v, s=2.0, color="#333333", alpha=0.35, lw=0, zorder=3)
            ax.scatter([x], [v.mean()], s=30, color="black", edgecolor="white", lw=0.9, zorder=5)
            tops.append(hi)
        s = stars(float(lmm.loc[MM2.CYCLE_EN[cy], "p_value"]))
        if s:
            L.bracket(ax, 1, 2, max(tops) + 0.15, 0.25, s)
        ax.set_title(MM2.CYCLE_EN[cy], fontweight="bold", pad=1.5)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Control", "CCI-IoN"])
        ax.set_xlim(0.4, 2.6)
        L.finish(ax)
        ax.tick_params(axis="x", length=0)
        if k:
            ax.tick_params(axis="y", left=False, labelleft=False)
    axes[0].set_ylabel("Ramification index")
    axes[0].set_ylim(min(0.8, d[ri].min() - 0.2), axes[0].get_ylim()[1] + 0.3)
    axes[0].set_yticks(range(1, int(axes[0].get_ylim()[1]) + 1, 2))
    b.save(OUT / "morphology_ramification_index")


def panel_morphology_umap():
    u = pd.read_csv(ROOT / "data" / "MorphOMICs" / "saved_instances_cam" / "Morphomics.PID_CAM_CCI_260224.cam_reduced.csv",
                    index_col=0)
    u["Phase"] = u["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
    u = u[(u["Layer"].astype(str) == "12") & (u["Side"] == "IPSI")].copy()
    xy = u[["umap___1", "umap___2"]]
    med = xy.median()
    mad = (xy - med).abs().median().replace(0, 1)
    u["in_zoom"] = np.sqrt((((xy - med) / mad) ** 2).sum(axis=1)) <= 6.0

    def lim(v, f=0.04):
        lo, hi = float(v.min()), float(v.max())
        pad = max((hi - lo) * f, 0.25)
        return lo - pad, hi + pad
    zx, zy = lim(u.loc[u["in_zoom"], "umap___1"]), lim(u.loc[u["in_zoom"], "umap___2"])
    IN = (0.60, 0.58, 0.38, 0.40)
    x_in = zx[0] + IN[0] * (zx[1] - zx[0])
    y_need = u.loc[u["in_zoom"] & (u["umap___1"] >= x_in - 0.3), "umap___2"].max() + 0.3
    zy = (zy[0], max(zy[1], zy[0] + (y_need - zy[0]) / (IN[1] - 0.02)))
    fx, fy = lim(u["umap___1"]), lim(u["umap___2"])
    b = block("morphology_umap")
    left, gap = 9.0, 1.5
    aw = (b.w - left - 3 * gap - 0.5) / 4
    for k, cy in enumerate(["Male", "Proestrus", "Estrus", "Diestrus"]):
        ax = b.ax(left + k * (aw + gap), 3.5, aw, b.h - 3.5 - 7.0)
        if k == 0:
            first = ax
        sub = u[u["Phase"] == cy]
        ins = ax.inset_axes(list(IN))
        for g in ("Control", "CCI"):
            s_ = sub[sub["Condition"] == g]
            zs = s_[s_["in_zoom"]]
            ax.scatter(zs["umap___1"], zs["umap___2"], s=2.5, color=COND[g], lw=0, alpha=0.9, zorder=3 if g == "CCI" else 2)
            ins.scatter(s_["umap___1"], s_["umap___2"], s=0.8, color=COND[g], lw=0, alpha=0.9)
        ins.set_xlim(fx); ins.set_ylim(fy); ins.set_xticks([]); ins.set_yticks([])
        for sp in ins.spines.values():
            sp.set_visible(True); sp.set_linewidth(0.4)
        ins.add_patch(Rectangle((zx[0], zy[0]), zx[1] - zx[0], zy[1] - zy[0], fill=False, ec="#333333", lw=0.4))
        ax.set_title(MM2.CYCLE_EN[cy], fontweight="bold", pad=1.5)
        ax.set_xlim(zx); ax.set_ylim(zy)
        ax.set_xlabel("UMAP 1")
        L.finish(ax)
        if k:
            ax.tick_params(axis="y", left=False, labelleft=False)
    first.set_ylabel("UMAP 2")
    b.save(OUT / "morphology_umap")


def panel_morphology_persistence_difference():
    from morphomics_plot_utils import compute_persistence_centroids, pick_representative_swc, read_swc_dataframe
    from scipy.ndimage import gaussian_filter
    data = ROOT / "data" / "MorphOMICs"
    with np.load(data / "saved_instances_cam" / "cam_persistence_vectors.npz", allow_pickle=False) as z:
        cols = ["Condition", "Phase", "Sex", "Side", "Layer"]
        dv = pd.DataFrame({c: z[c] for c in cols})
        dv["pi"] = list(z["persistence_images"])
    dv["Phase"] = dv["Phase"].replace({"Metestrus": "Diestrus"})
    cen = compute_persistence_centroids(dv, group_columns=["Condition", "Phase", "Side", "Layer"], vector_column="pi")
    man = pd.read_csv(data / "reports_cam" / "cam_manifest_selected.csv")
    man["Phase"] = man["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
    m12 = man[(man["side"].astype(str) == "IPSI") & (man["layer"].astype(str) == "12")]
    cmap = LinearSegmentedColormap.from_list("ctrl_cci", ["#3C3C3C", "#F4F4F4", COND["CCI"]])
    keyed = {(v["combo"]["Condition"], v["combo"]["Phase"]): v for v in cen.values()
             if str(v["combo"].get("Side")) == "IPSI" and str(v["combo"].get("Layer")) == "12"}
    b = block("morphology_persistence_difference")
    cw, chh = b.w / 2, b.h / 2
    for k, cy in enumerate(["Male", "Proestrus", "Estrus", "Diestrus"]):
        r, c = divmod(k, 2)
        ox, oy = c * cw, r * chh
        a, bb_ = keyed[("CCI", cy)], keyed[("Control", cy)]
        ai, bi = gaussian_filter(a["centroid"], 1.2), gaussian_filter(bb_["centroid"], 1.2)
        diff = gaussian_filter(ai / (ai.max() + 1e-12) - bi / (bi.max() + 1e-12), 1.0)
        vmax = float(np.abs(diff).max()) or 1.0
        side = chh - 4.0 - 7.5
        ax = b.ax(ox + 8.5, oy + 4.0, side, side)
        cax = b.ax(ox + 8.5 + side + 1.0, oy + 4.0, 1.6, side)
        im = ax.imshow(diff, origin="lower", cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_xlabel("Birth")
        ax.set_ylabel("Death")
        n = diff.shape[0]
        ticks = np.linspace(0, n - 1, 3)
        ax.set_xticks(ticks); ax.set_yticks(ticks)
        ax.set_xticklabels(["0", "40", "80"]); ax.set_yticklabels(["0", "40", "80"])
        ax.tick_params(labelsize=L.FS_TICK, length=L.TICK_LEN, width=0.5, pad=1)
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        cb = b.fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=L.FS_TICK - 0.5, length=1.5, width=0.4, pad=1)
        cb.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        cb.outline.set_linewidth(0.4)
        ph = MM2.CYCLE_EN[cy]
        b.text(ox + 8.5, oy + 3.2, f"{ph} CCI-IoN / {ph} C", ha="left", va="bottom", fontsize=L.FS_TITLE, fontweight="bold")
        tx0 = ox + 8.5 + side + 7.5
        tax = b.ax(tx0, oy + 5.0, cw - (tx0 - ox) - 0.5, side - 1.0)
        for cond, col in (("CCI", COND["CCI"]), ("Control", "#8C8C8C")):
            sub = m12[(m12["Condition"].astype(str) == cond) & (m12["Phase"] == cy)]
            paths = [Path(p) if Path(str(p)).is_absolute() else ROOT / str(p) for p in sub["target_path"].dropna()]
            sw = read_swc_dataframe(pick_representative_swc([p for p in paths if p.exists()]))
            ctr = sw.loc[sw["parent"] == -1, ["x", "y"]].iloc[0].to_numpy(float) if (sw["parent"] == -1).any() \
                else sw[["x", "y"]].mean().to_numpy(float)
            pts = {int(q.n): np.array([q.x, q.y]) - ctr for q in sw.itertuples(index=False)}
            for q in sw.itertuples(index=False):
                if int(q.parent) in pts:
                    p1, p2 = pts[int(q.n)], pts[int(q.parent)]
                    tax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=col, lw=0.6)
        tax.set_aspect("equal"); tax.axis("off")
    hs = [Line2D([], [], color=COND["CCI"], lw=1.6, label="CCI-IoN"), Line2D([], [], color="#8C8C8C", lw=1.6, label="Control")]
    lg = b.fig.legend(handles=hs, loc="lower right", ncol=2, bbox_to_anchor=(1.0, 0.0), handlelength=1.0,
                      borderaxespad=0.1)
    b.save(OUT / "morphology_persistence_difference")


def panel_morphology_ramification_pain_partial_spearman():
    corr = _load("corr", CODE / "MorphOMICs" / "plot_ramification_pain_correlation_layer12_english.py")
    src = GEN / "MorphOMICs" / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
    if not src.exists():
        src = ROOT / "data" / "MorphOMICs" / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
    d = pd.read_csv(src)
    d = d[(d["side"] == "IPSI") & (d["lamina"] == 12)][["animal_id", "condition", "sex_or_phase",
                                                         "ShollRamificationIndex", "test_0.16g"]].dropna()
    d["sex_or_phase"] = d["sex_or_phase"].replace({"Metestrus": "Diestrus"})
    d = d.drop_duplicates("animal_id").reset_index(drop=True)
    rho, p, dof, fx, py = corr.partial_spearman(d, "ShollRamificationIndex", "test_0.16g", ["sex_or_phase"])
    marks = {"Male": "s", "Proestrus": "o", "Estrus": "^", "Diestrus": "D"}
    b = block("morphology_ramification_pain_partial_spearman")
    ax = b.ax(10.0, 2.0, b.w - 10.0 - 24.0, b.h - 2.0 - 8.0)
    ax.grid(True, color="#E3E3E3", lw=0.5, zorder=0)
    bb, a = np.polyfit(fx, py, 1)
    xs = np.linspace(fx.min(), fx.max(), 50)
    ax.plot(xs, a + bb * xs, color="#202020", lw=0.8, zorder=2)
    for i, r in d.iterrows():
        ax.scatter(fx[i], py[i], marker=marks[r["sex_or_phase"]], s=L.SCATTER + 4,
                   facecolor=COND["CCI"] if r["condition"] == "CCI" else COND["Control"], edgecolor="#333333",
                   lw=0.4, zorder=3)
    ax.set_xlabel("Ramification residual rank")
    ax.set_ylabel("Pain residual rank")
    L.finish(ax)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=COND[c], mec="#333333", mew=0.4, ms=L.LEG_MS, label=l)
          for c, l in (("Control", "Control"), ("CCI", "CCI-IoN"))]
    h2 = [Line2D([0], [0], marker=m, ls="", mfc="white", mec="black", mew=0.5, ms=L.LEG_MS, label=k) for k, m in marks.items()]
    _two_legends(ax, h1, "Condition", h2, "Sex/estrous phase", 0.62)
    b.save(OUT / "morphology_ramification_pain_partial_spearman")
    return rho, p, len(d)


# ============================================================================ Candidate-window map
def panel_medip_candidate_window_map():
    mp = _load("medip_figure_panels", CODE / "MeDIP" / "figure_panels.py")
    d = pd.read_csv(GEN / "Figures_v2" / "panels" / "medip_candidate_window_map_display.csv")
    comps = ["F_CSF_vs_F_C", "M_CSF_vs_M_C", "F_C_vs_M_C", "F_CSF_vs_M_CSF"]
    markers = dict(zip(comps, ["o", "^", "s", "D"]))
    regions = ["TSS", "EXON", "GENEBODY"]
    genes = sorted(d["gene_label"].unique(), key=lambda s: s.lower(), reverse=True)
    b = block("medip_candidate_window_map")
    left = max(text_mm(g, L.FS_LABEL, italic=True)[0] for g in genes) + 2.0
    ax = b.ax(left, 1.0, b.w - left - 1.0, b.h - 1.0 - 17.0)
    rng = np.random.default_rng(7391)
    xm = {r: i for i, r in enumerate(regions)}
    ymap = {g: i for i, g in enumerate(genes)}
    for _, r in d.iterrows():
        if r["region"] not in xm:
            continue
        ax.scatter(xm[r["region"]] + rng.uniform(-0.2, 0.2), ymap[r["gene_label"]] + rng.uniform(-0.12, 0.12), s=8,
                   color=st.HYPER if r["direction_display"] == "Hyper" else st.HYPO, marker=markers[r["comparison"]],
                   lw=0, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(regions, fontsize=L.FS_LABEL, fontweight="bold")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=L.FS_LABEL, fontstyle="italic")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.tick_params(length=0, pad=1)
    ax.grid(axis="y", color="#E6E6E6", lw=0.4, zorder=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=L.LEG_MS, label=k)
          for k, c in (("Hyper", st.HYPER), ("Hypo", st.HYPO))]
    h2 = [Line2D([0], [0], marker=markers[c], ls="", mfc="#777777", mec="none", ms=L.LEG_MS, label=st.medip_label(c))
          for c in comps]
    l1 = b.fig.legend(handles=h1, title="Methylation", loc="upper left", bbox_to_anchor=(b.fx(1.0), b.fy(b.h - 10.8)))
    l1._legend_box.align = "left"
    l2 = b.fig.legend(handles=h2, title="Comparison", loc="upper left", ncol=2,
                      bbox_to_anchor=(b.fx(20.0), b.fy(b.h - 10.8)))
    l2._legend_box.align = "left"
    b.save(OUT / "medip_candidate_window_map")


# ============================================================================ main
def main():
    L.rc()
    OUT.mkdir(parents=True, exist_ok=True)
    vitro = _gate("in vitro RT-qPCR blocks", QPCR_VITRO)
    vivo = _gate("in vivo RT-qPCR blocks", QPCR_VIVO, VON_FREY)
    behaviour = _gate("behaviour block", BEHAVIOUR)
    pain_morph = _gate("ramification-pain block", VON_FREY, LINKAGE)
    # The first blocks use the block style set by L.rc(); every later block is
    # drawn after the shared journal style (journal_style.set_style) is applied.
    panel_medip_dmr_directionality(); panel_medip_go_rank_distribution()
    if vitro and vivo:
        bubbles()
    panel_medip_candidate_gene_heatmap()
    importlib.import_module("journal_style").set_style()
    panel_medip_csf1r_locus_dmr_ere_map(); panel_medip_candidate_window_map()
    panel_morphology_ramification_index(); panel_morphology_umap(); panel_morphology_persistence_difference()
    if vitro:
        panel_invitro_qpcr_csf1r_violin(); panel_invitro_qpcr_methylation_regulators_csf1_violins()
        panel_invitro_estrogenic_score_and_pca()
        iv6, _ = _iv6()
        vitro_grid("invitro_qpcr_methylation_regulators_violins", list(iv6.METHYLATION_GENES), 3)
        vitro_grid("invitro_qpcr_estrogen_genes_violins", list(iv6.ESTROGEN_GENES), 4)
    if behaviour:
        panel_behaviour_change_score_violins()
    if vivo:
        vivo_grid("invivo_qpcr_pain_microglia_markers", ["Slc12a5", "Csf1r", "Bdnf", "Tnf", "Il10", "Tlr4"], 3)
        vivo_grid("invivo_qpcr_methylation_regulators", ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45g"], 2, legend=False)
        vivo_grid("invivo_qpcr_gadd45ab_violins", ["Gadd45a", "Gadd45b"], 2, legend=True)
        vivo_grid("invivo_qpcr_estrogen_genes_violins", ["Esr1", "Esr2", "Gper1", "Cyp19a1"], 2)
        panel_invivo_marker_pain_spearman(); panel_estrogenic_score_vs_pain(); panel_invivo_mfa_phase_condition()
    if pain_morph:
        rho, p, n = panel_morphology_ramification_pain_partial_spearman()
        print(f"[OK] ramification-pain partial rho={rho:.3f}, P={p:.4f}, n={n}")
    print(f"[OK] figure blocks -> {OUT}")


if __name__ == "__main__":
    main()
