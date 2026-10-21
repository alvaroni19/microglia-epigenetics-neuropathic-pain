"""Regenerate the figure panels from the release tables.

No statistic is recomputed: log2FC, P values, stars and contrasts are read
from the tables written by the release pipeline (code/reproduce_all_figures.py).
The only derived quantities are display transformations:
  * MeDIP contrasts written CSF1 first (F-CSF1 vs F-C, M-CSF1 vs M-C) and
    male first (M-C vs F-C, M-CSF1 vs F-CSF1); sex contrasts and control-first
    keys have the sign of log2FC inverted and Hyper/Hypo exchanged (see
    original_style.py);
  * in vivo female fold changes expressed relative to the diestrus-control
    mean dCt (males relative to the male-control mean dCt), which rescales
    each panel without changing any test;
  * genome-wide ESR1/ESR2 PWM hits along the rn6 Csf1r locus (Csf1r locus track),
    using the same matrices and threshold (relative score >= 0.92) as the
    window-level scan.

Panels are written to generated_outputs/Figures_v2/panels/ as SVG, PDF and
600-dpi PNG.  build_figure_blocks.py redraws every data panel at its final
printed size (generated_outputs/Figures_v2/blocks/).
"""
from __future__ import annotations

import glob
import importlib
import importlib.util
import shutil
import sys
from pathlib import Path
from textwrap import wrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).absolute().parents[2]
CODE = ROOT / "code"
GEN = ROOT / "generated_outputs"
OUT = GEN / "Figures_v2"
DATA_MEDIP = ROOT / "data" / "MeDIP" / "workspace_and_tables"
for sub in ("Figures_v2", "_helpers", "MeDIP", "PCR_260514"):
    sys.path.insert(0, str(CODE / sub))

for font in glob.glob(str(Path(sys.prefix) / "fonts" / "arial*.ttf")):
    font_manager.fontManager.addfont(font)

import journal_style as ns  # noqa: E402
import original_style as st  # noqa: E402

ns.save_panel = lambda fig, stem, *a, **k: st.save(fig, stem)

PHASES = ["Proestrus", "Estrus", "Diestrus"]
GREEK = {"Gadd45a": "Gadd45α", "Gadd45b": "Gadd45β", "Gadd45g": "Gadd45γ", "Tnf": "Tnf-α"}


def _style():
    ns.set_style()
    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42,
                         "axes.linewidth": 0.6})


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _stars(p):
    if not np.isfinite(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""


# ============================================================================ MeDIP overview panels
def panel_medip_dmr_directionality(mp, out_dir):
    """DMR directionality bars, manuscript colours and contrast labels."""
    stt = mp._directionality_stats().set_index("comparison")
    stt.reset_index().to_csv(out_dir / "medip_dmr_directionality_counts.csv", index=False)
    specs = [
        ("medip_dmr_directionality_csf1_effect", "CSF1 effect within sex",
         [("F_CSF_vs_F_C", "F-CSF1 vs\nF-C"), ("M_CSF_vs_M_C", "M-CSF1 vs\nM-C")],
         # left (blue) = hypomethylated after CSF1, right (orange) = hypermethylated after CSF1
         lambda r: (r["prop_hypo"], r["prop_hyper"]),
         (st.CSF1_HYPO, "Hypomethylated after CSF1"), (st.CSF1_HYPER, "Hypermethylated after CSF1")),
        ("medip_dmr_directionality_sex_effect", "Sex effect (male vs female)",
         [("F_C_vs_M_C", "M-C vs\nF-C"), ("F_CSF_vs_M_CSF", "M-CSF1 vs\nF-CSF1")],
         # release Hyper = higher in female; left (green) = higher in female
         lambda r: (r["prop_hyper"], r["prop_hypo"]),
         (st.SEX_FEMALE_HIGH, "Higher methylation in female"),
         (st.SEX_MALE_HIGH, "Higher methylation in male")),
    ]
    for stem, title, rows, split, left_spec, right_spec in specs:
        fig = plt.figure(figsize=(64 * st.MM, 38 * st.MM))
        ax = fig.add_axes([0.25, 0.42, 0.72, 0.45])
        for i, (key, lab) in enumerate(rows):
            lft, rgt = split(stt.loc[key])
            y = len(rows) - 1 - i
            ax.barh(y, -lft, color=left_spec[0], height=0.62)
            ax.barh(y, rgt, color=right_spec[0], height=0.62)
            for val, xc in ((lft, -lft / 2), (rgt, rgt / 2)):
                ax.text(xc, y, f"{100 * val:.0f}%", ha="center", va="center", color="white",
                        fontsize=7, fontweight="bold")
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xticklabels(["100%", "50%", "0", "50%", "100%"], fontsize=6)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([lab for _, lab in rows][::-1], fontsize=7, fontweight="bold")
        ax.set_xlabel("Fraction of DMR regions", fontsize=7, fontweight="bold", labelpad=1)
        ax.set_title(title, fontsize=7.5, fontweight="bold", pad=3)
        st.finish(ax)
        ax.tick_params(axis="y", length=0)
        ax.legend(handles=[Patch(fc=left_spec[0], label=left_spec[1]),
                           Patch(fc=right_spec[0], label=right_spec[1])],
                  loc="upper center", bbox_to_anchor=(0.42, -0.30), ncol=1, frameon=False,
                  fontsize=6.5, handlelength=0.9, handleheight=0.7, borderpad=0, labelspacing=0.15)
        st.save(fig, out_dir / stem)


OLD_GO_TERMS = {
    # (sex, displayed direction) -> {theme: preferred GO term to label}
    ("M", "Hypo"): {"Morphology": "negative regulation of cytoskeleton organization",
                     "Mechanical stimulus / pain": "response to mechanical stimulus",
                     "Inflammation / neuroinflammation": "positive regulation of phagocytosis, engulfment",
                     "Sex hormones / estrogen": "negative regulation of intracellular estrogen receptor signaling pathway"},
    ("F", "Hypo"): {"Morphology": "cytoskeletal motor activity",
                     "Inflammation / neuroinflammation": "regulation of mitotic cytokinesis",
                     "Sex hormones / estrogen": "estrogen response element binding"},
    ("M", "Hyper"): {"Morphology": "cytoskeletal motor activity",
                    "Mechanical stimulus / pain": "response to mechanical stimulus",
                    "Inflammation / neuroinflammation": "regulation of inflammatory response to wounding",
                    "Sex hormones / estrogen": "intracellular estrogen receptor signaling pathway"},
    ("F", "Hyper"): {"Morphology": "cytoskeletal motor activity",
                    "Sex hormones / estrogen": "nuclear estrogen receptor binding",
                    "Inflammation / neuroinflammation": "positive regulation of immune effector process"},
}


def panel_medip_go_rank_distribution(mp, out_dir):
    """GO rank distributions, contrasts written CSF1 first.

    Displayed 'Hyper' (higher methylation after CSF1 than in the same-sex
    control) is the release direction 'Hyper' of the CSF1-vs-control contrast.
    Labels: the preferred term is shown when it is enriched
    in that panel; otherwise the release representative term of the same
    theme is shown.
    """
    from collections import Counter
    df = pd.read_excel(DATA_MEDIP / "GO_relevance_rank_distribution_ws400.xlsx")
    df["biological_theme"] = df["biological_theme"].replace({"Steroid / sex hormones": "Sex hormones / estrogen"})
    df = df[np.isfinite(df["rank_within_group"]) & np.isfinite(df["neglog10_pvalue"])].copy()
    tokens = df["Description"].map(mp._label_tokens)
    doc_freq = Counter(t for ts in tokens for t in ts)
    n_docs = len(df)
    panels = [("M", "GO M CSF1 vs M C", "M-CSF1 vs M-C"), ("F", "GO F CSF1 vs F C", "F-CSF1 vs F-C")]
    go_key = {"M": "M_CSF_vs_M_C", "F": "F_CSF_vs_F_C"}
    fig = plt.figure(figsize=(104 * st.MM, 84 * st.MM))
    gs = fig.add_gridspec(2, 2, wspace=0.22, hspace=0.30, left=0.10, right=0.99, top=0.93, bottom=0.15)
    used = []
    for r, disp_dir in enumerate(["Hyper", "Hypo"]):
        for c, (sex, rel_comp, title) in enumerate(panels):
            rel_dir = st.medip_direction(go_key[sex], disp_dir)
            ax = fig.add_subplot(gs[r, c])
            sub = df[(df["comparison_label"] == rel_comp) & (df["direction"] == rel_dir)].copy()
            ax.scatter(sub["rank_within_group"], sub["neglog10_pvalue"], s=1.2, color="#C9C9C9",
                       lw=0, zorder=2, rasterized=False)
            for theme, col in mp.GO_THEME_COLORS.items():
                th = sub[sub["biological_theme"] == theme]
                ax.scatter(th["rank_within_group"], th["neglog10_pvalue"], s=3.5, color=col, lw=0, zorder=3)
            picks = []
            for theme, term in OLD_GO_TERMS[(sex, disp_dir)].items():
                hit = sub[sub["Description"].str.lower() == term.lower()]
                if len(hit):
                    row = hit.sort_values("neglog10_pvalue", ascending=False).iloc[0]
                    picks.append((row, "preferred term"))
                else:
                    rep = mp._representative_theme_term(sub[sub["biological_theme"] == theme], theme,
                                                        doc_freq, n_docs)
                    if len(rep):
                        picks.append((rep.iloc[0], f"representative term ('{term}' not enriched)"))
            picks = sorted(picks, key=lambda p: -p[0]["neglog10_pvalue"])
            x_max = max(1.0, sub["rank_within_group"].max())
            y_top = max(1.0, sub["neglog10_pvalue"].max())
            ax.set_xlim(-0.02 * x_max, x_max * 1.04)
            ax.set_ylim(0, y_top * 1.12)
            lanes = np.linspace(y_top * 1.02, y_top * 0.30, max(len(picks), 1))
            for lane, (row, src) in zip(lanes, picks):
                col = mp.GO_THEME_COLORS.get(row["biological_theme"], "#666666")
                ax.annotate("\n".join(wrap(str(row["Description"]), 30)),
                            xy=(row["rank_within_group"], row["neglog10_pvalue"]),
                            xytext=(x_max * 0.30, lane), fontsize=5.6, color=col, ha="left", va="center",
                            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec=col, lw=0.5),
                            arrowprops=dict(arrowstyle="-", color="#777777", lw=0.35, ls=(0, (2, 1.5))),
                            zorder=6)
                used.append(dict(panel=f"{title} {disp_dir}", release_contrast=rel_comp,
                                 release_direction=rel_dir, term=row["Description"],
                                 theme=row["biological_theme"], region=row["region"],
                                 rank=int(row["rank_within_group"]),
                                 neglog10_p=round(float(row["neglog10_pvalue"]), 3), source=src))
            if r == 0:
                ax.set_title(title, fontsize=7, fontweight="bold", pad=3)
            if c == 0:
                ax.set_ylabel(f"{disp_dir}\n-log$_{{10}}$(P value)", fontsize=6.5, fontweight="bold", labelpad=1)
            if r == 1:
                ax.set_xlabel("Rank", fontsize=6.5, fontweight="bold", labelpad=1)
            ax.grid(True, color="#EEEEEE", lw=0.3, zorder=0)
            st.finish(ax, labelsize=5.5)
    handles = [Line2D([0], [0], marker="o", ls="", mfc=col, mec="none", ms=3.5, label=mp.GO_THEME_LEGEND[t])
               for t, col in mp.GO_THEME_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=6,
               title="Biological theme", title_fontsize=6, bbox_to_anchor=(0.54, -0.01),
               handletextpad=0.1, columnspacing=0.8)
    st.save(fig, out_dir / "medip_go_rank_distribution")
    pd.DataFrame(used).to_csv(out_dir / "medip_go_rank_distribution_labels_used.csv", index=False)


CANDIDATE_GENE_ORDER = ["Bdnf", "Esr1", "Irf3", "Tnf", "Irf5", "Irf9", "Cyp19a1", "Irf4", "Esr2", "Irf6",
               "Spi1", "Dnmt1", "Dnmt3b", "Csf1r", "Il10", "Irf2", "Dnmt3a", "Irf1", "Irf7", "Irf8"]


def panel_medip_candidate_gene_heatmap(out_dir):
    """Candidate-gene heatmap: mean window log2FC, contrasts CSF1 first."""
    hits = pd.read_excel(DATA_MEDIP / "methylation_map_candidate_genes_hits_ws400.xlsx")
    keep = ["F_CSF_vs_F_C", "M_CSF_vs_M_C"]
    hits = hits[hits["comparison"].isin(keep)].copy()
    hits["log2fc_display"] = hits["edgeR.logFC"] * hits["comparison"].map(lambda k: st.MEDIP_CONTRAST[k][1])
    gc = hits.groupby(["symbol", "comparison"]).agg(v=("log2fc_display", "mean"),
                                                     p=("edgeR.p.value", "min"),
                                                     n_windows=("log2fc_display", "size")).reset_index()
    # genes without a DMR window in either within-sex contrast are not shown
    order = [g for g in CANDIDATE_GENE_ORDER if g in set(gc["symbol"])]
    z = gc.pivot(index="symbol", columns="comparison", values="v").reindex(index=order, columns=keep)
    pv = gc.pivot(index="symbol", columns="comparison", values="p").reindex(index=order, columns=keep)
    tab = gc.copy()
    tab["contrast"] = tab["comparison"].map(st.medip_label)
    tab["stars"] = tab["p"].map(_stars)
    tab.rename(columns={"v": "mean_log2FC_display", "p": "min_window_P"}).to_csv(
        out_dir / "medip_candidate_gene_heatmap_values.csv", index=False)
    M = z.to_numpy(float)
    vmax = float(np.nanmax(np.abs(M)))
    nr, nc = M.shape
    cell_w, cell_h = 11.0, 3.4
    left, bottom = 14.0, 13.0
    fw, fh = left + nc * cell_w + 18.0, bottom + nr * cell_h + 2.0
    fig = plt.figure(figsize=(fw * st.MM, fh * st.MM))
    ax = fig.add_axes([left / fw, bottom / fh, nc * cell_w / fw, nr * cell_h / fh])
    ax.add_patch(Rectangle((-0.5, -0.5), nc, nr, facecolor="#F7F7F7", edgecolor="none", zorder=0))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    ax.pcolormesh(np.arange(nc + 1) - 0.5, np.arange(nr + 1) - 0.5, np.ma.masked_invalid(M),
                  cmap=st.LOGFC_CMAP, norm=norm, shading="flat", edgecolors="none", zorder=1)
    for i in range(nr):
        for j in range(nc):
            s = _stars(float(pv.iloc[i, j])) if np.isfinite(pv.iloc[i, j]) else ""
            if s:
                ax.text(j, i, s, ha="center", va="center", fontsize=6.5, fontweight="bold",
                        color="white" if abs(M[i, j]) > 0.6 * vmax else "black", zorder=3)
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.set_yticks(range(nr))
    ax.set_yticklabels([GREEK.get(g, g) for g in order], fontsize=7, fontstyle="italic")
    ax.set_xticks(range(nc))
    ax.set_xticklabels([st.medip_label(k) for k in keep], rotation=30, ha="right",
                       rotation_mode="anchor", fontsize=6.5)
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = fig.add_axes([(left + nc * cell_w + 2.0) / fw, (bottom + 0.1 * nr * cell_h) / fh,
                        2.6 / fw, 0.8 * nr * cell_h / fh])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=st.LOGFC_CMAP), cax=cax)
    cb.solids.set_rasterized(False)
    cb.set_label("log$_2$ fold change", fontsize=7, fontweight="bold", labelpad=2)
    cb.ax.tick_params(labelsize=6, width=0.4, length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)
    st.save(fig, out_dir / "medip_candidate_gene_heatmap")


# ============================================================================ Csf1r locus map
def panel_medip_csf1r_locus(out_dir):
    import logomaker
    import map_csf1r_esr_motifs_integrated as mc
    gene_start, gene_stop, strand, exons = mc.read_gene_model()
    rec = mc.read_fasta(DATA_MEDIP / "csf1r_context_rn6.fasta")[0]
    x0, x1 = int(rec.start), int(rec.start) + len(rec.sequence) - 1
    # genome-wide ESR PWM hits on the context sequence
    motifs = [mc.read_meme(ROOT / "data" / "MeDIP" / "motifs" / f) for f in ("MA0112.3.meme", "MA0258.2.meme")]
    pwm_hits = pd.concat([mc.scan_record(rec, m, 0.92) for m in motifs], ignore_index=True)
    pwm_hits.to_csv(out_dir / "medip_csf1r_locus_esr_pwm_hits_rs092.csv", index=False)

    scan = pd.read_excel(DATA_MEDIP / "DMR_candidate_genes_ERE_scan_ws400.xlsx", sheet_name="ERE_scan")
    scan = scan[(scan["chr"].astype(str) == "chr18") & (scan["stop"] >= x0) & (scan["start"] <= x1)].copy()
    scan = scan.drop_duplicates(["comparison", "start", "stop"])
    scan["contrast"] = scan["comparison"].map(st.medip_label)
    scan["log2FC_display"] = scan["edgeR.logFC"] * scan["comparison"].map(lambda k: st.MEDIP_CONTRAST[k][1])
    scan["direction_display"] = [st.medip_direction(k, d) for k, d in zip(scan["comparison"], scan["methylation_direction"])]
    win = (scan.sort_values("log2FC_display", key=lambda s: s.abs(), ascending=False)
               .drop_duplicates(["start", "stop"]).sort_values("start").reset_index(drop=True))
    win["window_id"] = [f"W{i + 1}" for i in range(len(win))]
    wid = dict(zip(win["start"], win["window_id"]))
    scan["window_id"] = scan["start"].map(wid)
    scan["plotted_colour_from_this_contrast"] = scan.set_index(["start", "comparison"]).index.isin(
        win.set_index(["start", "comparison"]).index)
    cols = ["window_id", "chr", "start", "stop", "contrast", "direction_display", "log2FC_display",
            "edgeR.p.value", "plotted_colour_from_this_contrast"]
    scan.sort_values(["window_id", "contrast"])[cols].to_csv(out_dir / "medip_csf1r_locus_windows_display.csv", index=False)

    fig = plt.figure(figsize=(178 * st.MM, 92 * st.MM))
    gene_ax = fig.add_axes([0.06, 0.86, 0.92, 0.10])
    cpg_ax = fig.add_axes([0.06, 0.66, 0.92, 0.18], sharex=gene_ax)
    trk_ax = fig.add_axes([0.06, 0.52, 0.92, 0.12], sharex=gene_ax)
    gene_ax.set_xlim(x0, x1)
    gene_ax.add_patch(Rectangle((gene_start, -0.06), gene_stop - gene_start, 0.12, facecolor="#2CA02C",
                                edgecolor="none", zorder=2))
    for _, e in exons.iterrows():
        gene_ax.add_patch(Rectangle((e["start"], -0.35), max(int(e["stop"] - e["start"]), 120), 0.70,
                                    facecolor="#9B59B6", edgecolor="#6C3483", lw=0.3, zorder=3))
    tss = gene_start if strand == 1 else gene_stop
    gene_ax.add_patch(Rectangle((tss - 420, -0.42), 380, 0.84, facecolor="#F39C12", edgecolor="none", zorder=4))
    gene_ax.text(tss - 230, -0.62, "TSS", ha="center", va="top", fontsize=7.5, color="#F39C12", fontweight="bold")
    gene_ax.text(gene_start + 0.06 * (gene_stop - gene_start), -0.62, "GENE BODY", ha="left", va="top",
                 fontsize=7.5, color="#2CA02C", fontweight="bold")
    gene_ax.text(gene_start + 0.24 * (gene_stop - gene_start), -0.62, "EXONS", ha="left", va="top",
                 fontsize=7.5, color="#9B59B6", fontweight="bold")
    gene_ax.text(gene_stop, 0.55, "Csf1r (rn6)", ha="right", va="bottom", fontsize=7, fontstyle="italic")
    gene_ax.set_ylim(-1.1, 0.9)
    gene_ax.axis("off")

    seq = rec.sequence.upper()
    step = 400
    starts = np.arange(0, len(seq), step)
    dens = [1000 * sum(seq[i + k:i + k + 2] == "CG" for k in range(min(step, len(seq) - i) - 1)) / max(1, min(step, len(seq) - i))
            for i in starts]
    cpg_ax.bar(x0 + starts, dens, width=step * 0.92, align="edge", color="#C8C8C8", edgecolor="none")
    cpg_ax.text(x0 + 0.5 * (x1 - x0), max(dens) * 0.95, "CpG / kb", ha="center", va="top", fontsize=7)
    cpg_ax.set_yticks([0, int(round(max(dens) / 10.0) * 10) // 2, int(round(max(dens) / 10.0) * 10)])
    cpg_ax.spines[["top", "right", "bottom"]].set_visible(False)
    cpg_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    cpg_ax.tick_params(axis="y", labelsize=6.5, length=2, width=0.5)

    trk_ax.axhline(0.5, color="#BBBBBB", lw=1.2, zorder=1)
    for _, h in pwm_hits.iterrows():
        inside = ((win["start"] <= h["stop"]) & (win["stop"] >= h["start"])).any()
        trk_ax.add_patch(Rectangle((h["start"], 0.32), max(int(h["stop"] - h["start"]), 110), 0.36,
                                   facecolor="black" if inside else "white", edgecolor="black",
                                   lw=0.5, zorder=5 if inside else 3))
    for _, w in win.iterrows():
        col = st.HYPER if str(w["direction_display"]).upper().startswith("HYPER") else st.HYPO
        trk_ax.add_patch(Rectangle((w["start"], 0.25), w["stop"] - w["start"], 0.50, facecolor="white",
                                   edgecolor=col, lw=1.0, zorder=4))
        mid = (w["start"] + w["stop"]) / 2
        nb = win[(win["start"] - w["start"]).abs().between(1, 1300)]
        dx = 0 if nb.empty else (-450 if (nb["start"] > w["start"]).any() else 450)
        trk_ax.text(mid + dx, 0.10, w["window_id"], ha="center", va="top", fontsize=7.5, color=col,
                    fontweight="bold")
    trk_ax.set_ylim(-0.55, 1.0)
    trk_ax.set_yticks([])
    trk_ax.spines[["top", "right", "left"]].set_visible(False)
    trk_ax.ticklabel_format(style="plain", axis="x", useOffset=False)
    trk_ax.tick_params(axis="x", labelsize=7, length=2, width=0.5)

    fig.text(0.50, 0.425, "MeDIP 400 bp windows", ha="right", va="center", fontsize=7.5, fontweight="bold")
    leg = fig.legend(handles=[Patch(fc="white", ec=st.HYPO, lw=1.0, label="Hypomethylated DMR"),
                              Patch(fc="white", ec=st.HYPER, lw=1.0, label="Hypermethylated DMR")],
                     loc="center left", bbox_to_anchor=(0.505, 0.425), ncol=2, frameon=False,
                     fontsize=7.5, handlelength=0.9, handletextpad=0.3, columnspacing=0.8)
    fig.text(0.50, 0.375, "ERE", ha="right", va="center", fontsize=7.5, fontweight="bold")
    t_ere = fig.text(0.505, 0.375, "PWM relative score ≥ 0.92", ha="left", va="center", fontsize=7.5)
    bb = t_ere.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.patches.append(Rectangle((bb.x1 + 0.006, 0.362), 0.009, 0.026, transform=fig.transFigure,
                                 fc="white", ec="black", lw=0.6))
    for i, (m, name) in enumerate(zip(motifs, ["ERα", "ERβ"])):
        lax = fig.add_axes([0.20 + i * 0.36, 0.03, 0.26, 0.24])
        prob = pd.DataFrame(m.pwm, columns=list("ACGT"))
        info = logomaker.transform_matrix(prob, from_type="probability", to_type="information")
        logomaker.Logo(info, ax=lax, color_scheme="classic", font_name="Arial", vpad=0.02)
        lax.set_ylim(0, 2)
        lax.set_ylabel("Bits", fontsize=7)
        lax.set_xticks([])
        lax.tick_params(axis="y", labelsize=6.5, length=2, width=0.5)
        lax.spines[["top", "right", "bottom"]].set_visible(False)
        lax.text(0.62, 1.06, name, transform=lax.transAxes, fontsize=8.5, ha="left", va="bottom")
        fig.add_artist(matplotlib.patches.FancyArrowPatch(
            (0.50, 0.345), (0.20 + i * 0.36 + 0.15, 0.30), transform=fig.transFigure,
            arrowstyle="-|>", mutation_scale=7, lw=0.7, color="black"))
    st.save(fig, out_dir / "medip_csf1r_locus_dmr_ere_map")
    return win, pwm_hits


# ============================================================================ violins
VIVO_ALIAS = {"Slc12a5": "KCC2", "Bdnf": "BDNF", "Csf1r": "CSF1R", "Tlr4": "TLR4", "Tnf": "TNFa",
              "Il10": "IL10", "Esr1": "ERa", "Esr2": "ERb", "Gper1": "GPER", "Cyp19a1": "ARO"}
TITLE = {"CSF1R": "Csf1r", "ERa": "Esr1", "ERb": "Esr2", "GPER": "Gper1", "ARO": "Cyp19a1", **GREEK}


def invivo_fold_changes(animal):
    """Males relative to the male-control mean dCt; females relative to the
    diestrus-control mean dCt (manuscript display)."""
    out = animal[["animal_id", "phase", "condition"]].copy()
    genes = [c[:-9] for c in animal.columns if c.endswith("_delta_ct")]
    for g in genes:
        dct = animal[f"{g}_delta_ct"]
        mref = dct[(animal.phase == "Male") & (animal.condition == "Control")].mean()
        fref = dct[(animal.phase == "Diestrus") & (animal.condition == "Control")].mean()
        ref = np.where(animal.phase == "Male", mref, fref)
        out[g] = np.power(2.0, -(dct - ref))
    return out


def _bracket_level(his, i, span):
    return his[i] + 0.03 * span


def invivo_gene_panel(gene, fc, stats):
    alias = VIVO_ALIAS.get(gene, gene)
    stx = stats[stats.gene == alias].set_index("comparison")["stars"]
    fig, (axm, axf) = plt.subplots(1, 2, figsize=(72 * st.MM, 42 * st.MM),
                                   gridspec_kw={"width_ratios": [1, 2.75], "wspace": 0.30})
    for ax, groups in ((axm, ["Male"]), (axf, PHASES)):
        his, ticks = [], []
        for i, ph in enumerate(groups):
            xc, xi = i * 2.3, i * 2.3 + 1
            sub = fc[fc.phase == ph]
            c = sub[sub.condition == "Control"][gene].dropna().values
            t = sub[sub.condition == "CCI"][gene].dropna().values
            _, h1 = st.violin(ax, c, xc, st.LIGHT[ph], width=0.85)
            _, h2 = st.violin(ax, t, xi, st.SAT[ph], width=0.85)
            his.append(np.nanmax([h1, h2]))
            ticks += [xc, xi]
        ymin, ymax = ax.get_ylim()
        span = ymax - ymin
        for i, ph in enumerate(groups):
            s = stx.get(f"{ph} CCI vs control")
            if isinstance(s, str) and s:
                st.bracket(ax, i * 2.3, i * 2.3 + 1, his[i] + 0.03 * span, 0.03 * span, s)
        ax.set_ylim(ymin, max(ymax, max(his) + 0.16 * span))
        ax.set_xlim(-0.65, (len(groups) - 1) * 2.3 + 1.65)
        ax.set_xticks(ticks)
        ax.set_xticklabels(["Control", "CCI-IoN"] * len(groups), rotation=45, ha="right",
                           rotation_mode="anchor", fontsize=6)
        ax.set_title(TITLE.get(gene, gene), fontstyle="italic", fontsize=7, pad=2)
        st.finish(ax)
    axm.set_ylabel("Fold change", fontweight="bold", fontsize=7, labelpad=1)
    return fig


def invitro_panel(gene, data, brackets, order, labels, figsize=(42, 38)):
    fig, ax = plt.subplots(figsize=(figsize[0] * st.MM, figsize[1] * st.MM))
    pos = {g: i + 1.0 for i, g in enumerate(order)}
    his = {}
    for g in order:
        sex = "male" if g.startswith("male") or g.endswith("_male") else "female"
        treat = ("CSF1E2" if "E2" in g else "CSF1" if "csf1" in g.lower() else "control")
        _, his[g] = st.violin(ax, data[g], pos[g], st.TREAT_FILL[treat], edge=st.SEX_EDGE[sex],
                              lw=1.3, width=0.82)
    ymin, ymax = ax.get_ylim()
    span = max(his.values()) - ymin
    tiers = []
    for g1, g2, s in sorted(brackets, key=lambda b: abs(pos[b[1]] - pos[b[0]])):
        a, b = sorted((pos[g1], pos[g2]))
        level = max(his[g] for g in order if a <= pos[g] <= b) + 0.03 * span
        for (ta, tb, ty) in tiers:
            if not (b < ta or a > tb):
                level = max(level, ty + 0.11 * span)
        tiers.append((a, b, level))
        st.bracket(ax, a, b, level, 0.03 * span, s)
    top = max([t[2] for t in tiers], default=max(his.values())) + 0.10 * span
    ax.set_ylim(ymin, max(ymax, top))
    ax.set_xlim(0.4, len(order) + 0.6)
    ax.set_xticks(list(pos.values()))
    ax.set_xticklabels([labels[g] for g in order], rotation=45, ha="right", rotation_mode="anchor", fontsize=6)
    ax.set_title(TITLE.get(gene, gene), fontstyle="italic", fontsize=7, pad=2)
    ax.set_ylabel("Fold change", fontweight="bold", fontsize=7, labelpad=1)
    st.finish(ax)
    return fig


def legend_panel(stem, kind):
    fig, ax = plt.subplots(figsize=(90 * st.MM, 6 * st.MM))
    ax.axis("off")
    if kind == "phase":
        hs = [Line2D([], [], marker="o", ls="", mfc=st.SAT[k], mec="black", mew=0.5, ms=6, label=k)
              for k in ["Male"] + PHASES]
        title = "Biological phase"
    else:
        hs = [Line2D([], [], marker="o", ls="", mfc="white", mec=st.SEX_EDGE[k], mew=1.6, ms=6,
                     label=k.capitalize()) for k in ("male", "female")]
        title = "Sex"
    hs = [Line2D([], [], ls="", marker="", label=title)] + hs
    ax.legend(handles=hs, ncol=len(hs), loc="center", frameon=False, fontsize=7,
              handletextpad=0.1, columnspacing=0.7, handlelength=0.9)
    st.save(fig, stem)


# ============================================================================ candidate-window map
def panel_medip_candidate_window_map(mp, out_dir):
    df = mp._load_sig()
    cg = pd.read_excel(DATA_MEDIP / "DMR_candidate_genes_ERE_scan_ws400.xlsx", sheet_name="candidate_genes")[["symbol", "family"]]
    cg["sl"] = cg["symbol"].str.lower()
    cg = cg[~cg["sl"].isin({f"irf{i}" for i in range(1, 10)})]
    sym = dict(zip(cg["sl"], cg["symbol"]))
    d = df[["chr", "start", "stop", "external_transcript_gene_name", "region", "diffexpressed", "comparison",
            "edgeR.logFC", "edgeR.p.value"]].copy()
    d["sl"] = d["external_transcript_gene_name"].astype(str).str.lower()
    d = d[d["sl"].isin(set(cg["sl"])) & d["diffexpressed"].isin(["Hyper", "Hypo"])]
    d = mp._reannotate_candidate_windows(d, set(cg["sl"]))
    d["region"] = d["current_region"]
    d = d.drop_duplicates(["chr", "start", "stop", "sl", "diffexpressed", "comparison", "region"])
    d["gene_label"] = d["sl"].map(lambda s: GREEK.get(sym[s][0].upper() + sym[s][1:].lower(),
                                                      sym[s][0].upper() + sym[s][1:].lower()))
    d["contrast"] = d["comparison"].map(st.medip_label)
    d["direction_display"] = [st.medip_direction(k, v) for k, v in zip(d["comparison"], d["diffexpressed"])]
    d.to_csv(out_dir / "medip_candidate_window_map_display.csv", index=False)
    comps = ["F_CSF_vs_F_C", "M_CSF_vs_M_C", "F_C_vs_M_C", "F_CSF_vs_M_CSF"]
    markers = dict(zip(comps, ["o", "^", "s", "D"]))
    regions = ["TSS", "EXON", "GENEBODY"]
    genes = sorted(d["gene_label"].unique(), key=lambda s: s.lower(), reverse=True)
    rng = np.random.default_rng(7391)
    fig, ax = plt.subplots(figsize=(78 * st.MM, (12 + 4.4 * len(genes)) * st.MM))
    xm = {r: i for i, r in enumerate(regions)}
    ymap = {g: i for i, g in enumerate(genes)}
    for _, r in d.iterrows():
        if r["region"] not in xm:
            continue
        ax.scatter(xm[r["region"]] + rng.uniform(-0.2, 0.2), ymap[r["gene_label"]] + rng.uniform(-0.12, 0.12),
                   s=11, color=st.HYPER if r["direction_display"] == "Hyper" else st.HYPO,
                   marker=markers[r["comparison"]], lw=0, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(regions, fontsize=7.5, fontweight="bold")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6.5, fontstyle="italic")
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.tick_params(length=0, pad=1)
    ax.grid(axis="y", color="#E6E6E6", lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=c, mec="none", ms=3.5, label=k)
          for k, c in (("Hyper", st.HYPER), ("Hypo", st.HYPO))]
    h2 = [Line2D([0], [0], marker=markers[c], ls="", mfc="#777777", mec="none", ms=4, label=st.medip_label(c))
          for c in comps]
    l1 = ax.legend(handles=h1, title="Methylation", loc="upper left", bbox_to_anchor=(-0.02, -0.06),
                   frameon=False, fontsize=6.5, title_fontsize=7, handletextpad=0.1, labelspacing=0.1)
    ax.add_artist(l1)
    ax.legend(handles=h2, title="Comparison", loc="upper left", bbox_to_anchor=(0.28, -0.06), ncol=2,
              frameon=False, fontsize=6.5, title_fontsize=7, handletextpad=0.1, columnspacing=0.6,
              labelspacing=0.1)
    st.save(fig, out_dir / "medip_candidate_window_map")


# ============================================================================ main
def _bubble(fp, csv, out, comps, genes, labels):
    df = pd.read_csv(GEN / "PCR_260514" / csv)
    df["gene_label"] = df["gene_label"].map(fp._disp_gene).map(lambda g: GREEK.get(g, g)).replace({"Tnfα": "Tnf-α"})
    genes = [GREEK.get(fp._disp_gene(g), fp._disp_gene(g)) for g in genes]
    st.save(st.bubble_matrix(df, comps, genes, labels), OUT / out)


def _copy(src_glob, dst_dir):
    dst_dir.mkdir(parents=True, exist_ok=True)
    for f in glob.glob(str(src_glob)):
        shutil.copy2(f, dst_dir)


QPCR_VITRO = ROOT / "data" / "PCR_260514" / "invitro_qpcr_relative_expression.xlsx"
QPCR_VIVO = ROOT / "data" / "PCR_260514" / "invivo_qpcr_delta_ct.xlsx"
BEHAVIOUR = ROOT / "data" / "Allodynia_0_16" / "allodynia_change_scores.xlsx"
VON_FREY = ROOT / "data" / "Allodynia_0_16" / "allodynia_von_frey_scores.xlsx"
RESTRICTED = ("behaviour and RT-qPCR inputs are not distributed "
              "(available from the corresponding author on request)")


def _gate(label, *paths):
    """True when every input exists; otherwise report the panel group as skipped."""
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        print(f"[SKIP] {label}: {missing[0].name} not found; {RESTRICTED}")
    return not missing


def main():
    _style()
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "panels").mkdir(parents=True, exist_ok=True)

    mp = _load("medip_figure_panels", CODE / "MeDIP" / "figure_panels.py")
    panel_medip_dmr_directionality(mp, OUT / "panels")
    _style()
    panel_medip_go_rank_distribution(mp, OUT / "panels")
    panel_medip_candidate_gene_heatmap(OUT / "panels")
    panel_medip_csf1r_locus(OUT / "panels")
    panel_medip_candidate_window_map(mp, OUT / "panels")
    _style()

    have_vitro = _gate("in vitro RT-qPCR panels", QPCR_VITRO)
    have_vivo = _gate("in vivo RT-qPCR panels", QPCR_VIVO, VON_FREY)
    if have_vitro or have_vivo:
        fp = _load("pcr_figure_panels", CODE / "PCR_260514" / "figure_panels.py")
        vitro = ["Male CSF1 vs control", "Male CSF1+E2 vs CSF1", "Female CSF1 vs control", "Female CSF1+E2 vs CSF1"]
        vitro_labels = {"Male CSF1 vs control": "M-CSF1 vs M-C", "Male CSF1+E2 vs CSF1": "M-CSF1+E2 vs M-CSF1",
                        "Female CSF1 vs control": "F-CSF1 vs F-C", "Female CSF1+E2 vs CSF1": "F-CSF1+E2 vs F-CSF1"}
        vivo = ["Male CCI vs control", "Proestrus CCI vs control", "Estrus CCI vs control", "Diestrus CCI vs control"]
        vivo_labels = {"Male CCI vs control": "M-CCI-IoN vs M-C",
                       "Proestrus CCI vs control": "Proestrus-CCI-IoN vs Proestrus-C",
                       "Estrus CCI vs control": "Estrus-CCI-IoN vs Estrus-C",
                       "Diestrus CCI vs control": "Diestrus-CCI-IoN vs Diestrus-C"}
        meth = ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"]

    if have_vitro:
        _bubble(fp, "invitro_methylation_regulators_stats.csv", "panels/invitro_qpcr_methylation_regulators_bubble", vitro, meth, vitro_labels)
        _bubble(fp, "invitro_estrogen_genes_stats.csv", "panels/invitro_estrogen_genes_bubble", vitro, fp.ER_ARO_GENES, vitro_labels)
        iv6 = importlib.import_module("build_invitro_6group_violins")
        iv2 = importlib.import_module("build_invitro_2way_violins")
        _style()
        lab6 = {g: iv6.GROUP_LABELS[g].replace("Ctrl", "C") for g in iv6.GROUPS}
        for gene, dst in [("CSF1R", "panels/invitro_qpcr_csf1r_violin")] + \
                [(g, f"panels/invitro_qpcr_methylation_regulators_violins_{g}") for g in iv6.METHYLATION_GENES] + \
                [(g, f"panels/invitro_qpcr_estrogen_genes_violins_{TITLE.get(g, g)}") for g in iv6.ESTROGEN_GENES]:
            size = (62, 50) if gene == "CSF1R" else (42, 38)
            st.save(invitro_panel(gene, iv6.load_data(gene), iv6.load_brackets(gene), iv6.GROUPS, lab6, size), OUT / dst)
        order2 = ["control_male", "control_female", "csf1_male", "csf1_female"]   # manuscript order
        lab2 = {"control_male": "M-C", "csf1_male": "M-CSF1", "control_female": "F-C", "csf1_female": "F-CSF1"}
        for gene in iv2.METHYLATION_GENES:
            br = [(b[0], b[1], b[2]) for b in iv2.load_brackets(gene)]
            st.save(invitro_panel(gene, iv2.load_data(gene), br, order2, lab2, (34, 38)),
                    OUT / f"panels/invitro_qpcr_methylation_regulators_csf1_violins_{gene}")


    if have_vivo:
        _bubble(fp, "invivo_estrogen_genes_stats.csv", "panels/invivo_estrogen_genes_bubble", vivo, fp.ER_ARO_GENES, vivo_labels)
        animal = pd.read_csv(GEN / "PCR_260514" / "invivo_qpcr_animal_level.csv")
        fc = invivo_fold_changes(animal)
        fc.to_csv(OUT / "panels" / "invivo_fold_change_display_values.csv", index=False)
        stats = pd.concat([pd.read_csv(GEN / "PCR_260514" / f) for f in (
            "invivo_pain_microglia_markers_stats.csv", "invivo_methylation_regulators_stats.csv",
            "invivo_estrogen_genes_stats.csv")], ignore_index=True)
        groups = {"panels/invivo_qpcr_pain_microglia_markers": ["Slc12a5", "Csf1r", "Bdnf", "Tnf", "Il10", "Tlr4"],
                  "panels/invivo_qpcr_methylation_regulators": ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45g"],
                  "panels/invivo_qpcr_gadd45ab_violins": ["Gadd45a", "Gadd45b"],
                  "panels/invivo_qpcr_estrogen_genes_violins": ["Esr1", "Esr2", "Gper1", "Cyp19a1"]}
        for prefix, genes in groups.items():
            for g in genes:
                st.save(invivo_gene_panel(g, fc, stats), OUT / f"{prefix}_{g}")
        legend_panel(OUT / "panels" / "legend_biological_phase", "phase")
        legend_panel(OUT / "panels" / "legend_sex", "sex")

        # estrogenic score vs pain scatter
        d5 = pd.read_csv(GEN / "PCR_MFA" / "invivo_estrogenic_score_delta_pain_per_animal.csv")
        fig, ax = plt.subplots(figsize=(68 * st.MM, 46 * st.MM))
        x, y = d5["delta_estrogenic_score"].to_numpy(float), d5["delta_pain_016"].to_numpy(float)
        b, a = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, a + b * xs, color="black", lw=1.0, zorder=2)
        for ph in PHASES:
            s = d5[d5["phase"] == ph]
            ax.scatter(s["delta_estrogenic_score"], s["delta_pain_016"], s=12, color=st.SAT[ph],
                       edgecolor="black", lw=0.4, zorder=3, label=ph)
        ax.set_xlabel("Δ estrogenic expression score", fontweight="bold", fontsize=7, labelpad=1)
        ax.set_ylabel("Δ mechanical response (0.16 g)", fontweight="bold", fontsize=7, labelpad=1)
        ax.legend(title="Biological phase", fontsize=6.5, title_fontsize=7, frameon=False, loc="upper left",
                  bbox_to_anchor=(0.62, 1.05), handletextpad=0.1, labelspacing=0.2)
        st.finish(ax)
        st.save(fig, OUT / "panels" / "estrogenic_score_vs_pain")


    if _gate("behaviour panel", BEHAVIOUR):
        # behaviour panel (built by code/Behaviour/build_behaviour_violins.py)
        _copy(GEN / "Behaviour_reanalysis" / "behaviour_change_score_violins.*", OUT / "panels")
    # correlation, MFA, estrogenic-score and morphology panels
    _load("build_panels_v2_mfa_morphology", CODE / "Figures_v2" / "build_panels_v2_mfa_morphology.py").main()
    print("[OK] manuscript panels v2 ->", OUT)


if __name__ == "__main__":
    main()
