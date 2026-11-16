# -*- coding: utf-8 -*-
"""Build the MeDIP panels (heatmap, directionality, GO ranks, candidate-window map)."""
from __future__ import annotations

import re
import gzip
import math
import sys
from collections import Counter
from pathlib import Path
from textwrap import wrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

ROOT = Path(__file__).absolute().parents[2]
sys.path.append(str(ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402

DATA = ROOT / "data" / "MeDIP" / "workspace_and_tables"
OUT = ROOT / "generated_outputs" / "MeDIP"
SIG_XLSX = DATA / "GENES_TOTAL_sig_contig_colap_ws400.xlsx"
CAND_XLSX = DATA / "DMR_candidate_genes_ERE_scan_ws400.xlsx"
RN6_GTF = (
    ROOT / "data" / "MeDIP" / "annotations" / "rn6"
    / "Rattus_norvegicus.Rnor_6.0.104.gtf.gz"
)

COMPARISON_LABELS = {
    "F_CSF_vs_F_C": "F-CSF1 vs F-C", "M_CSF_vs_M_C": "M-CSF1 vs M-C",
    "F_C_vs_M_C": "F-C vs M-C", "F_CSF_vs_M_CSF": "F-CSF1 vs M-CSF1",
}


def _load_sig():
    return pd.read_excel(io_path(SIG_XLSX), sheet_name=0)


def _candidate_symbols():
    cg = pd.read_excel(io_path(CAND_XLSX), sheet_name="candidate_genes")
    return cg["symbol"].astype(str).tolist()


def _parse_gtf_attributes(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def _candidate_gene_models(symbols):
    """Read Ensembl release-104 rn6 gene spans and exon intervals."""
    models = {}
    with gzip.open(io_path(RN6_GTF), "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = _parse_gtf_attributes(fields[8])
            if fields[2] not in {"gene", "exon"}:
                continue
            symbol = attrs.get("gene_name", "").lower()
            if symbol not in symbols:
                continue
            model = models.setdefault(symbol, {"exons": []})
            if fields[2] == "gene":
                model.update(
                    chrom=fields[0].removeprefix("chr"),
                    start=int(fields[3]),
                    stop=int(fields[4]),
                    strand=fields[6],
                    gene_id=attrs.get("gene_id", ""),
                )
            else:
                model["exons"].append((int(fields[3]), int(fields[4])))
    return models


def _overlaps(start, stop, interval_start, interval_stop):
    return min(stop, interval_stop) >= max(start, interval_start)


def _reannotate_candidate_windows(df, candidate_symbols):
    """Assign regions against the same rn6 assembly used for alignment."""
    models = _candidate_gene_models(candidate_symbols)
    rows = []
    for row in df.itertuples(index=False):
        symbol = str(row.external_transcript_gene_name).lower()
        model = models.get(symbol)
        chrom = str(row.chr).removeprefix("chr")
        start, stop = int(row.start), int(row.stop)
        current_region = None
        status = "gene_not_resolved"
        if model and chrom == model.get("chrom"):
            tss = model["start"] if model["strand"] == "+" else model["stop"]
            if any(_overlaps(start, stop, a, b) for a, b in model["exons"]):
                current_region, status = "EXON", "retained"
            elif _overlaps(start, stop, tss - 1000, tss + 1000):
                current_region, status = "TSS", "retained"
            elif _overlaps(start, stop, model["start"], model["stop"]):
                current_region, status = "GENEBODY", "retained"
            else:
                status = "excluded_no_current_gene_overlap"
        elif model:
            status = "excluded_chromosome_mismatch"
        record = row._asdict()
        record.update(
            old_region=str(row.region).replace("GENE", "GENEBODY"),
            current_region=current_region,
            annotation_status=status,
            reference_gene_id=(model or {}).get("gene_id", ""),
            reference_gene_start=(model or {}).get("start", np.nan),
            reference_gene_stop=(model or {}).get("stop", np.nan),
            reference_gene_strand=(model or {}).get("strand", ""),
            reference_assembly="rn6 / Rnor_6.0",
            reference_annotation="Ensembl release 104",
        )
        rows.append(record)
    audit = pd.DataFrame(rows)
    audit.to_csv(io_path(OUT / "medip_candidate_window_reannotation_audit.csv"), index=False)
    return audit[audit["annotation_status"].eq("retained")].copy()


def _cluster_order(mat):
    """Average-linkage order on 1 - Pearson correlation (rows of `mat`)."""
    if mat.shape[0] < 3:
        return np.arange(mat.shape[0])
    cc = np.corrcoef(mat)
    cc[~np.isfinite(cc)] = 0.0
    np.fill_diagonal(cc, 1.0)
    d = 1.0 - cc
    d = (d + d.T) / 2.0
    np.fill_diagonal(d, 0.0)
    d[d < 0] = 0.0
    Z = linkage(squareform(d, checks=False), method="average")
    return leaves_list(Z)


# ----------------------------------------------------------------------------- candidate-gene heatmap
HITS_XLSX = DATA / "methylation_map_candidate_genes_hits_ws400.xlsx"
HEATMAP_ORDERS_XLSX = DATA / "heatmap_candidate_genes_logFC_ws400_orders.xlsx"


def build_candidate_gene_heatmap():
    """Candidate-gene mean edgeR logFC heatmap.

    Source: methylation_map_candidate_genes_hits_ws400.xlsx. Cells show mean
    window-level logFC and the row order is read from
    heatmap_candidate_genes_logFC_ws400_orders.xlsx.
    """
    ns.set_style()
    hits = pd.read_excel(io_path(HITS_XLSX))
    gc = (hits.groupby(["symbol", "comparison"], as_index=False)
          .agg(v=("edgeR.logFC", "mean"), p_value=("edgeR.p.value", "min"))
          .rename(columns={"symbol": "gene"}))
    wide = gc.pivot(index="gene", columns="comparison", values="v")
    p_wide = gc.pivot(index="gene", columns="comparison", values="p_value")
    keep = ["F_CSF_vs_F_C", "M_CSF_vs_M_C"]
    z = wide.reindex(columns=keep)
    p_wide = p_wide.reindex(index=z.index, columns=keep)
    # restore original R gene order
    orders = pd.read_excel(io_path(HEATMAP_ORDERS_XLSX), sheet_name="row_order").sort_values("order")
    gene_order = [g for g in orders["row_id"].tolist() if g in z.index]
    remaining = [g for g in z.index if g not in gene_order]
    z = z.loc[gene_order + remaining]
    p_wide = p_wide.loc[z.index]

    genes = list(z.index)
    cols = [COMPARISON_LABELS.get(c, c) for c in z.columns]
    M = np.ma.masked_invalid(z.to_numpy(dtype=float))
    vmax = float(np.nanmax(np.abs(M))) or 1.0

    fig, ax = ns.fig_mm(76, 86)
    cmap = ns.diverging_logfc_cmap()
    cmap.set_bad("white")
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    x_edges = np.arange(M.shape[1] + 1) - 0.5
    y_edges = np.arange(M.shape[0] + 1) - 0.5
    im = ax.pcolormesh(
        x_edges, y_edges, M, cmap=cmap, norm=norm, shading="flat",
        edgecolors="none",
    )
    ax.set_aspect("auto")
    ax.invert_yaxis()
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=ns.FS_TICK, rotation=30, ha="right")
    ax.set_yticks(range(len(genes)))
    # Title-case italic gene symbols; Tnf is displayed as Tnfalpha.
    gene_disp = ["Tnfα" if str(g).lower() == "tnf" else str(g) for g in genes]
    ax.set_yticklabels(gene_disp, fontsize=ns.FS_AXIS, fontstyle="italic")
    ax.tick_params(length=0)
    def significance_stars(p_value):
        if not np.isfinite(p_value):
            return ""
        if p_value < 1e-4:
            return "****"
        if p_value < 1e-3:
            return "***"
        if p_value < 1e-2:
            return "**"
        if p_value < 5e-2:
            return "*"
        return ""
    for row_index, gene in enumerate(z.index):
        for column_index, comparison in enumerate(z.columns):
            stars = significance_stars(float(p_wide.loc[gene, comparison]))
            if stars:
                value = float(z.loc[gene, comparison])
                text_color = "white" if abs(value) > 0.55 * vmax else "black"
                ax.text(column_index, row_index, stars, ha="center", va="center",
                        fontsize=ns.FS_SMALL, fontweight="bold", color=text_color)
    for sp in ax.spines.values():
        sp.set_linewidth(0.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.03, aspect=18)
    cb.solids.set_rasterized(False)
    cb.set_label("Mean log$_2$FC", fontsize=ns.FS_AXIS, fontweight="bold")
    cb.ax.tick_params(labelsize=ns.FS_TICK, width=0.4, length=2)
    cb.outline.set_linewidth(0.4)
    ax.set_title("Candidate-gene mean edgeR log$_2$FC",
                 fontsize=ns.FS_TITLE, fontweight="bold")
    ns.save_panel(fig, OUT / "medip_candidate_gene_heatmap_logFC", vector=True, tiff=True)
    print(f"[OK] candidate-gene heatmap ({len(genes)} genes x {len(cols)} contrasts)")


# ----------------------------------------------------------------------------- DMR directionality
def _collapse_regions(df):
    """Collapse contiguous 400 bp windows of the same direction into DMR regions."""
    df = df.sort_values(["comparison", "chr", "direction", "start"]).copy()
    key = list(zip(df["comparison"], df["chr"], df["direction"]))
    prev_stop = df["stop"].shift(1)
    prev_key = pd.Series(key).shift(1)
    cur_key = pd.Series(key, index=df.index)
    brk = (prev_key.values != pd.Series(key).values) | (df["start"].values > (prev_stop.values + 1))
    df["region_id"] = np.cumsum(brk)
    reg = df.groupby(["comparison", "chr", "direction", "region_id"], sort=False).agg(
        start=("start", "min"), stop=("stop", "max")).reset_index()
    return reg


def _directionality_stats():
    df = _load_sig()[["comparison", "chr", "start", "stop", "MSets1.rms.mean", "MSets2.rms.mean"]].copy()
    df = df.rename(columns={"MSets1.rms.mean": "ms1", "MSets2.rms.mean": "ms2"})
    df["direction"] = np.where(df["ms1"] > df["ms2"], "Hyper",
                               np.where(df["ms1"] < df["ms2"], "Hypo", None))
    df = df.dropna(subset=["direction"])
    df = df.drop_duplicates(["comparison", "chr", "start", "stop", "direction"])
    reg = _collapse_regions(df)
    rows = []
    for comp, sub in reg.groupby("comparison"):
        n_total = len(sub)
        n_hyper = int((sub["direction"] == "Hyper").sum())
        n_hypo = int((sub["direction"] == "Hypo").sum())
        rows.append(dict(comparison=comp, n_total=n_total, n_hyper=n_hyper, n_hypo=n_hypo,
                         prop_hyper=n_hyper / n_total, prop_hypo=n_hypo / n_total,
                         interpretation="descriptive DMR catalogue; no region-level inference"))
    return pd.DataFrame(rows)


def _directionality_figure(sub, lab_blue, lab_red, title, out_base,
                            col_hypo_override=None, col_hyper_override=None):
    ns.set_style()
    col_hyper = col_hyper_override if col_hyper_override is not None else ns.METHYL_COLORS["Hyper"]
    col_hypo = col_hypo_override if col_hypo_override is not None else ns.METHYL_COLORS["Hypo"]
    n = len(sub)
    fig = ns.new_fig_mm(92, 16 + 9 * n)
    axA = fig.add_subplot(111)
    y = np.arange(n)
    labels = list(sub["panel_label"])
    for i, (_, r) in enumerate(sub.iterrows()):
        axA.barh(i, r["prop_hyper"], color=col_hyper, height=0.6, zorder=2)
        axA.barh(i, -r["prop_hypo"], color=col_hypo, height=0.6, zorder=2)
        axA.text(r["prop_hyper"] / 2, i, f"{100*r['prop_hyper']:.0f}%\n(n={r['n_hyper']})",
                 ha="center", va="center", color="white", fontsize=ns.FS_SMALL, fontweight="bold")
        axA.text(-r["prop_hypo"] / 2, i, f"{100*r['prop_hypo']:.0f}%\n(n={r['n_hypo']})",
                 ha="center", va="center", color="white", fontsize=ns.FS_SMALL, fontweight="bold")
    axA.axvline(0, color="#444444", lw=0.6)
    axA.set_xlim(-1, 1); axA.set_ylim(-0.6, n - 0.4)
    axA.set_xticks([-1, -0.5, 0, 0.5, 1])
    axA.set_xticklabels(["100%", "50%", "0", "50%", "100%"], fontsize=ns.FS_TICK)
    axA.set_yticks(y); axA.set_yticklabels(labels, fontsize=ns.FS_AXIS)
    axA.set_xlabel("Fraction of DMR regions", fontsize=ns.FS_AXIS, fontweight="bold")
    axA.set_title(title, fontsize=ns.FS_TITLE, fontweight="bold")
    # Reserve explicit space for the full contrast labels and the legend.
    fig.subplots_adjust(left=0.27, right=0.98, top=0.84, bottom=0.32)
    # legend (hypo / hyper meaning)
    from matplotlib.patches import Patch
    axA.legend(handles=[Patch(fc=col_hypo, label=lab_blue), Patch(fc=col_hyper, label=lab_red)],
               loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1, fontsize=ns.FS_LEGEND,
               handlelength=1.0, borderpad=0.2)
    ns.save_panel(fig, out_base, vector=True, tiff=True)


def build_dmr_directionality():
    """DMR directionality, CSF1 effect (within-sex) and sex effect."""
    st = _directionality_stats().set_index("comparison")
    # Hyphenated group names, males listed first.
    csf1 = st.loc[["M_CSF_vs_M_C", "F_CSF_vs_F_C"]].copy()
    csf1["panel_label"] = ["M-CSF1 vs M-C", "F-CSF1 vs F-C"]
    # Use distinct palettes for methylation direction and sex effects.
    _directionality_figure(csf1, "Hypomethylated after CSF1", "Hypermethylated after CSF1",
                           "CSF1 effect within sex", OUT / "global_directionality" / "medip_dmr_directionality_csf1_effect",
                           col_hypo_override="#44AA99", col_hyper_override="#B2182B")
    sex = st.loc[["F_C_vs_M_C", "F_CSF_vs_M_CSF"]].copy()
    sex["panel_label"] = ["F-C vs M-C", "F-CSF1 vs M-CSF1"]
    # Sex-effect bars use the unified sex palette (female=blue, male=vermillion).
    _directionality_figure(sex, "Higher methylation in male", "Higher methylation in female",
                           "Sex effect (female vs male)", OUT / "global_directionality" / "medip_dmr_directionality_sex_effect",
                           col_hypo_override="#D55E00", col_hyper_override="#0072B2")
    print("[OK] DMR directionality (CSF1 + sex)")


# ----------------------------------------------------------------------------- 1E
GO_THEME_COLORS = {
    "Morphology": ns.THEME_COLORS["Morphology"],
    "Inflammation / neuroinflammation": ns.THEME_COLORS["Inflammation / neuroinflammation"],
    "Mechanical stimulus / pain": ns.THEME_COLORS["Mechanical stimulus / pain"],
    "Sex hormones / estrogen": ns.THEME_COLORS["Sex hormones / estrogen"],
}
GO_THEME_LEGEND = {
    "Morphology": "Morphology",
    "Inflammation / neuroinflammation": "Inflammation",
    "Mechanical stimulus / pain": "Mechanical stimulus / pain",
    "Sex hormones / estrogen": "Sex hormones / estrogen",
}


GO_LABEL_STOPWORDS = {
    "a", "an", "and", "at", "by", "cellular", "for", "from", "in",
    "into", "of", "on", "or", "pathway", "process", "regulation",
    "response", "the", "to", "via", "with",
}


def io_path(path: Path) -> Path:
    """Return a Windows extended-length path when the bundle is deeply nested."""
    if sys.platform == "win32" and not str(path).startswith("\\\\?\\"):
        return Path("\\\\?\\" + str(path.absolute()))
    return path


def _label_tokens(value):
    """Return normalized content words used for automatic label selection."""
    tokens = []
    for token in re.findall(r"[a-z]+", str(value).lower()):
        if token in GO_LABEL_STOPWORDS:
            continue
        if len(token) > 4 and token.endswith("s"):
            token = token[:-1]
        tokens.append(token)
    return set(tokens)


def _representative_theme_term(frame, theme, document_frequency, n_documents):
    """Select a representative term from a theme without term-specific rules.

    Terms sharing a rarer content word with the theme name are preferred;
    lexical coverage and enrichment strength provide deterministic tie-breaks.
    If no description shares a theme word, the strongest enriched term is used.
    """
    if frame.empty:
        return frame
    theme_tokens = _label_tokens(theme)
    scored = frame.copy()

    def score(description):
        description_tokens = _label_tokens(description)
        shared = theme_tokens & description_tokens
        if not shared:
            return (0.0, 0.0)
        specificity = max(
            math.log((n_documents + 1) / (document_frequency[token] + 1)) + 1
            for token in shared
        )
        coverage = len(shared) / max(len(description_tokens), 1)
        return (specificity, coverage)

    scores = scored["Description"].map(score)
    scored["_theme_specificity"] = scores.map(lambda value: value[0])
    scored["_theme_coverage"] = scores.map(lambda value: value[1])
    scored = scored.sort_values(
        ["_theme_specificity", "_theme_coverage", "neglog10_pvalue", "Count", "Description"],
        ascending=[False, False, False, False, True],
        kind="mergesort",
    )
    return scored.head(1).drop(columns=["_theme_specificity", "_theme_coverage"])


def _plot_individual_points(ax, frame, color, marker_size, alpha, zorder=2):
    """One artist per observation so PowerPoint exposes one editable point."""
    marker_size_points = np.sqrt(marker_size)
    for row in frame.itertuples(index=False):
        ax.plot(
            [row.rank_within_group],
            [row.neglog10_pvalue],
            marker="o",
            linestyle="none",
            markersize=marker_size_points,
            markerfacecolor=color,
            markeredgecolor="none",
            alpha=alpha,
            zorder=zorder,
        )


def build_go_rank_distribution():
    """GO relevance rank distribution; highlighted terms are labelled
    to the right of each panel with no overlap."""
    ns.set_style()
    df = pd.read_excel(io_path(DATA / "GO_relevance_rank_distribution_ws400.xlsx"))
    df["biological_theme"] = df["biological_theme"].replace(
        {"Steroid / sex hormones": "Sex hormones / estrogen"}
    )
    df = df[np.isfinite(df["rank_within_group"]) & np.isfinite(df["neglog10_pvalue"])].copy()
    description_tokens = df["Description"].map(_label_tokens)
    document_frequency = Counter(token for tokens in description_tokens for token in tokens)
    n_documents = len(df)
    # Keep only the two within-sex CSF1-effect contrasts (males first),
    # relabelled with the hyphenated group names.
    CSF1_CONTRASTS = ["GO M CSF1 vs M C", "GO F CSF1 vs F C"]
    DISPLAY_LABEL = {"GO M CSF1 vs M C": "M-CSF1 vs M-C", "GO F CSF1 vs F C": "F-CSF1 vs F-C"}
    df = df[df["comparison_label"].isin(CSF1_CONTRASTS)].copy()
    comparisons = [c for c in CSF1_CONTRASTS if c in set(df["comparison_label"])]
    directions = ["Hyper", "Hypo"]

    fig = ns.new_fig_mm(110, 92)
    gs = fig.add_gridspec(len(directions), len(comparisons), wspace=0.34, hspace=0.32,
                          left=0.085, right=0.995, top=0.90, bottom=0.14)
    plot_axes = []
    for r, direction in enumerate(directions):
        for c, comp in enumerate(comparisons):
            ax = fig.add_subplot(gs[r, c])
            plot_axes.append(ax)
            sub = df[(df["direction"] == direction) & (df["comparison_label"] == comp)].copy()
            _plot_individual_points(
                ax, sub, color="#C9C9C9", marker_size=2.0, alpha=0.5, zorder=2
            )
            for theme, color in GO_THEME_COLORS.items():
                th = sub[sub["biological_theme"] == theme]
                if len(th):
                    _plot_individual_points(
                        ax, th, color=color, marker_size=6.0, alpha=0.95, zorder=3
                    )
            # Select one representative term per highlighted biological theme.
            parts = [
                _representative_theme_term(
                    sub[sub["biological_theme"] == theme],
                    theme,
                    document_frequency,
                    n_documents,
                )
                for theme in GO_THEME_COLORS
            ]
            labels = pd.concat(parts, ignore_index=True).drop_duplicates(["Description"]).head(4)

            x_max = max(1, sub["rank_within_group"].max())
            y_top = max(1.0, sub["neglog10_pvalue"].max())
            ax.set_xlim(0, x_max * 3.0)
            ax.set_ylim(-0.05 * y_top, y_top * 1.10)
            texts = []
            # Edge-crossing minimisation: order the GO boxes top->bottom by the
            # vertical position of the point they annotate. All boxes share one
            # right-hand x-column, so a y-monotone box order is provably free
            # of connector-line crossings.
            top_labels = (labels.head(4)
                          .sort_values("neglog10_pvalue", ascending=False)
                          .reset_index(drop=True))
            # Position annotations in ordered lanes beside each panel.
            lanes = np.linspace(y_top * 0.93, y_top * 0.22, max(len(top_labels), 1))
            for lane, (_, row) in zip(lanes, top_labels.iterrows()):
                color = GO_THEME_COLORS.get(row["biological_theme"], "#666666")
                # Wrap annotation text within the right margin.
                txt = "\n".join(wrap(str(row["Description"]), 34))
                t = ax.annotate(txt, xy=(row["rank_within_group"], row["neglog10_pvalue"]),
                                xytext=(x_max * 1.16, lane), fontsize=5.5, color=color,
                                ha="left", va="center",
                                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=color, lw=0.5, alpha=0.98),
                                arrowprops=dict(arrowstyle="-", color="#666666", lw=0.4, ls="--"))
                texts.append(t)
            if r == 0:
                ax.set_title(DISPLAY_LABEL.get(comp, comp.replace("GO ", "")),
                             fontsize=ns.FS_AXIS, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{direction}\n-log$_{{10}}$(P value)",
                              fontsize=ns.FS_TICK, fontweight="bold")
            if r == len(directions) - 1:
                ax.set_xlabel("Rank", fontsize=ns.FS_TICK)
            ax.tick_params(labelsize=ns.FS_SMALL)
            ax.grid(True, color="#EEEEEE", lw=0.3)

    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=col,
                          markersize=4, label=GO_THEME_LEGEND[t])
               for t, col in GO_THEME_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=ns.FS_LEGEND, title="Biological theme", title_fontsize=ns.FS_LEGEND_TITLE,
               bbox_to_anchor=(0.5, 0.0))
    # Export the complete panel and each component axis.
    for index, target_ax in enumerate(plot_axes, start=1):
        for candidate_ax in plot_axes:
            candidate_ax.set_visible(candidate_ax is target_ax)
        for legend in fig.legends:
            legend.set_visible(False)
        part_base = OUT / f"medip_go_rank_distribution_part{index}"
        fig.savefig(
            part_base.with_suffix(".svg"),
            format="svg",
            bbox_inches="tight",
            pad_inches=0.08,
            facecolor="white",
            transparent=False,
        )
        fig.savefig(
            part_base.with_suffix(".pdf"),
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.08,
            facecolor="white",
            transparent=False,
        )
    for candidate_ax in plot_axes:
        candidate_ax.set_visible(True)
    for legend in fig.legends:
        legend.set_visible(True)
    ns.save_panel(fig, OUT / "medip_go_rank_distribution", vector=False, tiff=True)
    print("[OK] GO relevance rank distribution")


# ----------------------------------------------------------------------------- candidate-window map
def build_candidate_window_map():
    """Candidate-gene methylation map (significant 400-bp windows),
    region (x) by gene (y), colour = methylation direction, shape = comparison.
    Gene labels are enlarged and italicised."""
    ns.set_style()
    df = _load_sig()
    cg = pd.read_excel(io_path(CAND_XLSX), sheet_name="candidate_genes")[["symbol", "family"]]
    cg["sl"] = cg["symbol"].str.lower()
    irf_remove = {f"irf{i}" for i in range(2, 10)}
    cg = cg[~cg["sl"].isin(irf_remove)]
    fam = dict(zip(cg["sl"], cg["family"]))
    sym = dict(zip(cg["sl"], cg["symbol"]))

    d = df[["chr", "start", "stop", "external_transcript_gene_name", "region",
            "diffexpressed", "comparison", "edgeR.logFC", "edgeR.p.value"]].copy()
    d["sl"] = d["external_transcript_gene_name"].astype(str).str.lower()
    d = d[d["sl"].isin(set(cg["sl"]))]
    d = d[d["diffexpressed"].isin(["Hyper", "Hypo"])]
    d = _reannotate_candidate_windows(d, set(cg["sl"]))
    d["region"] = d["current_region"]
    d = d.drop_duplicates(["chr", "start", "stop", "sl", "diffexpressed", "comparison", "region"])
    # Display title-case italic gene symbols on the y-axis.
    def _sym_disp(s):
        g = sym[s]
        g = (g[0].upper() + g[1:].lower()) if g else g
        return "Tnfα" if g == "Tnf" else g
    d["gene_label"] = d["sl"].map(_sym_disp)
    comps = ["F_CSF_vs_F_C", "M_CSF_vs_M_C", "F_C_vs_M_C", "F_CSF_vs_M_CSF"]
    markers = {c: m for c, m in zip(comps, ["o", "^", "s", "D"])}
    regions = ["TSS", "EXON", "GENEBODY"]
    genes = sorted(d["gene_label"].unique(), reverse=True)
    rng = np.random.default_rng(7391)
    fig, ax = ns.fig_mm(92, max(60, 4 + 4.6 * len(genes)))
    xmap = {r: i for i, r in enumerate(regions)}
    ymap = {g: i for i, g in enumerate(genes)}
    for _, r in d.iterrows():
        if r["region"] not in xmap:
            continue
        x = xmap[r["region"]] + rng.uniform(-0.16, 0.16)
        y = ymap[r["gene_label"]] + rng.uniform(-0.16, 0.16)
        ax.scatter(x, y, s=14, color=ns.METHYL_COLORS.get(r["diffexpressed"], "#999"),
                   marker=markers.get(r["comparison"], "o"), alpha=0.8,
                   edgecolor="none", zorder=3)
    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(regions, fontsize=ns.FS_AXIS, fontweight="bold")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7.0)
    # italicising only the gene-symbol portion of the label is non-trivial, so
    # the whole label (symbol only, no suffix) is set to italic instead.
    for t in ax.get_yticklabels():
        t.set_fontstyle("italic")
    ax.set_xlim(-0.5, len(regions) - 0.5)
    ax.set_ylim(-0.6, len(genes) - 0.4)
    ax.tick_params(length=0)
    ax.grid(axis="y", color="#EEEEEE", lw=0.3)
    ax.set_title("Candidate-gene methylation (significant windows)", fontsize=ns.FS_TITLE, fontweight="bold")
    dir_h = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=ns.METHYL_COLORS[k],
                        markersize=4, label=k) for k in ["Hyper", "Hypo"]]
    comp_h = [plt.Line2D([0], [0], marker=markers[c], color="w", markerfacecolor="#777",
                         markersize=4, label=COMPARISON_LABELS[c]) for c in comps]
    leg1 = ax.legend(handles=dir_h, title="Methylation", fontsize=ns.FS_LEGEND,
                     title_fontsize=ns.FS_LEGEND_TITLE, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.add_artist(leg1)
    ax.legend(handles=comp_h, title="Comparison", fontsize=ns.FS_LEGEND,
              title_fontsize=ns.FS_LEGEND_TITLE, loc="upper left", bbox_to_anchor=(1.02, 0.6))
    ns.save_panel(fig, OUT / "medip_candidate_window_map", vector=True, tiff=True)
    print(f"[OK] candidate-gene methylation map ({len(genes)} genes)")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    if a.only in ("", "heatmap"):
        build_candidate_gene_heatmap()
    if a.only in ("", "directionality"):
        build_dmr_directionality()
    if a.only in ("", "go"):
        build_go_rank_distribution()
    if a.only in ("", "window_map"):
        build_candidate_window_map()
