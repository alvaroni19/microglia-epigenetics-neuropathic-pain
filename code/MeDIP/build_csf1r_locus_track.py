# -*- coding: utf-8 -*-
"""Build the rn6 Csf1r gene, CpG, DMR and ERE tracks."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402
import map_csf1r_esr_motifs_integrated as mc  # noqa: E402 (read_gene_model, paths)

ns.set_style()

ERE_SCAN = (
    ROOT
    / "data"
    / "MeDIP"
    / "workspace_and_tables"
    / "DMR_candidate_genes_ERE_scan_ws400.xlsx"
)
GENERATED_HITS = ROOT / "generated_outputs" / "MeDIP" / "csf1r_dmr_window_esr_pwm_hits_ws400.csv"
FROZEN_HITS = ROOT / "data" / "MeDIP" / "csf1r_dmr_window_esr_pwm_hits_ws400.csv"
CONTEXT_FASTA = (
    ROOT / "data" / "MeDIP" / "workspace_and_tables" / "csf1r_context_rn6.fasta"
)
OUT = ROOT / "generated_outputs" / "MeDIP"

HYPER, HYPO = "#B2182B", "#2166AC"     # red up / blue down (methylation track)
ERE_FILL, ERE_EDGE = "#F0E442", "#7A6E00"
GENE_COL, EXON_COL, TSS_COL = "#117733", "#6A3D9A", "#E69F00"

PAD_GENE = 3000


def io_path(path: Path) -> Path:
    """Return a Windows extended-length path when the bundle is deeply nested."""
    if sys.platform == "win32" and not str(path).startswith("\\\\?\\"):
        return Path("\\\\?\\" + str(path.absolute()))
    return path


def _is_hyper(direction: str) -> bool:
    d = str(direction).upper()
    return d.startswith("HYPER")


def _read_context_sequence() -> str:
    if not io_path(CONTEXT_FASTA).exists():
        raise FileNotFoundError(
            f"Missing {CONTEXT_FASTA}. Re-run the MEDIPS notebook to export "
            "the rn6 Csf1r context sequence."
        )
    lines = io_path(CONTEXT_FASTA).read_text(encoding="ascii").splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


def _cpg_density(sequence: str, genomic_start: int, window: int = 400) -> pd.DataFrame:
    rows = []
    sequence = sequence.upper()
    for offset in range(0, len(sequence), window):
        fragment = sequence[offset:offset + window]
        if not fragment:
            continue
        cpg_count = sum(fragment[index:index + 2] == "CG" for index in range(len(fragment) - 1))
        rows.append(
            {
                "start": genomic_start + offset,
                "stop": genomic_start + offset + len(fragment) - 1,
                "cpg_per_kb": 1000.0 * cpg_count / len(fragment),
            }
        )
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hits_file = GENERATED_HITS if io_path(GENERATED_HITS).exists() else FROZEN_HITS
    if not io_path(hits_file).exists():
        raise FileNotFoundError(
            "Missing Csf1r motif-hit table. Run map_csf1r_esr_motifs_integrated.py "
            f"or provide {FROZEN_HITS}."
        )
    gene_start, gene_stop, strand, exons = mc.read_gene_model()
    x0, x1 = int(gene_start) - PAD_GENE, int(gene_stop) + PAD_GENE

    scan = pd.read_excel(io_path(ERE_SCAN), sheet_name="ERE_scan")
    # The ERE scan table stores within-sex CSF1 contrasts control-first
    # (F_C_vs_F_CSF) and sex contrasts female-first (F_C_vs_M_C). MeDIP panels
    # write CSF1 contrasts CSF1 first and sex contrasts male first, so both
    # kinds are re-expressed (log2FC sign and direction inverted).
    csf1_rev = scan["comparison"].astype(str).str.fullmatch(r"([FM])_C_vs_\1_CSF")
    scan.loc[csf1_rev, "edgeR.logFC"] = -scan.loc[csf1_rev, "edgeR.logFC"]
    scan.loc[csf1_rev, "methylation_direction"] = scan.loc[csf1_rev, "methylation_direction"].map(
        lambda d: "HYPO" if str(d).upper().startswith("HYPER") else "HYPER")
    scan.loc[csf1_rev, "comparison"] = scan.loc[csf1_rev, "comparison"].str.replace(
        r"^([FM])_C_vs_\1_CSF$", r"\1_CSF_vs_\1_C", regex=True)
    reverse = scan["comparison"].astype(str).str.fullmatch(r"F_(C|CSF)_vs_M_\1")
    scan.loc[reverse, "edgeR.logFC"] = -scan.loc[reverse, "edgeR.logFC"]
    scan.loc[reverse, "methylation_direction"] = scan.loc[reverse, "methylation_direction"].map(
        lambda d: "HYPO" if str(d).upper().startswith("HYPER") else "HYPER")
    scan.loc[reverse, "comparison"] = scan.loc[reverse, "comparison"].str.replace(
        r"^F_(C|CSF)_vs_M_\1$", r"M_\1_vs_F_\1", regex=True)
    # One bar per genomic window: collapse duplicate (start, stop) across the
    # pairwise comparisons (direction/logFC are consistent per window here).
    scan = (scan.sort_values("edgeR.logFC", key=lambda s: s.abs(), ascending=False)
                .drop_duplicates(["start", "stop"], keep="first")
                .sort_values("start")
                .reset_index(drop=True))
    hits = pd.read_csv(io_path(hits_file))
    scan = scan[(scan["chr"].astype(str).eq("chr18"))
                & (scan["stop"] >= x0)
                & (scan["start"] <= x1)].copy()
    if scan.empty:
        raise RuntimeError(
            "No rn6 MeDIP windows overlap the rn6 Csf1r context. Regenerate the "
            "MeDIP gene annotation from the bundled rn6 GTF before building the Csf1r locus track."
        )

    hits = hits[(hits["stop"] >= x0) & (hits["start"] <= x1)].copy()
    sequence = _read_context_sequence()
    density = _cpg_density(sequence, x0)

    fig = ns.new_fig_mm(180, 43)
    gs = fig.add_gridspec(
        3, 1, height_ratios=[0.85, 1.7, 0.9], hspace=0.04,
        left=0.075, right=0.99, top=0.96, bottom=0.25,
    )
    ax_gene = fig.add_subplot(gs[0, 0])
    ax_cpg = fig.add_subplot(gs[1, 0], sharex=ax_gene)
    ax_dmr = fig.add_subplot(gs[2, 0], sharex=ax_gene)
    ax_gene.set_xlim(x0, x1)

    ax_gene.hlines(0, gene_start, gene_stop, color=GENE_COL, lw=1.5, zorder=2)
    for _, exon in exons.iterrows():
        ax_gene.add_patch(
            Rectangle(
                (exon["start"], -0.32), max(int(exon["stop"] - exon["start"]), 60), 0.64,
                facecolor=EXON_COL, edgecolor="#3D2350", lw=0.25, zorder=3,
            )
        )
    tss = gene_stop if strand == -1 else gene_start
    ax_gene.axvline(tss, ymin=0.1, ymax=0.92, color=TSS_COL, lw=1.0)
    ax_gene.text(tss, 0.70, "TSS", ha="center", va="bottom", fontsize=ns.FS_SMALL,
                 color=TSS_COL, fontweight="bold")
    ax_gene.text(gene_start + 1700, -0.67, "GENE BODY", ha="left", va="center",
                 fontsize=ns.FS_SMALL, color=GENE_COL, fontweight="bold")
    ax_gene.text(gene_start + 7200, -0.67, "EXONS", ha="left", va="center",
                 fontsize=ns.FS_SMALL, color=EXON_COL, fontweight="bold")
    ax_gene.set_ylim(-0.95, 0.95)
    ax_gene.axis("off")

    ax_cpg.bar(
        density["start"], density["cpg_per_kb"],
        width=density["stop"] - density["start"] + 1,
        align="edge", color="#BDBDBD", edgecolor="none",
    )
    ax_cpg.set_ylabel("CpG / kb", fontsize=ns.FS_AXIS, fontweight="bold")
    ax_cpg.tick_params(labelsize=ns.FS_TICK, length=2)
    ax_cpg.spines[["top", "right"]].set_visible(False)
    ax_cpg.tick_params(axis="x", labelbottom=False)

    ax_dmr.set_ylim(-0.38, 1.05)
    ax_dmr.axhline(0.30, color="#777777", lw=0.45)
    ordered_scan = scan.sort_values("start").reset_index(drop=True)
    window_midpoints = ((ordered_scan["start"] + ordered_scan["stop"]) / 2).to_numpy()
    for index, row in ordered_scan.iterrows():
        window_number = index + 1
        hyper = _is_hyper(row["methylation_direction"])
        color = HYPER if hyper else HYPO
        midpoint = (row["start"] + row["stop"]) / 2
        label_midpoint = midpoint
        if index + 1 < len(window_midpoints) and window_midpoints[index + 1] - midpoint < 1200:
            label_midpoint -= 520
        elif index > 0 and midpoint - window_midpoints[index - 1] < 1200:
            label_midpoint += 520
        ax_dmr.add_patch(
            Rectangle(
                (row["start"], 0.20), max(int(row["stop"] - row["start"]), 90), 0.20,
                facecolor="white", edgecolor=color, lw=0.7, zorder=3,
            )
        )
        ax_dmr.text(
            label_midpoint, 0.02, f"W{window_number}",
            ha="center", va="top", fontsize=ns.FS_SMALL,
            color=color, fontweight="bold",
        )
    for _, hit in hits.iterrows():
        midpoint = (int(hit["start"]) + int(hit["stop"])) / 2.0
        ax_dmr.add_patch(
            Rectangle((midpoint - 55, 0.55), 110, 0.18, facecolor="black",
                      edgecolor="black", lw=0.25, zorder=5)
        )
    ax_dmr.text(x0 + 0.50 * (x1 - x0), 0.64, "ERE", ha="center", va="center",
                fontsize=ns.FS_SMALL, fontweight="bold")
    ax_dmr.set_yticks([])
    ax_dmr.spines[["top", "right", "left"]].set_visible(False)
    ax_dmr.ticklabel_format(style="plain", axis="x", useOffset=False)
    ax_dmr.tick_params(axis="x", labelsize=ns.FS_TICK, length=2)

    legend_handles = [
        Patch(facecolor="white", edgecolor=HYPO, label="Hypomethylated DMR"),
        Patch(facecolor="white", edgecolor=HYPER, label="Hypermethylated DMR"),
    ]
    fig.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.47, 0.09),
        ncol=2, frameon=False, fontsize=ns.FS_LEGEND,
        title="MeDIP 400-bp windows", title_fontsize=ns.FS_LEGEND_TITLE,
        handlelength=0.8, handletextpad=0.35, columnspacing=0.8,
    )

    ns.save_panel(fig, OUT / "medip_csf1r_directional_track", vector=True, tiff=True)
    print(f"[OK] medip_csf1r_directional_track ({len(scan)} rn6 windows, "
          f"{len(hits)} ERE hits)")


if __name__ == "__main__":
    main()
