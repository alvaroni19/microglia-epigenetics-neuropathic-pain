"""Validate the regenerated MeDIP tables before figure construction."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).absolute().parents[2]
TABLES = ROOT / "data" / "MeDIP" / "workspace_and_tables"
CSF1R_START = 56_414_488
CSF1R_STOP = 56_458_300
CSF1R_MARGIN = 1_500
CONTEXT_START = 56_411_488
CONTEXT_STOP = 56_461_300


def io_path(path: Path) -> Path:
    """Return a Windows extended-length path when the bundle is deeply nested."""
    if sys.platform == "win32" and not str(path).startswith("\\\\?\\"):
        return Path("\\\\?\\" + str(path.absolute()))
    return path


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise RuntimeError(f"{label} is missing columns: {', '.join(missing)}")


def validate_csf1r_coordinates(frame: pd.DataFrame, label: str) -> None:
    require_columns(frame, {"gene", "chr", "start", "stop"}, label)
    csf1r = frame.loc[frame["gene"].astype(str).str.lower() == "csf1r"].copy()
    if csf1r.empty:
        raise RuntimeError(f"{label} contains no Csf1r rows")
    coordinates_are_rn6 = (
        (csf1r["chr"].astype(str) == "chr18")
        & (pd.to_numeric(csf1r["start"], errors="coerce") <= CSF1R_STOP + CSF1R_MARGIN)
        & (pd.to_numeric(csf1r["stop"], errors="coerce") >= CSF1R_START - CSF1R_MARGIN)
    )
    if not coordinates_are_rn6.all():
        raise RuntimeError(f"{label} contains Csf1r coordinates outside the rn6 locus")


def validate_fasta() -> set[tuple[str, int, int]]:
    fasta = TABLES / "csf1r_dmr_windows_rn6.fasta"
    lines = io_path(fasta).read_text(encoding="utf-8").splitlines()
    if len(lines) == 0 or len(lines) % 2 != 0:
        raise RuntimeError("Csf1r FASTA must contain header-sequence pairs")
    windows: set[tuple[str, int, int]] = set()
    for header, sequence in zip(lines[0::2], lines[1::2]):
        match = re.fullmatch(r">(chr[^:]+):(\d+)-(\d+)", header)
        if match is None:
            raise RuntimeError(f"Invalid FASTA header: {header}")
        chromosome, start_text, stop_text = match.groups()
        start, stop = int(start_text), int(stop_text)
        if chromosome != "chr18" or start > CSF1R_STOP + CSF1R_MARGIN or stop < CSF1R_START - CSF1R_MARGIN:
            raise RuntimeError(f"FASTA window is outside the rn6 Csf1r locus: {header}")
        if len(sequence) != stop - start + 1:
            raise RuntimeError(f"FASTA sequence length does not match coordinates: {header}")
        windows.add((chromosome, start, stop))
    return windows


def validate_context_fasta() -> None:
    lines = io_path(TABLES / "csf1r_context_rn6.fasta").read_text(
        encoding="utf-8"
    ).splitlines()
    if not lines or lines[0] != f">rn6_chr18:{CONTEXT_START}-{CONTEXT_STOP}":
        raise RuntimeError("Invalid rn6 Csf1r context FASTA header")
    sequence = "".join(lines[1:]).upper()
    if len(sequence) != CONTEXT_STOP - CONTEXT_START + 1:
        raise RuntimeError("Invalid rn6 Csf1r context FASTA length")
    if set(sequence) - set("ACGTN"):
        raise RuntimeError("Invalid characters in rn6 Csf1r context FASTA")


def main() -> None:
    genes = pd.read_excel(io_path(TABLES / "GENES_TOTAL_sig_contig_colap_ws400.xlsx"))
    require_columns(
        genes,
        {
            "chr", "start", "stop", "edgeR.logFC", "edgeR.p.value",
            "diffexpressed", "external_transcript_gene_name", "comparison", "region",
        },
        "significant-gene table",
    )
    if genes.empty or not set(genes["diffexpressed"].dropna().astype(str)) <= {"Hyper", "Hypo"}:
        raise RuntimeError("Significant-gene directions are empty or non-English")

    go = pd.read_excel(io_path(TABLES / "GO_relevance_rank_distribution_ws400.xlsx"))
    require_columns(
        go,
        {"comparison", "direction", "Description", "p_value", "neglog10_pvalue"},
        "GO rank table",
    )
    if go.empty or not set(go["direction"].dropna().astype(str)) <= {"Hyper", "Hypo"}:
        raise RuntimeError("GO directions are empty or invalid")
    p_values = pd.to_numeric(go["p_value"], errors="coerce")
    if p_values.isna().any() or not ((p_values >= 0) & (p_values < 0.05)).all():
        raise RuntimeError("GO p_value contains missing values or values outside [0, 0.05)")

    candidate = pd.read_excel(
        io_path(TABLES / "DMR_candidate_genes_ERE_scan_ws400.xlsx"), sheet_name="ERE_scan"
    )
    validate_csf1r_coordinates(candidate, "candidate-gene ERE table")
    fasta_windows = validate_fasta()
    validate_context_fasta()

    rms = pd.read_csv(io_path(TABLES / "significant_gene_rms_matrix_ws400.csv"))
    if rms.empty or rms.shape[1] != 13 or rms.columns[0] != "gene":
        raise RuntimeError("Significant-gene RMS matrix must contain gene plus 12 samples")
    metadata = pd.read_csv(io_path(TABLES / "medip_sample_metadata.csv"))
    if len(metadata) != 12 or set(metadata["sample_label"]) != set(rms.columns[1:]):
        raise RuntimeError("MeDIP sample metadata does not match the RMS matrix")

    motif_hits_path = ROOT / "generated_outputs" / "MeDIP" / "csf1r_dmr_window_esr_pwm_hits_ws400.csv"
    motif_hits = pd.read_csv(io_path(motif_hits_path))
    require_columns(motif_hits, {"chr", "start", "stop"}, "ESR PWM table")
    for row in motif_hits.itertuples(index=False):
        contained = any(
            str(row.chr) == chromosome
            and int(row.start) >= window_start
            and int(row.stop) <= window_stop
            for chromosome, window_start, window_stop in fasta_windows
        )
        if not contained:
            raise RuntimeError("ESR PWM output contains a hit outside the rn6 Csf1r FASTA")

    print(
        f"[OK] MeDIP outputs validated: {len(genes)} DMR rows, "
        f"{len(go)} GO rows, {len(fasta_windows)} Csf1r windows"
    )


if __name__ == "__main__":
    main()
