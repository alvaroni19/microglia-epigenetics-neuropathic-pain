"""Scan rn6 Csf1r DMR windows for ESR1 and ESR2 PWM matches."""

from __future__ import annotations

import gzip
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
MEDIP_OUTPUT = ROOT / "generated_outputs" / "MeDIP"
MEDIP_DATA = ROOT / "data" / "MeDIP"
DMR_XLSX = MEDIP_DATA / "workspace_and_tables" / "DMR_candidate_genes_ERE_scan_ws400.xlsx"
DMR_FASTA = MEDIP_DATA / "workspace_and_tables" / "csf1r_dmr_windows_rn6.fasta"
ESR1_MEME = MEDIP_DATA / "motifs" / "MA0112.3.meme"
ESR2_MEME = MEDIP_DATA / "motifs" / "MA0258.2.meme"
DMR_HITS_CSV = MEDIP_OUTPUT / "csf1r_dmr_window_esr_pwm_hits_ws400.csv"
GTF = MEDIP_DATA / "annotations" / "rn6" / "Rattus_norvegicus.Rnor_6.0.104.gtf.gz"

GENE_ID = "ENSRNOG00000018414"
BASES = "ACGT"
RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass
class FastaRecord:
    header: str
    sequence: str
    chromosome: str | None
    start: int | None
    stop: int | None


@dataclass
class Motif:
    motif_id: str
    name: str
    pwm: np.ndarray


def make_record(header: str, sequence: str) -> FastaRecord:
    match = re.search(r"(chr[\w]+):([0-9,]+)-([0-9,]+)", header)
    if match is None:
        return FastaRecord(header, sequence, None, None, None)
    return FastaRecord(
        header,
        sequence,
        match.group(1),
        int(match.group(2).replace(",", "")),
        int(match.group(3).replace(",", "")),
    )


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(make_record(header, "".join(sequence).upper()))
                header, sequence = line[1:], []
            else:
                sequence.append(line)
    if header is not None:
        records.append(make_record(header, "".join(sequence).upper()))
    return records


def read_meme(path: Path) -> Motif:
    motif_id = motif_name = path.stem
    matrix: list[list[float]] = []
    reading_matrix = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("MOTIF"):
                parts = line.split()
                motif_id = parts[1]
                motif_name = parts[2] if len(parts) > 2 else motif_id
            elif line.startswith("letter-probability matrix"):
                reading_matrix = True
            elif reading_matrix:
                values = line.split()
                if len(values) != 4:
                    break
                matrix.append([float(value) for value in values])
    if not matrix:
        raise RuntimeError(f"No PWM matrix found in {path}")
    pwm = np.clip(np.asarray(matrix, dtype=float), 1e-6, 1.0)
    pwm /= pwm.sum(axis=1, keepdims=True)
    return Motif(motif_id, motif_name, pwm)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1].upper()


def score_sequence(sequence: str, motif: Motif) -> float:
    indices = [BASES.find(base) for base in sequence.upper()]
    if any(index < 0 for index in indices):
        return float("-inf")
    return float(
        sum(math.log2(motif.pwm[position, index] / 0.25)
            for position, index in enumerate(indices))
    )


def scan_record(record: FastaRecord, motif: Motif, threshold: float) -> pd.DataFrame:
    width = motif.pwm.shape[0]
    log_odds = np.log2(motif.pwm / 0.25)
    maximum = float(log_odds.max(axis=1).sum())
    minimum = float(log_odds.min(axis=1).sum())
    rows = []
    for offset in range(len(record.sequence) - width + 1):
        forward = record.sequence[offset:offset + width]
        for strand, candidate in (("+", forward), ("-", reverse_complement(forward))):
            score = score_sequence(candidate, motif)
            relative_score = (score - minimum) / (maximum - minimum)
            if relative_score >= threshold:
                start = int(record.start) + offset
                rows.append({
                    "chr": record.chromosome,
                    "start": start,
                    "stop": start + width - 1,
                    "motif_id": motif.motif_id,
                    "motif_name": motif.name,
                    "strand": strand,
                    "sequence": candidate,
                    "log_odds_score": score,
                    "relative_score": relative_score,
                })
    return pd.DataFrame(rows)


def read_gene_model() -> tuple[int, int, int, pd.DataFrame]:
    gene_start = gene_stop = None
    strand_value = 1
    exons: list[tuple[int, int]] = []
    with gzip.open(GTF, "rt", encoding="utf-8") as handle:
        for line in handle:
            if GENE_ID not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            start, stop = int(fields[3]), int(fields[4])
            if fields[2] == "gene":
                gene_start, gene_stop = start, stop
                strand_value = -1 if fields[6] == "-" else 1
            elif fields[2] == "exon":
                exons.append((start, stop))
    if gene_start is None or gene_stop is None:
        raise RuntimeError(f"Could not locate {GENE_ID} in {GTF}")
    merged: list[list[int]] = []
    for start, stop in sorted(exons):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, stop])
        else:
            merged[-1][1] = max(merged[-1][1], stop)
    return gene_start, gene_stop, strand_value, pd.DataFrame(merged, columns=["start", "stop"])


def read_csf1r_windows() -> pd.DataFrame:
    windows = pd.read_excel(DMR_XLSX, sheet_name="ERE_scan")
    windows = windows[windows["gene"].astype(str).str.lower().eq("csf1r")].copy()
    windows["region_biological"] = windows["region"].replace({"GENE": "GENEBODY"})
    windows = windows.sort_values(
        ["start", "stop", "comparison", "region_biological"]
    ).drop_duplicates(["start", "stop", "comparison"], keep="first")
    windows = windows.reset_index(drop=True)
    windows["window_id"] = [f"W{index}" for index in range(1, len(windows) + 1)]
    if windows.empty:
        raise RuntimeError("The rn6 DMR table contains no Csf1r-associated windows.")
    return windows


def main() -> None:
    threshold = 0.92
    windows = read_csf1r_windows()
    fasta_by_coordinate = {
        (record.chromosome, record.start, record.stop): record
        for record in read_fasta(DMR_FASTA)
    }
    motifs = [read_meme(ESR1_MEME), read_meme(ESR2_MEME)]
    hit_tables = []
    for row in windows.itertuples(index=False):
        key = (str(row.chr), int(row.start), int(row.stop))
        record = fasta_by_coordinate.get(key)
        if record is None:
            raise RuntimeError(f"Missing rn6 FASTA sequence for DMR window {key}")
        for motif in motifs:
            hits = scan_record(record, motif, threshold)
            if hits.empty:
                continue
            hits["window_id"] = row.window_id
            hits["comparison"] = row.comparison
            hits["region_biological"] = row.region_biological
            hits["methylation_direction"] = row.methylation_direction
            hit_tables.append(hits)

    columns = [
        "window_id", "comparison", "region_biological", "methylation_direction",
        "chr", "start", "stop", "motif_id", "motif_name", "strand",
        "sequence", "log_odds_score", "relative_score",
    ]
    hits = pd.concat(hit_tables, ignore_index=True) if hit_tables else pd.DataFrame(columns=columns)
    hits = hits.reindex(columns=columns)
    MEDIP_OUTPUT.mkdir(parents=True, exist_ok=True)
    hits.to_csv(DMR_HITS_CSV, index=False)
    print(f"[OK] {len(hits)} PWM matches written to {DMR_HITS_CSV}")


if __name__ == "__main__":
    main()
