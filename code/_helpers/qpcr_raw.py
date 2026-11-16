"""Canonical RT-qPCR loaders for the manuscript.

In vivo values are reconstructed from the animal-level Delta-Ct sheets in
``data/PCR_260514/invivo_qpcr_delta_ct.xlsx``.  For every gene, male controls are the
single cross-sex calibrator:

    DeltaDeltaCt = DeltaCt(sample) - mean(DeltaCt(male controls))
    relative_expression = 2 ** (-DeltaDeltaCt)

Using one calibrator per gene preserves between-sex and between-phase
differences. Multiplying every observation for a gene by another positive
constant would not change t statistics, ANOVA tests, Spearman correlations or
z-scored multivariate analyses, but this convention also makes the plotted
fold changes reproducible.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


GENE_SHEETS = {
    "Dnmt1": "Dnmt1 new",
    "Dnmt3a": "Dnmt3a new",
    "Dnmt3b": "Dnmt3b new",
    "Gadd45a": "Gadd45a new",
    "Gadd45b": "Gadd45b new",
    "Gadd45g": "Gadd45g new",
    "Il10": "IL10 new",
    "Bdnf": "BDNF new",
    "Tlr4": "TLR4 new",
    "Slc12a5": "KCC2 new",
    "Tnf": "TNFa new",
    "Csf1r": "CSF1R new",
    "Esr1": "ERa new 2",
    "Esr2": "ERb new",
    "Gper1": "GPER new",
    "Cyp19a1": "ARO new",
}

SCRIPT_TO_CANONICAL = {
    "Dnmt1": "Dnmt1",
    "Dnmt3a": "Dnmt3a",
    "Dnmt3b": "Dnmt3b",
    "Gadd45a": "Gadd45a",
    "Gadd45b": "Gadd45b",
    "Gadd45g": "Gadd45g",
    "IL10": "Il10",
    "BDNF": "Bdnf",
    "TLR4": "Tlr4",
    "KCC2": "Slc12a5",
    "TNFa": "Tnf",
    "CSF1R": "Csf1r",
    "ERa": "Esr1",
    "ERb": "Esr2",
    "GPER": "Gper1",
    "ARO": "Cyp19a1",
}

# The final MFA/correlation analysis deliberately excludes prkcg, p2rx4, AR,
# P4R and hsd17b3; they are outside the final marker set for this manuscript.
MFA_GENES = [
    "Slc12a5",
    "Bdnf",
    "Csf1r",
    "Tlr4",
    "Esr1",
    "Esr2",
    "Gper1",
    "Cyp19a1",
    "Dnmt1",
    "Dnmt3a",
    "Dnmt3b",
    "Gadd45a",
    "Gadd45b",
    "Gadd45g",
]

ESTROGEN_SCORE_GENES = ["Esr1", "Esr2", "Gper1", "Cyp19a1"]

# Measurements excluded for technical reasons before any analysis. These are
# the same exclusions applied in the original GraphPad Prism analysis of the
# in vivo RT-qPCR; the corresponding Delta-Ct values are set to missing for
# the listed gene only (all other genes of these animals are retained).
TECHNICAL_EXCLUSIONS = {
    ("Cyp19a1", "RH36"): "proestrus control; technical exclusion (as in the GraphPad analysis)",
    ("Cyp19a1", "RH37"): "diestrus CCI-IoN; technical exclusion (as in the GraphPad analysis)",
}


def _default_root() -> Path:
    # Keep the mapped-drive path on Windows. Resolving it to the longer UNC
    # target can exceed the legacy path-length limit for generated filenames.
    return Path(__file__).absolute().parents[2]


def _normalize_group(value: object) -> tuple[str, str]:
    text = re.sub(r"\s+", " ", str(value).strip().lower())
    phase_map = {
        "male": "Male",
        "proestrus": "Proestrus",
        "diestrus": "Diestrus",
        "estrus": "Estrus",
    }
    condition = "Control" if "control" in text else "CCI" if "cci" in text else ""
    phase = next((english for spanish, english in phase_map.items() if spanish in text), "")
    if not phase or not condition:
        raise ValueError(f"Unrecognized qPCR group label: {value!r}")
    return phase, condition


def _animal_id(sample_id: int, phase: str) -> str:
    prefix = "RM" if phase == "Male" else "RH"
    return f"{prefix}{sample_id}"


def load_invivo_long(
    workbook: Path | None = None,
) -> pd.DataFrame:
    """Return one row per animal and gene with Delta-Ct and relative expression."""
    root = _default_root()
    workbook = workbook or root / "data" / "PCR_260514" / "invivo_qpcr_delta_ct.xlsx"
    frames: list[pd.DataFrame] = []

    available = set(pd.ExcelFile(workbook).sheet_names)
    missing_sheets = sorted(set(GENE_SHEETS.values()) - available)
    if missing_sheets:
        raise ValueError(f"Missing final Delta-Ct sheets: {missing_sheets}")

    for gene, sheet in GENE_SHEETS.items():
        frame = pd.read_excel(workbook, sheet_name=sheet, header=2).iloc[:, :3].copy()
        frame.columns = ["source_group", "id", "delta_ct"]
        frame = frame.dropna(subset=["id"]).copy()
        frame["id"] = pd.to_numeric(frame["id"], errors="raise").astype(int)
        frame["delta_ct"] = pd.to_numeric(frame["delta_ct"], errors="coerce")
        parsed = frame["source_group"].map(_normalize_group)
        frame["phase"] = parsed.map(lambda item: item[0])
        frame["condition"] = parsed.map(lambda item: item[1])
        frame["animal_id"] = [
            _animal_id(sample_id, phase)
            for sample_id, phase in zip(frame["id"], frame["phase"])
        ]
        frame["gene"] = gene
        for (excluded_gene, excluded_animal), _reason in TECHNICAL_EXCLUSIONS.items():
            if excluded_gene == gene:
                frame.loc[frame["animal_id"] == excluded_animal, "delta_ct"] = np.nan
        frame["source_sheet"] = sheet

        duplicated = frame.duplicated(["animal_id"], keep=False)
        if duplicated.any():
            ids = sorted(frame.loc[duplicated, "animal_id"].unique())
            raise ValueError(f"Duplicate animals in {sheet}: {ids}")

        calibrator = frame.loc[
            (frame["phase"] == "Male")
            & (frame["condition"] == "Control"),
            "delta_ct",
        ].dropna()
        if len(calibrator) < 2:
            raise ValueError(
                f"{sheet}: fewer than two non-missing male-control calibrators"
            )
        calibrator_mean = float(calibrator.mean())
        frame["calibrator_mean_delta_ct"] = calibrator_mean
        frame["delta_delta_ct"] = frame["delta_ct"] - calibrator_mean
        frame["relative_expression"] = np.power(2.0, -frame["delta_delta_ct"])
        frames.append(frame)

    output = pd.concat(frames, ignore_index=True)
    return output[
        [
            "animal_id",
            "id",
            "phase",
            "condition",
            "gene",
            "delta_ct",
            "calibrator_mean_delta_ct",
            "delta_delta_ct",
            "relative_expression",
            "source_sheet",
        ]
    ]


def _canonical_behavior_id(value: object) -> str:
    text = re.sub(r"\s+", "", str(value).upper())
    male = re.fullmatch(r"R(\d+)M", text)
    if male:
        return f"RM{int(male.group(1))}"
    female = re.fullmatch(r"R(?:H)?(\d+)", text)
    if female:
        return f"RH{int(female.group(1))}"
    raise ValueError(f"Unrecognized behavioural animal identifier: {value!r}")


def load_behavior_016(
    workbook: Path | None = None,
) -> pd.DataFrame:
    """Return baseline, test and test-minus-baseline 0.16-g values by animal."""
    root = _default_root()
    workbook = (
        workbook
        or root / "data" / "Allodynia_0_16" / "allodynia_von_frey_scores.xlsx"
    )
    frame = pd.read_excel(workbook, sheet_name="allodynia_plot")
    required = {"n", "measurement", "Surgery", "pain"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing behavioural columns: {sorted(missing)}")
    frame = frame[list(required)].copy()
    frame["animal_id"] = frame["n"].map(_canonical_behavior_id)
    frame["measurement"] = frame["measurement"].astype(str).str.strip().str.lower()
    frame["condition_behavior"] = (
        frame["Surgery"].astype(str).str.strip().str.capitalize()
    )
    frame["pain_016"] = pd.to_numeric(frame["pain"], errors="raise")
    frame = frame[frame["measurement"].isin(["bl", "test"])].copy()

    duplicated = frame.duplicated(["animal_id", "measurement"], keep=False)
    if duplicated.any():
        keys = (
            frame.loc[duplicated, ["animal_id", "measurement"]]
            .drop_duplicates()
            .to_dict("records")
        )
        raise ValueError(f"Duplicate behavioural measurements: {keys}")

    values = frame.pivot(
        index="animal_id", columns="measurement", values="pain_016"
    )
    conditions = frame.groupby("animal_id")["condition_behavior"].agg(
        lambda values_: values_.iloc[0]
        if values_.nunique() == 1
        else (_ for _ in ()).throw(ValueError("Inconsistent behaviour condition"))
    )
    output = values.join(conditions).reset_index()
    missing_measurements = output[["bl", "test"]].isna().any(axis=1)
    if missing_measurements.any():
        ids = output.loc[missing_measurements, "animal_id"].tolist()
        raise ValueError(f"Animals missing baseline or test 0.16-g values: {ids}")
    return output.rename(
        columns={
            "bl": "pain_016_baseline",
            "test": "pain_016_test",
        }
    ).assign(
        pain_016_delta=lambda data: (
            data["pain_016_test"] - data["pain_016_baseline"]
        )
    )


def load_invivo_wide() -> pd.DataFrame:
    """Return the complete 34-animal qPCR table joined to raw behaviour."""
    long = load_invivo_long()
    metadata = (
        long[["animal_id", "id", "phase", "condition"]]
        .drop_duplicates()
        .sort_values(["phase", "condition", "id"])
    )
    if metadata["animal_id"].duplicated().any():
        raise ValueError("Animal metadata differ between final qPCR sheets")

    expression = long.pivot(
        index="animal_id", columns="gene", values="relative_expression"
    ).reset_index()
    delta_ct = long.pivot(
        index="animal_id", columns="gene", values="delta_ct"
    ).add_suffix("_delta_ct").reset_index()
    output = metadata.merge(expression, on="animal_id", validate="one_to_one")
    output = output.merge(delta_ct, on="animal_id", validate="one_to_one")
    behavior = load_behavior_016()
    output = output.merge(
        behavior, on="animal_id", how="left", validate="one_to_one", indicator=True
    )
    if not output["_merge"].eq("both").all():
        ids = output.loc[output["_merge"] != "both", "animal_id"].tolist()
        raise ValueError(f"qPCR animals missing raw behaviour: {ids}")
    output = output.drop(columns="_merge")
    mismatch = (
        output["condition"].str.lower()
        != output["condition_behavior"].str.lower()
    )
    if mismatch.any():
        rows = output.loc[
            mismatch, ["animal_id", "condition", "condition_behavior"]
        ].to_dict("records")
        raise ValueError(f"qPCR/behaviour condition mismatch: {rows}")
    output["group"] = output["phase"] + " " + output["condition"]
    return output.sort_values(
        ["phase", "condition", "id"], kind="stable"
    ).reset_index(drop=True)
