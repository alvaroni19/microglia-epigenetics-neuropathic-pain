#!/usr/bin/env python3
"""Build the animal-linked ramification and pain source table."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(os.path.abspath(__file__)).parent
ROOT = SCRIPT_DIR.parents[1]
SHOLL = (
    ROOT
    / "data"
    / "MorphOMICs"
    / "reports_cam"
    / "cam_sholl_per_cell.csv"
)
BEHAVIOUR = (
    ROOT
    / "data"
    / "Allodynia_0_16"
    / "allodynia_von_frey_scores.xlsx"
)
OUTPUT = (
    ROOT
    / "generated_outputs"
    / "MorphOMICs"
    / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
)
LEGACY_LINKAGE = (
    ROOT
    / "data"
    / "MorphOMICs"
    / "morphology_behaviour_linkage.xlsx"
)

METRICS = [
    "ShollTotalIntersections",
    "ShollMaxIntersections",
    "ShollMeanIntersections",
    "ShollEnclosingRadius",
    "ShollRamificationIndex",
]
PHASES = {
    "M": "Male",
    "Male": "Male",
    "Proestrus": "Proestrus",
    "Estrus": "Estrus",
    "Diestrus": "Diestrus",
    "Metestrus": "Diestrus",
}


def io_path(path: Path) -> Path:
    """Return a Windows extended-length path when the bundle is deeply nested."""
    if sys.platform == "win32" and not str(path).startswith("\\\\?\\"):
        return Path("\\\\?\\" + str(path.absolute()))
    return path


def canonical_animal_id(value: object) -> str:
    text = str(value).strip().upper()
    male = re.fullmatch(r"R(\d+[A-Z]?)M", text)
    if male:
        return f"RM{male.group(1)}"
    female = re.fullmatch(r"R(?:H)?(\d+[A-Z]?)", text)
    if female:
        return f"RH{female.group(1)}"
    raise ValueError(f"Unrecognized behavioural animal ID: {value!r}")


def main() -> None:
    sholl = pd.read_csv(io_path(SHOLL))
    sholl["animal_id"] = sholl["AnimalID"].astype(str).str.strip().str.upper()
    sholl["condition"] = sholl["Condition"].astype(str).str.strip()
    sholl["sex_or_phase"] = sholl["Phase"].replace(PHASES)
    sholl["side"] = sholl["Side"].astype(str).str.strip()
    sholl["lamina"] = pd.to_numeric(sholl["Layer"], errors="raise").astype(int)

    grouped = (
        sholl.groupby(
            ["animal_id", "condition", "sex_or_phase", "side", "lamina"],
            as_index=False,
        )
        .agg(**{metric: (metric, "mean") for metric in METRICS}, n_cells=("AnimalID", "size"))
    )

    legacy = pd.read_excel(io_path(LEGACY_LINKAGE))
    aliases: dict[str, str] = {}
    for value in legacy["ID"].dropna().astype(str):
        match = re.fullmatch(
            r"\s*(RH\d+[A-Z]?)\s*\((RH\d+)\)\s*",
            value.upper(),
        )
        if match:
            aliases[match.group(1)] = match.group(2)
    if aliases != {"RH55A": "RH55", "RH55B": "RH56"}:
        raise ValueError(f"Unexpected morphology-behaviour aliases: {aliases}")
    grouped["behaviour_animal_id"] = grouped["animal_id"].replace(aliases)
    grouped["linkage_method"] = grouped["animal_id"].map(
        lambda value: "validated_alias" if value in aliases else "exact_id"
    )

    behaviour = pd.read_excel(io_path(BEHAVIOUR), sheet_name="allodynia_statistics")
    behaviour["behaviour_animal_id"] = behaviour["n"].map(
        canonical_animal_id
    )
    behaviour_out = pd.DataFrame(
        {
            "behaviour_animal_id": behaviour["behaviour_animal_id"],
            "behaviour_condition": behaviour["Surgery"].astype(str).str.strip(),
            "phase_at_surgery": behaviour["phase_surgery"].replace(PHASES),
            "phase_at_test": behaviour["phase_test"].replace(PHASES),
            "baseline_0.16g": pd.to_numeric(behaviour["BL016"], errors="raise"),
            "test_0.16g": pd.to_numeric(behaviour["Test016"], errors="raise"),
        }
    )
    if behaviour_out["behaviour_animal_id"].duplicated().any():
        raise ValueError("Behavioural animal IDs are not unique.")
    behaviour_out["delta_0.16g"] = (
        behaviour_out["test_0.16g"] - behaviour_out["baseline_0.16g"]
    )

    output = grouped.merge(
        behaviour_out,
        on="behaviour_animal_id",
        how="inner",
        validate="many_to_one",
    )
    output["condition_matches_behaviour"] = (
        output["condition"].str.casefold()
        == output["behaviour_condition"].str.casefold()
    )
    output = output.sort_values(
        ["animal_id", "side", "lamina"],
        kind="stable",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(io_path(OUTPUT), index=False)

    mismatches = output.loc[
        ~output["condition_matches_behaviour"],
        ["animal_id", "condition", "behaviour_condition"],
    ].drop_duplicates()
    animal_count = output["animal_id"].nunique()
    print(
        f"[OK] Wrote {len(output)} rows for {animal_count} animals: {OUTPUT}"
    )
    if len(mismatches):
        print("[WARNING] Condition conflicts retained for author review:")
        print(mismatches.to_string(index=False))


if __name__ == "__main__":
    main()
