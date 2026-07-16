"""Rebuild the selected-cell Sholl report directly from archived SWC files."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from ramification_index.ramification_index import (  # noqa: E402
    compute_ramification_index_for_swc_file,
)


MANIFEST = ROOT / "data" / "MorphOMICs" / "reports_cam" / "cam_manifest_selected.csv"
ARCHIVED_REPORT = ROOT / "data" / "MorphOMICs" / "reports_cam" / "cam_sholl_per_cell.csv"
OUTPUT = ROOT / "generated_outputs" / "MorphOMICs" / "cam_sholl_per_cell.csv"
METADATA_COLUMNS = [
    "AnimalID", "Condicion", "Ciclo", "Sexo", "Lado", "Capa",
    "CellIndex", "ImageIndex", "SWCPath",
]
METRIC_COLUMNS = [
    "ShollTotalIntersections", "ShollMaxIntersections",
    "ShollMeanIntersections", "ShollEnclosingRadius",
    "ShollRamificationIndex",
]


def build_report() -> pd.DataFrame:
    manifest = pd.read_csv(MANIFEST)
    selected = manifest.loc[manifest["selected"].astype(bool)].copy()
    records: list[dict[str, object]] = []
    for row in selected.itertuples(index=False):
        relative_path = Path(str(row.target_path))
        metrics = compute_ramification_index_for_swc_file(ROOT / relative_path)
        records.append(
            {
                "AnimalID": row.animal_id,
                "Condicion": row.Condicion,
                "Ciclo": row.Ciclo,
                "Sexo": row.Sexo,
                "Lado": row.lado,
                "Capa": row.capa,
                "CellIndex": row.cell_index,
                "ImageIndex": row.image_index,
                "SWCPath": relative_path.as_posix(),
                **metrics,
            }
        )
    return pd.DataFrame.from_records(records, columns=METADATA_COLUMNS + METRIC_COLUMNS)


def verify_against_archived(regenerated: pd.DataFrame) -> None:
    archived = pd.read_csv(ARCHIVED_REPORT)
    if archived[METADATA_COLUMNS].astype(str).to_dict("records") != regenerated[
        METADATA_COLUMNS
    ].astype(str).to_dict("records"):
        raise RuntimeError("Regenerated Sholl metadata does not match the archived report.")
    if not np.allclose(
        archived[METRIC_COLUMNS].to_numpy(float),
        regenerated[METRIC_COLUMNS].to_numpy(float),
        rtol=0.0,
        atol=1e-12,
        equal_nan=True,
    ):
        raise RuntimeError("Regenerated Sholl metrics do not match the archived report.")


def main() -> None:
    regenerated = build_report()
    verify_against_archived(regenerated)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    regenerated.to_csv(OUTPUT, index=False)
    print(f"[OK] Rebuilt and verified {len(regenerated)} SWC-derived Sholl rows")


if __name__ == "__main__":
    main()
