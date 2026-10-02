"""Calculate Sholl profiles and ramification indices from SWC reconstructions."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

SWC_COLUMNS = ["n", "type", "x", "y", "z", "radius", "parent"]


def read_swc_dataframe(filepath: Path) -> pd.DataFrame:
    """Read an SWC file into a numeric dataframe, dropping malformed rows."""
    swc = pd.read_csv(
        filepath,
        sep=r"\s+",
        comment="#",
        names=SWC_COLUMNS,
        engine="python",
    )
    if swc.empty:
        return swc
    swc = swc.apply(pd.to_numeric, errors="coerce").dropna()
    swc["n"] = swc["n"].astype(int)
    swc["parent"] = swc["parent"].astype(int)
    return swc


def compute_sholl_profile_from_swc(
    swc: pd.DataFrame,
    radius_step: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Count process intersections with concentric radii centred on the soma.

    The soma is taken as the root node (parent == -1); if no root is flagged,
    the centroid of all points is used instead. Radii are stepped from
    `radius_step` up to the maximum node distance from the soma. An edge
    (parent-child pair) is counted as crossing a radius shell if it straddles
    that radius.
    """
    if swc.empty:
        return np.array([], dtype=float), np.array([], dtype=float)

    if (swc["parent"] == -1).any():
        soma = swc.loc[swc["parent"] == -1, ["x", "y", "z"]].iloc[0].to_numpy(dtype=float)
    else:
        soma = swc[["x", "y", "z"]].to_numpy(dtype=float).mean(axis=0)

    coords = swc[["x", "y", "z"]].to_numpy(dtype=float)
    distances = np.linalg.norm(coords - soma, axis=1)
    node_dist = {int(row.n): float(distances[idx]) for idx, row in enumerate(swc.itertuples(index=False))}

    max_distance = float(distances.max()) if len(distances) else 0.0
    if max_distance <= 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    radii = np.arange(radius_step, max_distance + radius_step, radius_step, dtype=float)
    intersections: List[int] = []
    for radius in radii:
        count = 0
        for row in swc.itertuples(index=False):
            parent_id = int(row.parent)
            node_id = int(row.n)
            if parent_id == -1 or parent_id not in node_dist or node_id not in node_dist:
                continue
            d_node = node_dist[node_id]
            d_parent = node_dist[parent_id]
            if (d_node >= radius > d_parent) or (d_parent >= radius > d_node):
                count += 1
        intersections.append(count)
    return radii, np.array(intersections, dtype=float)


def compute_sholl_metrics(radii: np.ndarray, intersections: np.ndarray) -> Dict[str, float]:
    """Summarize a Sholl profile, including the ramification index.

    ShollRamificationIndex = sum(intersections across all radii) / max(intersections at any single radius)
    """
    if len(intersections) == 0:
        return {
            "ShollTotalIntersections": 0.0,
            "ShollMaxIntersections": 0.0,
            "ShollMeanIntersections": 0.0,
            "ShollEnclosingRadius": 0.0,
            "ShollRamificationIndex": 0.0,
        }
    max_intersections = float(intersections.max())
    enclosing_radius = float(radii[np.argmax(intersections)]) if len(radii) else 0.0
    return {
        "ShollTotalIntersections": float(intersections.sum()),
        "ShollMaxIntersections": max_intersections,
        "ShollMeanIntersections": float(intersections.mean()),
        "ShollEnclosingRadius": enclosing_radius,
        "ShollRamificationIndex": float(intersections.sum() / max_intersections) if max_intersections > 0 else 0.0,
    }


def compute_ramification_index_for_swc_file(filepath: Path, radius_step: float = 5.0) -> Dict[str, float]:
    """Convenience wrapper: SWC file path -> Sholl/ramification-index metrics."""
    swc = read_swc_dataframe(Path(filepath))
    radii, intersections = compute_sholl_profile_from_swc(swc, radius_step=radius_step)
    return compute_sholl_metrics(radii, intersections)


def compute_ramification_index_for_swc_folder(
    folder: Path,
    radius_step: float = 5.0,
) -> pd.DataFrame:
    """Batch helper: compute Sholl/ramification-index metrics for every .swc file in a folder."""
    folder = Path(folder)
    records: List[Dict[str, object]] = []
    for filepath in sorted(folder.glob("*.swc")):
        metrics = compute_ramification_index_for_swc_file(filepath, radius_step=radius_step)
        record: Dict[str, object] = {"SWCPath": str(filepath)}
        record.update(metrics)
        records.append(record)
    return pd.DataFrame(records)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("swc_folder", type=Path, help="Folder containing .swc files")
    parser.add_argument("--radius-step", type=float, default=5.0, help="Sholl radius increment (default 5.0)")
    parser.add_argument("--out", type=Path, default=None, help="Optional output CSV path")
    args = parser.parse_args()

    result = compute_ramification_index_for_swc_folder(args.swc_folder, radius_step=args.radius_step)
    if args.out:
        result.to_csv(args.out, index=False)
        print(f"Wrote {len(result)} rows to {args.out}")
    else:
        print(result.to_string())
