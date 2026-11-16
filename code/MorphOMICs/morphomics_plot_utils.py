"""Small plotting-only utilities for archived MorphOMICs objects and SWC files."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


SWC_COLUMNS = ["n", "type", "x", "y", "z", "radius", "parent"]


def read_swc_dataframe(filepath: Path) -> pd.DataFrame:
    swc = pd.read_csv(
        Path(filepath), sep=r"\s+", comment="#", names=SWC_COLUMNS, engine="python"
    )
    if swc.empty:
        return swc
    swc = swc.apply(pd.to_numeric, errors="coerce").dropna()
    swc["n"] = swc["n"].astype(int)
    swc["parent"] = swc["parent"].astype(int)
    return swc


def _total_length(filepath: Path) -> float:
    swc = read_swc_dataframe(filepath)
    coordinates = {
        int(row.n): np.array([row.x, row.y, row.z], dtype=float)
        for row in swc.itertuples(index=False)
    }
    lengths = []
    for row in swc.itertuples(index=False):
        node, parent = int(row.n), int(row.parent)
        if parent != -1 and node in coordinates and parent in coordinates:
            lengths.append(np.linalg.norm(coordinates[node] - coordinates[parent]))
    return float(np.sum(lengths))


def pick_representative_swc(filepaths: Iterable[Path]) -> Path | None:
    paths = [Path(path) for path in filepaths]
    if not paths:
        return None
    pairs = [(path, _total_length(path)) for path in paths]
    positive = [(path, length) for path, length in pairs if length > 0]
    if not positive:
        return paths[0]
    median = float(np.median([length for _, length in positive]))
    return min(positive, key=lambda item: abs(item[1] - median))[0]


def compute_persistence_centroids(
    frame: pd.DataFrame,
    group_columns: Iterable[str],
    vector_column: str = "pi",
) -> dict[str, dict[str, object]]:
    columns = [column for column in group_columns if column in frame.columns]
    if frame.empty or vector_column not in frame.columns or not columns:
        return {}
    centroids = {}
    for combination, subset in frame.groupby(columns):
        vectors = [
            np.asarray(value, dtype=float)
            for value in subset[vector_column].tolist()
            if value is not None
        ]
        if not vectors:
            continue
        stack = np.vstack(vectors)
        mean_vector = stack.mean(axis=0)
        side = int(round(np.sqrt(mean_vector.size)))
        if side * side != mean_vector.size:
            continue
        values = combination if isinstance(combination, tuple) else (combination,)
        metadata = dict(zip(columns, values))
        key = "/".join(map(str, values))
        centroids[key] = {
            "centroid": mean_vector.reshape(side, side),
            "n_samples": int(len(subset)),
            "combo": metadata,
        }
    return centroids
