from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


VALID_LAYERS = {"12", "34"}
SWC_COLUMNS = ["n", "type", "x", "y", "z", "radius", "parent"]
_FILENAME_INDEX_CACHE: Dict[str, Dict[str, List[Path]]] = {}

CYCLE_LABELS_EN = {
    "Proestrus": "Proestrus",
    "Estrus": "Estrus",
    "Diestrus": "Diestrus",
    "Metestrus": "Diestrus",
    "Male": "Male",
}
SEX_LABELS_EN = {"Female": "Female", "Male": "Male"}
SIDE_LABELS_EN = {"IPSI": "Ipsilateral", "CONTRA": "Contralateral"}
LAYER_LABELS_EN = {"12": "Layer 1/2", "34": "Layer 3/4"}


# Metadata map of the experimental groups.
ANIMAL_METADATA: Dict[str, Dict[str, str]] = {
    "RH12": {"Condition": "CCI", "Phase": "Proestrus", "Sex": "Female"},
    "RH15": {"Condition": "CCI", "Phase": "Diestrus", "Sex": "Female"},
    "RH16": {"Condition": "CCI", "Phase": "Proestrus", "Sex": "Female"},
    "RH17": {"Condition": "Control", "Phase": "Estrus", "Sex": "Female"},
    "RH21": {"Condition": "Control", "Phase": "Estrus", "Sex": "Female"},
    "RH22": {"Condition": "Control", "Phase": "Estrus", "Sex": "Female"},
    "RH25": {"Condition": "Control", "Phase": "Proestrus", "Sex": "Female"},
    "RH27": {"Condition": "Control", "Phase": "Diestrus", "Sex": "Female"},
    "RH30": {"Condition": "CCI", "Phase": "Diestrus", "Sex": "Female"},
    "RH31": {"Condition": "Control", "Phase": "Diestrus", "Sex": "Female"},
    "RH40": {"Condition": "CCI", "Phase": "Diestrus", "Sex": "Female"},
    "RH41": {"Condition": "CCI", "Phase": "Diestrus", "Sex": "Female"},
    "RH43": {"Condition": "Control", "Phase": "Proestrus", "Sex": "Female"},
    "RH44": {"Condition": "CCI", "Phase": "Diestrus", "Sex": "Female"},
    "RH45": {"Condition": "CCI", "Phase": "Diestrus", "Sex": "Female"},
    "RH51": {"Condition": "Control", "Phase": "Proestrus", "Sex": "Female"},
    "RH53": {"Condition": "CCI", "Phase": "Proestrus", "Sex": "Female"},
    "RH55a": {"Condition": "Control", "Phase": "Estrus", "Sex": "Female"},
    "RH55b": {"Condition": "CCI", "Phase": "Estrus", "Sex": "Female"},
    "RH59": {"Condition": "CCI", "Phase": "Estrus", "Sex": "Female"},
    "RH66": {"Condition": "CCI", "Phase": "Proestrus", "Sex": "Female"},
    "RH75": {"Condition": "CCI", "Phase": "Estrus", "Sex": "Female"},
    "RH78": {"Condition": "Control", "Phase": "Diestrus", "Sex": "Female"},
    "RH79": {"Condition": "Control", "Phase": "Diestrus", "Sex": "Female"},
    "RH83": {"Condition": "CCI", "Phase": "Estrus", "Sex": "Female"},
    "RH752": {"Condition": "CCI", "Phase": "Estrus", "Sex": "Female"},
    "RM2": {"Condition": "CCI", "Phase": "Male", "Sex": "Male"},
    "RM5": {"Condition": "CCI", "Phase": "Male", "Sex": "Male"},
    "RM19": {"Condition": "CCI", "Phase": "Male", "Sex": "Male"},
    "RM24": {"Condition": "CCI", "Phase": "Male", "Sex": "Male"},
    "RM25": {"Condition": "Control", "Phase": "Male", "Sex": "Male"},
    "RM37": {"Condition": "Control", "Phase": "Male", "Sex": "Male"},
    "RM41": {"Condition": "Control", "Phase": "Male", "Sex": "Male"},
    "RM43": {"Condition": "Control", "Phase": "Male", "Sex": "Male"},
}

# Known aliases in laboratory notes.
ANIMAL_ID_ALIASES = {
    "RH56": "RH55b",
    "RH55": "RH55a",
}

_META_BY_UPPER = {animal_id.upper(): meta for animal_id, meta in ANIMAL_METADATA.items()}
_ALIAS_BY_UPPER = {k.upper(): v for k, v in ANIMAL_ID_ALIASES.items()}


FILENAME_PATTERN = re.compile(
    r"^(?P<prefix>.+?)_CAM_(?P<ii>ii\d+)_(?P<sidechunk>[A-Z]+[-_]?\d+)_ProjMax(?P<proj>\d+)"
    r"_Image_(?P<image>\d+)_Filaments_(?P<rest>.+)\.swc$",
    re.IGNORECASE,
)


def _normalize_animal_id(animal_id_raw: str) -> str:
    clean = animal_id_raw.strip().replace(" ", "")
    upper = clean.upper()
    if upper in _ALIAS_BY_UPPER:
        return _ALIAS_BY_UPPER[upper]
    if upper in _META_BY_UPPER:
        return next(k for k in ANIMAL_METADATA if k.upper() == upper)
    return clean


def _metadata_for_animal(animal_id: str) -> Tuple[bool, Dict[str, str]]:
    upper = animal_id.upper()
    if upper in _META_BY_UPPER:
        return True, dict(_META_BY_UPPER[upper])
    inferred_sex = "Male" if upper.startswith("RM") else ("Female" if upper.startswith("RH") else "Unknown")
    inferred_cycle = "Male" if inferred_sex == "Male" else "Unknown"
    return False, {"Condition": "Unknown", "Phase": inferred_cycle, "Sex": inferred_sex}


def _parse_side(sidechunk_raw: str) -> Tuple[Optional[str], Optional[int], bool]:
    sidechunk = sidechunk_raw.upper().replace("-", "_")
    match = re.match(r"^(IPSI|CONTRA|CONTROL)_?(\d+)$", sidechunk)
    if match is None:
        return None, None, False
    side = match.group(1)
    side_image = int(match.group(2))
    if side == "CONTROL":
        side = "CONTRA"
    separator_present = "_" in sidechunk_raw
    return side, side_image, separator_present


def _parse_rest(rest: str) -> Tuple[Optional[int], Optional[str], Optional[str], List[str]]:
    issues: List[str] = []
    parts = rest.split("_")
    cell: Optional[int] = None
    layer: Optional[str] = None
    variant: Optional[str] = None

    if len(parts) == 1:
        token = parts[0]
        if token.isdigit() and token in VALID_LAYERS:
            layer = token
            issues.append("missing_cell_index")
        elif token.isdigit():
            cell = int(token)
            issues.append("missing_layer")
        else:
            issues.append("invalid_filaments_suffix")
    elif len(parts) >= 2:
        if parts[0].isdigit():
            cell = int(parts[0])
        else:
            issues.append("invalid_cell_index")
        if parts[1].isdigit():
            layer = parts[1]
        else:
            issues.append("invalid_layer")
        if len(parts) >= 3:
            variant = parts[2]
    else:
        issues.append("invalid_filaments_suffix")

    # Pattern observed in malformed entries such as "Filaments_34_02".
    if layer in {"01", "02", "03", "04", "05"} and cell is not None and str(cell) in VALID_LAYERS:
        issues.append("layer_variant_without_cell")
        variant = layer
        layer = str(cell)
        cell = None

    return cell, layer, variant, issues


def parse_swc_filename(filename: str) -> Dict[str, object]:
    name = Path(filename).name
    normalized_name = name.replace("-CAM_", "_CAM_")
    out: Dict[str, object] = {
        "source_name": name,
        "parse_ok": False,
        "parse_issue": "",
        "animal_id_raw": None,
        "animal_id": None,
        "marker": "CAM",
        "ii_label": None,
        "ii_index": None,
        "side": None,
        "side_image": None,
        "image_index": None,
        "cell_index": None,
        "layer": None,
        "variant": None,
        "canonical_key": None,
        "canonical_name": None,
        "prefix_token_count": None,
        "side_separator_present": None,
    }

    match = FILENAME_PATTERN.match(normalized_name)
    if match is None:
        out["parse_issue"] = "pattern_mismatch"
        return out

    groups = match.groupdict()
    prefix = groups["prefix"]
    prefix_tokens = prefix.split("_")
    animal_raw = prefix_tokens[0]
    animal_id = _normalize_animal_id(animal_raw)

    side, side_image, side_separator_present = _parse_side(groups["sidechunk"])
    if side is None:
        out["parse_issue"] = "invalid_side_chunk"
        out["animal_id_raw"] = animal_raw
        out["animal_id"] = animal_id
        return out

    cell_index, layer, variant, rest_issues = _parse_rest(groups["rest"])
    ii_label = groups["ii"].lower()
    ii_index = int(ii_label.replace("ii", ""))
    image_index = int(groups["image"])

    issues = list(rest_issues)
    if layer is not None and layer not in VALID_LAYERS:
        issues.append("layer_not_12_or_34")

    valid_cell_layer = cell_index is not None and layer in VALID_LAYERS
    canonical_key: Optional[str] = None
    canonical_name: Optional[str] = None
    if valid_cell_layer:
        canonical_key = f"{animal_id}|{ii_label}|{side}|{side_image}|{image_index}|{cell_index}|{layer}"
        canonical_name = (
            f"{animal_id}_CAM_{ii_label}_{side}_{side_image}_Image_{image_index}_"
            f"Cell_{cell_index}_Layer_{layer}.swc"
        )
    else:
        if cell_index is None:
            issues.append("missing_valid_cell")
        if layer not in VALID_LAYERS:
            issues.append("missing_valid_layer")

    out.update(
        {
            "parse_ok": True,
            "parse_issue": ";".join(sorted(set(issues))),
            "animal_id_raw": animal_raw,
            "animal_id": animal_id,
            "ii_label": ii_label,
            "ii_index": ii_index,
            "side": side,
            "side_image": side_image,
            "image_index": image_index,
            "cell_index": cell_index,
            "layer": layer,
            "variant": variant,
            "canonical_key": canonical_key,
            "canonical_name": canonical_name,
            "prefix_token_count": len(prefix_tokens),
            "side_separator_present": side_separator_present,
        }
    )
    return out


def parse_cam_filename(filename: str) -> Dict[str, object]:
    """Backward-compatible CAM parser name used by older notebooks."""
    return parse_swc_filename(filename)


def validate_swc_file(filepath: Path) -> Tuple[bool, str, int]:
    try:
        swc = pd.read_csv(
            filepath,
            sep=r"\s+",
            comment="#",
            names=SWC_COLUMNS,
            engine="python",
        )
    except Exception as exc:
        return False, f"read_error:{exc}", 0

    if swc.empty:
        return False, "empty_file", 0

    if swc.shape[1] != 7:
        return False, "invalid_column_count", int(len(swc))

    numeric = swc[SWC_COLUMNS].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        return False, "non_numeric_values", int(len(swc))

    if not np.isfinite(numeric[["x", "y", "z", "radius"]].to_numpy(dtype=float)).all():
        return False, "non_finite_values", int(len(swc))

    if (numeric["radius"] < 0).any():
        return False, "negative_radius", int(len(swc))

    node_ids = numeric["n"].astype(int)
    parents = numeric["parent"].astype(int)
    if node_ids.duplicated().any():
        return False, "duplicate_node_ids", int(len(swc))

    valid_nodes = set(node_ids.tolist())
    invalid_parent_count = int(((parents != -1) & (~parents.isin(valid_nodes))).sum())
    if invalid_parent_count > 0:
        return False, f"invalid_parent_refs:{invalid_parent_count}", int(len(swc))

    return True, "", int(len(swc))


def _rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    ranked = df.copy()
    ranked["q_variant"] = ranked["variant"].fillna("").ne("")
    ranked["q_prefix"] = ranked["prefix_token_count"].fillna(999).astype(int)
    ranked["q_side_sep"] = (~ranked["side_separator_present"].fillna(False)).astype(int)
    ranked["q_name_len"] = ranked["source_name"].astype(str).str.len()
    ranked = ranked.sort_values(
        ["q_variant", "q_prefix", "q_side_sep", "q_name_len", "source_name"],
        ascending=[True, True, True, True, True],
    )
    return ranked


def select_primary_rows(
    frame: pd.DataFrame,
    key_col: str = "canonical_key",
    eligible_col: Optional[str] = "eligible",
) -> pd.DataFrame:
    """Select one deterministic primary row per canonical morphology key.

    The ranking mirrors ``prepare_cam_dataset``: prefer non-variant names,
    short prefixes, explicit side separators, and finally lexical filename
    order. This avoids arbitrary duplicate resolution.
    """
    if frame.empty:
        return frame.copy()
    if key_col not in frame.columns:
        raise KeyError(f"Missing required key column: {key_col}")

    work = frame.copy()
    if eligible_col is not None and eligible_col in work.columns:
        work = work[work[eligible_col].fillna(False)].copy()
    work = work[work[key_col].notna()].copy()
    if work.empty:
        return work

    required_defaults = {
        "variant": None,
        "prefix_token_count": 999,
        "side_separator_present": False,
        "source_name": "",
    }
    for column, default in required_defaults.items():
        if column not in work.columns:
            work[column] = default

    ranked = _rank_candidates(work)
    return ranked.groupby(key_col, sort=False).head(1).copy()


def prepare_cam_dataset(
    source_dir: Path,
    output_dir: Path,
    report_dir: Path,
    include_unknown_animals: bool = False,
    overwrite_output: bool = True,
) -> Dict[str, object]:
    source_dir = Path(source_dir)
    output_dir = Path(output_dir)
    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(source_dir.glob("*.swc"))
    rows: List[Dict[str, object]] = []
    for filepath in files:
        row = parse_swc_filename(filepath.name)
        row["source_path"] = str(filepath)
        known_animal, animal_meta = _metadata_for_animal(str(row.get("animal_id", "")))
        row["known_animal"] = known_animal
        row.update(animal_meta)
        swc_valid, swc_issue, point_count = validate_swc_file(filepath)
        row["swc_valid"] = swc_valid
        row["swc_issue"] = swc_issue
        row["swc_points"] = point_count
        rows.append(row)

    manifest = pd.DataFrame(rows)
    if manifest.empty:
        raise ValueError(f"No SWC files were found in: {source_dir}")

    manifest["valid_cell_layer"] = (
        manifest["cell_index"].notna() & manifest["layer"].isin(sorted(VALID_LAYERS))
    )
    manifest["eligible"] = (
        manifest["parse_ok"]
        & manifest["valid_cell_layer"]
        & manifest["swc_valid"]
        & (manifest["known_animal"] | include_unknown_animals)
    )

    ranked = _rank_candidates(manifest[manifest["eligible"]].copy())
    selected_idx = ranked.groupby("canonical_key", sort=False).head(1).index

    manifest["selected"] = False
    manifest.loc[selected_idx, "selected"] = True

    manifest["exclude_reason"] = ""
    manifest.loc[~manifest["parse_ok"], "exclude_reason"] = manifest.loc[
        ~manifest["parse_ok"], "parse_issue"
    ].replace("", "parse_failed")
    mask_no_cell_layer = manifest["parse_ok"] & ~manifest["valid_cell_layer"]
    manifest.loc[mask_no_cell_layer, "exclude_reason"] = "invalid_cell_or_layer"
    mask_unknown = manifest["parse_ok"] & manifest["valid_cell_layer"] & ~manifest["known_animal"]
    if not include_unknown_animals:
        manifest.loc[mask_unknown, "exclude_reason"] = "unknown_animal_id"
    mask_invalid_swc = (
        manifest["parse_ok"]
        & manifest["valid_cell_layer"]
        & (manifest["known_animal"] | include_unknown_animals)
        & ~manifest["swc_valid"]
    )
    manifest.loc[mask_invalid_swc, "exclude_reason"] = "swc_validation_failed"
    duplicate_mask = manifest["eligible"] & ~manifest["selected"]
    manifest.loc[duplicate_mask, "exclude_reason"] = "duplicate_canonical_key"

    selected = manifest[manifest["selected"]].copy()

    if overwrite_output and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    target_paths: List[str] = []
    target_dirs: List[str] = []
    for idx, row in selected.iterrows():
        rel_dir = Path(
            str(row["Condition"]),
            str(row["Phase"]),
            str(row["Sex"]),
            str(row["side"]),
            str(row["layer"]),
            str(row["animal_id"]),
        )
        filename = str(row["canonical_name"])
        destination_dir = output_dir / rel_dir
        destination_dir.mkdir(parents=True, exist_ok=True)
        destination_file = destination_dir / filename
        shutil.copy2(Path(str(row["source_path"])), destination_file)

        target_dirs.append(str(rel_dir).replace("\\", "/"))
        target_paths.append(str(destination_file))
        selected.at[idx, "target_dir"] = target_dirs[-1]
        selected.at[idx, "target_path"] = target_paths[-1]

    manifest = manifest.merge(
        selected[["source_path", "target_dir", "target_path"]],
        on="source_path",
        how="left",
    )

    selected = manifest[manifest["selected"]].copy()
    excluded = manifest[~manifest["selected"]].copy()

    group_summary = (
        selected.groupby(["Condition", "Phase", "Sex", "side", "layer"], dropna=False)
        .size()
        .reset_index(name="n_swc")
        .rename(columns={"side": "Side", "layer": "Layer"})
        .sort_values(["Condition", "Phase", "Sex", "Side", "Layer"])
    )
    animal_summary = (
        selected.groupby(["animal_id", "Condition", "Phase", "Sex"], dropna=False)
        .size()
        .reset_index(name="n_swc")
        .rename(columns={"animal_id": "AnimalID"})
        .sort_values(["AnimalID"])
    )

    manifest_path = report_dir / "cam_manifest_all.csv"
    selected_path = report_dir / "cam_manifest_selected.csv"
    excluded_path = report_dir / "cam_manifest_excluded.csv"
    group_path = report_dir / "cam_group_summary.csv"
    animal_path = report_dir / "cam_animal_summary.csv"

    manifest.to_csv(manifest_path, index=False)
    selected.to_csv(selected_path, index=False)
    excluded.to_csv(excluded_path, index=False)
    group_summary.to_csv(group_path, index=False)
    animal_summary.to_csv(animal_path, index=False)

    return {
        "manifest": manifest,
        "selected": selected,
        "excluded": excluded,
        "group_summary": group_summary,
        "animal_summary": animal_summary,
        "paths": {
            "manifest_all": str(manifest_path),
            "manifest_selected": str(selected_path),
            "manifest_excluded": str(excluded_path),
            "group_summary": str(group_path),
            "animal_summary": str(animal_path),
            "organized_dataset": str(output_dir),
        },
        "counts": {
            "source_files": int(len(manifest)),
            "selected_files": int(len(selected)),
            "excluded_files": int(len(excluded)),
            "unknown_animals": int((~manifest["known_animal"]).sum()),
            "duplicates_removed": int((manifest["exclude_reason"] == "duplicate_canonical_key").sum()),
            "invalid_cell_or_layer": int((manifest["exclude_reason"] == "invalid_cell_or_layer").sum()),
        },
    }


def _infer_data_root_from_path(filepath: Path) -> Optional[Path]:
    parts = list(filepath.parts)
    for i, part in enumerate(parts):
        if part.lower() == "data_cam_morphomics":
            return Path(*parts[: i + 1])
    return None


def normalize_cycle_label(cycle: object) -> str:
    text = str(cycle)
    if text.lower() == "metestrus":
        return "Diestrus"
    return text


def spanish_to_english_labels(value: object) -> object:
    """Translate common CAM metadata labels while preserving unknown values."""
    if isinstance(value, pd.DataFrame):
        return add_display_columns(value)
    if isinstance(value, pd.Series):
        return value.map(spanish_to_english_labels)

    text = str(value)
    normalized_cycle = normalize_cycle_label(text)
    for mapping in (CYCLE_LABELS_EN, SEX_LABELS_EN, SIDE_LABELS_EN, LAYER_LABELS_EN):
        if normalized_cycle in mapping:
            return mapping[normalized_cycle]
        if text in mapping:
            return mapping[text]
    return value


def add_display_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Add English display columns without changing the analysis columns."""
    out = frame.copy()
    if "Phase" in out.columns:
        out["Phase"] = out["Phase"].map(normalize_cycle_label)
        out["Cycle_EN"] = out["Phase"].map(CYCLE_LABELS_EN).fillna(out["Phase"])
    if "Condition" in out.columns:
        out["Condition_EN"] = out["Condition"].astype(str)
    if "Sex" in out.columns:
        out["Sex_EN"] = out["Sex"].map(SEX_LABELS_EN).fillna(out["Sex"])
    if "Side" in out.columns:
        out["Side_EN"] = out["Side"].map(SIDE_LABELS_EN).fillna(out["Side"])
    elif "side" in out.columns:
        out["Side_EN"] = out["side"].map(SIDE_LABELS_EN).fillna(out["side"])
    if "Layer" in out.columns:
        out["Layer_EN"] = out["Layer"].astype(str).map(LAYER_LABELS_EN).fillna(out["Layer"].astype(str))
    elif "layer" in out.columns:
        out["Layer_EN"] = out["layer"].astype(str).map(LAYER_LABELS_EN).fillna(out["layer"].astype(str))

    group_parts = [c for c in ["Condition_EN", "Cycle_EN", "Side_EN", "Layer_EN"] if c in out.columns]
    if group_parts:
        out["Group_EN"] = out[group_parts].astype(str).agg(" / ".join, axis=1)
    return out


def _read_report_csv(report_dir: Path, filename: str) -> pd.DataFrame:
    path = Path(report_dir) / filename
    if not path.exists():
        raise FileNotFoundError(f"Required report CSV not found: {path}")
    out = pd.read_csv(path)
    if "Phase" in out.columns:
        out["Phase"] = out["Phase"].map(normalize_cycle_label)
    if "Layer" in out.columns:
        out["Layer"] = out["Layer"].astype(str)
    return out


def load_filtered_morphomics(report_dir: Path) -> pd.DataFrame:
    """Load per-cell morphometrics generated by notebook 2."""
    return _read_report_csv(Path(report_dir), "cam_morphometrics_per_cell.csv")


def load_sholl_tables(report_dir: Path) -> pd.DataFrame:
    """Load per-cell Sholl metrics generated by notebook 2."""
    return _read_report_csv(Path(report_dir), "cam_sholl_per_cell.csv")


def compute_ramification_index(frame: pd.DataFrame) -> pd.DataFrame:
    """Add a Sholl ramification index column when Sholl metrics are present."""
    out = frame.copy()
    if "ShollRamificationIndex" in out.columns:
        out["RamificationIndex"] = pd.to_numeric(out["ShollRamificationIndex"], errors="coerce")
        return out
    required = {"ShollTotalIntersections", "ShollMaxIntersections"}
    if not required.issubset(out.columns):
        missing = sorted(required - set(out.columns))
        raise KeyError(f"Missing Sholl columns for ramification index: {missing}")
    total = pd.to_numeric(out["ShollTotalIntersections"], errors="coerce")
    max_intersections = pd.to_numeric(out["ShollMaxIntersections"], errors="coerce")
    out["RamificationIndex"] = np.where(max_intersections > 0, total / max_intersections, np.nan)
    return out


def compute_layer_animal_summary(
    frame: pd.DataFrame,
    metric_cols: Optional[Iterable[str]] = None,
    group_cols: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Average cell-level measurements to animal-level units by layer."""
    if frame.empty:
        return frame.copy()
    if group_cols is None:
        group_cols = [c for c in ["AnimalID", "Condition", "Phase", "Sex", "Layer"] if c in frame.columns]
    group_cols = list(group_cols)
    if "AnimalID" not in group_cols:
        raise KeyError("AnimalID is required for animal-level summaries.")

    if metric_cols is None:
        exclude = set(group_cols) | {"CellIndex", "ImageIndex", "SWCPath", "Radii", "Intersections"}
        metric_cols = [
            c for c in frame.columns
            if c not in exclude and pd.api.types.is_numeric_dtype(frame[c])
        ]
    metric_cols = list(metric_cols)
    if not metric_cols:
        raise ValueError("No numeric metric columns were available for summarization.")

    out = (
        frame.groupby(group_cols, dropna=False)[metric_cols]
        .mean()
        .reset_index()
    )
    n_cells = frame.groupby(group_cols, dropna=False).size().reset_index(name="N_Cells")
    out = out.merge(n_cells, on=group_cols, how="left")
    if "Phase" in out.columns:
        out["Phase"] = out["Phase"].map(normalize_cycle_label)
    return out


def compute_animal_level_pca(
    frame: pd.DataFrame,
    metric_cols: Optional[Iterable[str]] = None,
    group_cols: Optional[Iterable[str]] = None,
    n_components: int = 2,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, object]:
    """Compute PCA after collapsing cells to animal-level summaries."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    summary = compute_layer_animal_summary(frame, metric_cols=metric_cols, group_cols=group_cols)
    if metric_cols is None:
        exclude = {"N_Cells"}
        metric_cols = [
            c for c in summary.columns
            if c not in exclude and pd.api.types.is_numeric_dtype(summary[c])
        ]
    metric_cols = list(metric_cols)
    usable = summary.dropna(subset=metric_cols).copy()
    if usable.empty:
        raise ValueError("No complete animal-level rows are available for PCA.")

    n_components_eff = min(int(n_components), len(metric_cols), len(usable))
    scaler = StandardScaler()
    matrix = scaler.fit_transform(usable[metric_cols].to_numpy(dtype=float))
    pca = PCA(n_components=n_components_eff, random_state=random_state)
    scores = pca.fit_transform(matrix)
    for idx in range(n_components_eff):
        usable[f"PC{idx + 1}"] = scores[:, idx]
    usable.attrs["pca_explained_variance_ratio"] = pca.explained_variance_ratio_
    usable.attrs["pca_metric_cols"] = metric_cols
    usable.attrs["pca_scaler"] = scaler
    return usable, pca


def _get_filename_index(data_root: Path) -> Dict[str, List[Path]]:
    key = str(data_root.resolve()).lower()
    if key not in _FILENAME_INDEX_CACHE:
        index: Dict[str, List[Path]] = {}
        for swc_path in data_root.rglob("*.swc"):
            index.setdefault(swc_path.name, []).append(swc_path)
        _FILENAME_INDEX_CACHE[key] = index
    return _FILENAME_INDEX_CACHE[key]


def resolve_swc_path(filepath: Path) -> Path:
    filepath = Path(filepath)
    if filepath.exists():
        return filepath

    # Compatibility fix for old manifests that still reference Metestrus.
    text_path = str(filepath)
    replacement_candidates = [
        text_path.replace("\\Metestrus\\", "\\Diestrus\\"),
        text_path.replace("/Metestrus/", "/Diestrus/"),
        text_path.replace("Metestrus", "Diestrus"),
    ]
    for candidate in replacement_candidates:
        candidate_path = Path(candidate)
        if candidate_path.exists():
            return candidate_path

    # Remap UNC / old-drive paths: find "data_cam_morphomics" anchor in the
    # stored path and rebuild from the current working directory.
    parts = list(filepath.parts)
    try:
        idx = next(i for i, p in enumerate(parts) if p.lower() == "data_cam_morphomics")
        relative_suffix = Path(*parts[idx:])   # data_cam_morphomics/CCI/…/file.swc
        for base in (Path.cwd(), Path.cwd().parent, Path.cwd().parent.parent):
            rebuilt = base / relative_suffix
            if rebuilt.exists():
                return rebuilt
    except (StopIteration, TypeError):
        pass

    data_root = _infer_data_root_from_path(filepath)
    if data_root is not None and data_root.exists():
        index = _get_filename_index(data_root)
        same_name = index.get(filepath.name, [])
        if len(same_name) == 1:
            return same_name[0]
        if len(same_name) > 1:
            # Pick the shortest path as a deterministic fallback.
            return sorted(same_name, key=lambda p: (len(str(p)), str(p)))[0]

    raise FileNotFoundError(f"No SWC file found for path: {filepath}")


def read_swc_dataframe(filepath: Path) -> pd.DataFrame:
    filepath = resolve_swc_path(filepath)
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


def _edge_lengths(swc: pd.DataFrame) -> np.ndarray:
    if swc.empty:
        return np.array([], dtype=float)
    coords = {
        int(row.n): np.array([row.x, row.y, row.z], dtype=float)
        for row in swc.itertuples(index=False)
    }
    lengths: List[float] = []
    for row in swc.itertuples(index=False):
        node_id = int(row.n)
        parent_id = int(row.parent)
        if parent_id == -1:
            continue
        if node_id not in coords or parent_id not in coords:
            continue
        lengths.append(float(np.linalg.norm(coords[node_id] - coords[parent_id])))
    return np.array(lengths, dtype=float)


def estimate_total_length(filepath: Path) -> float:
    swc = read_swc_dataframe(filepath)
    return float(_edge_lengths(swc).sum())


def compute_morphometrics_from_swc(swc: pd.DataFrame) -> Dict[str, float]:
    if swc.empty:
        return {
            "NumPoints": 0,
            "NumRoots": 0,
            "NumEdges": 0,
            "NumBranchPoints": 0,
            "NumTips": 0,
            "TotalLength": 0.0,
            "MeanSegmentLength": 0.0,
            "MaxSegmentLength": 0.0,
            "MaxRadialDistance": 0.0,
            "MeanRadius": 0.0,
        }

    node_ids = swc["n"].astype(int).tolist()
    root_ids = set(swc.loc[swc["parent"] == -1, "n"].astype(int).tolist())
    children_count = {node_id: 0 for node_id in node_ids}
    for parent_id in swc["parent"].astype(int).tolist():
        if parent_id in children_count:
            children_count[parent_id] += 1

    num_branch_points = int(sum(v > 1 for v in children_count.values()))
    num_tips = int(sum((v == 0) and (k not in root_ids) for k, v in children_count.items()))

    seg_lengths = _edge_lengths(swc)
    total_length = float(seg_lengths.sum()) if len(seg_lengths) else 0.0
    mean_seg_length = float(seg_lengths.mean()) if len(seg_lengths) else 0.0
    max_seg_length = float(seg_lengths.max()) if len(seg_lengths) else 0.0

    coords = swc[["x", "y", "z"]].to_numpy(dtype=float)
    if len(root_ids):
        soma_point = swc.loc[swc["n"].isin(root_ids), ["x", "y", "z"]].iloc[0].to_numpy(dtype=float)
    else:
        soma_point = coords.mean(axis=0)
    radial = np.linalg.norm(coords - soma_point, axis=1)

    return {
        "NumPoints": float(len(swc)),
        "NumRoots": float(len(root_ids)),
        "NumEdges": float(len(seg_lengths)),
        "NumBranchPoints": float(num_branch_points),
        "NumTips": float(num_tips),
        "TotalLength": total_length,
        "MeanSegmentLength": mean_seg_length,
        "MaxSegmentLength": max_seg_length,
        "MaxRadialDistance": float(radial.max()) if len(radial) else 0.0,
        "MeanRadius": float(swc["radius"].mean()),
    }


def compute_morphometrics_for_manifest(selected_manifest: pd.DataFrame) -> pd.DataFrame:
    records: List[Dict[str, object]] = []
    for row in selected_manifest.itertuples(index=False):
        filepath = resolve_swc_path(Path(str(row.target_path)))
        swc = read_swc_dataframe(filepath)
        metrics = compute_morphometrics_from_swc(swc)
        record = {
            "AnimalID": row.animal_id,
            "Condition": row.Condition,
            "Phase": normalize_cycle_label(row.Phase),
            "Sex": row.Sex,
            "Side": row.side,
            "Layer": row.layer,
            "CellIndex": row.cell_index,
            "ImageIndex": row.image_index,
            "SWCPath": str(filepath),
        }
        record.update(metrics)
        records.append(record)
    return pd.DataFrame(records)


def compute_sholl_profile_from_swc(
    swc: pd.DataFrame,
    radius_step: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray]:
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


def compute_sholl_for_manifest(
    selected_manifest: pd.DataFrame,
    radius_step: float = 5.0,
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []
    for row in selected_manifest.itertuples(index=False):
        filepath = resolve_swc_path(Path(str(row.target_path)))
        swc = read_swc_dataframe(filepath)
        radii, intersections = compute_sholl_profile_from_swc(swc, radius_step=radius_step)
        metrics = compute_sholl_metrics(radii, intersections)
        record = {
            "AnimalID": row.animal_id,
            "Condition": row.Condition,
            "Phase": normalize_cycle_label(row.Phase),
            "Sex": row.Sex,
            "Side": row.side,
            "Layer": row.layer,
            "CellIndex": row.cell_index,
            "ImageIndex": row.image_index,
            "SWCPath": str(filepath),
            "Radii": radii,
            "Intersections": intersections,
        }
        record.update(metrics)
        records.append(record)
    return pd.DataFrame(records)


def make_group_label(df: pd.DataFrame, columns: Iterable[str], new_column: str = "Group") -> pd.DataFrame:
    columns = list(columns)
    out = df.copy()
    out[new_column] = out[columns].astype(str).agg("/".join, axis=1)
    return out


def pick_representative_swc(filepaths: Iterable[Path]) -> Optional[Path]:
    paths = [Path(p) for p in filepaths]
    if not paths:
        return None
    length_pairs: List[Tuple[Path, float]] = []
    for filepath in paths:
        try:
            total_length = estimate_total_length(filepath)
            if total_length > 0:
                length_pairs.append((filepath, total_length))
        except Exception:
            continue
    if not length_pairs:
        return paths[0]
    lengths = np.array([x[1] for x in length_pairs], dtype=float)
    median = float(np.median(lengths))
    return min(length_pairs, key=lambda x: abs(x[1] - median))[0]


def compute_persistence_centroids(
    frame: pd.DataFrame,
    group_cols: Iterable[str],
    vector_col: str = "pi",
) -> Dict[str, Dict[str, object]]:
    if frame is None or frame.empty or vector_col not in frame.columns:
        return {}
    group_cols = list(group_cols)
    available_cols = [c for c in group_cols if c in frame.columns]
    if not available_cols:
        return {}

    centroids: Dict[str, Dict[str, object]] = {}
    for combo, sub in frame.groupby(available_cols):
        vectors = [np.asarray(v, dtype=float) for v in sub[vector_col].tolist() if v is not None]
        if not vectors:
            continue
        stack = np.vstack(vectors)
        mean_vec = stack.mean(axis=0)
        std_vec = stack.std(axis=0)
        side = int(round(np.sqrt(mean_vec.size)))
        if side * side != mean_vec.size:
            continue
        if isinstance(combo, tuple):
            key = "/".join(map(str, combo))
            combo_dict = {available_cols[i]: combo[i] for i in range(len(available_cols))}
        else:
            key = str(combo)
            combo_dict = {available_cols[0]: combo}
        centroids[key] = {
            "centroid": mean_vec.reshape(side, side),
            "std": std_vec.reshape(side, side),
            "n_samples": int(len(sub)),
            "combo": combo_dict,
        }
    return centroids
