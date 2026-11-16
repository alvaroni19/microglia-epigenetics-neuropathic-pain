#!/usr/bin/env python3
"""Build all animal-level in vivo RT-qPCR sources from final Delta-Ct data."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
from qpcr_raw import (  # noqa: E402
    ESTROGEN_SCORE_GENES,
    MFA_GENES,
    load_invivo_long,
    load_invivo_wide,
)

OUT_QPCR = ROOT / "generated_outputs" / "PCR_260514"
OUT_MFA = ROOT / "generated_outputs" / "PCR_MFA"


def bh_adjust(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def marker_correlations(wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # Pain/microglial and DNA-methylation markers only; the estrogen-related genes
    # are analysed separately through the estrogenic score.
    for gene in [g for g in MFA_GENES if g not in ESTROGEN_SCORE_GENES]:
        complete = wide[[gene, "pain_016_delta"]].dropna()
        result = stats.spearmanr(
            complete[gene], complete["pain_016_delta"]
        )
        rows.append(
            {
                "marker": gene,
                "n": len(complete),
                "rho": float(result.statistic),
                "p_value": float(result.pvalue),
                "gene_label": gene,
            }
        )
    output = pd.DataFrame(rows)
    output["q_BH"] = bh_adjust(output["p_value"])
    fisher_z = np.arctanh(output["rho"].clip(-0.999999, 0.999999))
    fisher_se = 1.0 / np.sqrt(output["n"] - 3)
    output["approx_ci95_low"] = np.tanh(fisher_z - 1.96 * fisher_se)
    output["approx_ci95_high"] = np.tanh(fisher_z + 1.96 * fisher_se)
    output["ci_method"] = "Fisher-z approximation"
    return output.sort_values("rho", ascending=False).reset_index(drop=True)


def estrogenic_score(wide: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scored = wide.copy()
    z_columns = []
    for gene in ESTROGEN_SCORE_GENES:
        column = f"{gene}_z"
        scored[column] = (
            scored[gene] - scored[gene].mean(skipna=True)
        ) / scored[gene].std(skipna=True, ddof=1)
        z_columns.append(column)
    scored["estrogenic_score"] = scored[z_columns].mean(axis=1, skipna=True)
    scored["estrogenic_score_genes_observed"] = scored[z_columns].notna().sum(axis=1)

    female_controls = scored[
        (scored["phase"] != "Male") & (scored["condition"] == "Control")
    ]
    reference = female_controls.groupby("phase").agg(
        control_median_estrogenic_score=("estrogenic_score", "median"),
        control_median_pain_016_delta=("pain_016_delta", "median"),
    )
    cci = scored[
        (scored["phase"] != "Male") & (scored["condition"] == "CCI")
    ].copy()
    cci = cci.merge(reference, left_on="phase", right_index=True, validate="many_to_one")
    cci["delta_estrogenic_score"] = (
        cci["estrogenic_score"] - cci["control_median_estrogenic_score"]
    )
    cci["delta_pain_016"] = (
        cci["pain_016_delta"] - cci["control_median_pain_016_delta"]
    )
    complete = cci[
        ["delta_estrogenic_score", "delta_pain_016"]
    ].dropna()
    result = stats.spearmanr(
        complete["delta_estrogenic_score"], complete["delta_pain_016"]
    )
    summary = pd.DataFrame(
        [
            {
                "method": "Spearman",
                "rho": float(result.statistic),
                "p_value": float(result.pvalue),
                "n": len(complete),
            }
        ]
    )
    cci = cci.rename(
        columns={
            "phase": "phase",
            "condition": "condition",
            "group": "group",
            "pain_016_delta": "pain016",
        }
    )
    return cci, summary


def main() -> None:
    OUT_QPCR.mkdir(parents=True, exist_ok=True)
    OUT_MFA.mkdir(parents=True, exist_ok=True)
    long = load_invivo_long()
    wide = load_invivo_wide()
    if len(wide) != 34 or wide["animal_id"].nunique() != 34:
        raise ValueError("Expected 34 unique in vivo qPCR animals")

    long.to_csv(OUT_QPCR / "invivo_qpcr_from_delta_ct_long.csv", index=False)
    wide.to_csv(OUT_QPCR / "invivo_qpcr_animal_level.csv", index=False)
    wide.to_csv(OUT_MFA / "invivo_qpcr_mfa_input.csv", index=False)

    correlations = marker_correlations(wide)
    correlations.to_csv(
        OUT_MFA / "invivo_marker_pain_spearman.csv", index=False
    )
    cci_score, score_summary = estrogenic_score(wide)
    cci_score.to_csv(
        OUT_MFA / "invivo_estrogenic_score_delta_pain_per_animal.csv",
        index=False,
    )
    score_summary.to_csv(
        OUT_MFA / "invivo_estrogenic_score_delta_pain_correlations.csv",
        index=False,
    )

    calibrators = (
        long[
            [
                "gene",
                "source_sheet",
                "calibrator_mean_delta_ct",
            ]
        ]
        .drop_duplicates()
        .sort_values("gene")
    )
    calibrators.to_csv(
        OUT_QPCR / "invivo_qpcr_delta_delta_ct_calibrators.csv", index=False
    )
    print(
        "[OK] final in vivo qPCR sources: "
        f"{len(wide)} animals, {len(calibrators)} genes, "
        f"{len(correlations)} MFA/correlation markers"
    )
    print(
        "[OK] excluded from all final MFA/correlation sources: "
        "prkcg, p2rx4, AR, P4R, hsd17b3"
    )


if __name__ == "__main__":
    main()
