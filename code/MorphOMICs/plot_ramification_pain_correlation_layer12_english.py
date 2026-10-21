from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).absolute().parents[2]
sys.path.append(str(ROOT / "code" / "_helpers"))
import journal_style as ns  # noqa: E402
SOURCE = ROOT / "generated_outputs" / "MorphOMICs" / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
OUT_TIFF = ROOT / "generated_outputs" / "MorphOMICs" / "cam_ramification_allodynia_correlation_layer12_IPSI_publication.tiff"
OUT_PDF = ROOT / "generated_outputs" / "MorphOMICs" / "cam_ramification_allodynia_correlation_layer12_IPSI_publication.pdf"
OUT_SVG = ROOT / "generated_outputs" / "MorphOMICs" / "cam_ramification_allodynia_correlation_layer12_IPSI_publication.svg"
OUT_STATS = ROOT / "generated_outputs" / "MorphOMICs" / "cam_ramification_allodynia_correlation_layer12_IPSI_stats.csv"

FEATURE = "ShollRamificationIndex"
PAIN = "test_0.16g"
GROUP = "sex_or_phase"
CONDITION = "condition"


def io_path(path: Path) -> Path:
    """Return a Windows extended-length path when the bundle is deeply nested."""
    if sys.platform == "win32" and not str(path).startswith("\\\\?\\"):
        return Path("\\\\?\\" + str(path.absolute()))
    return path


def residualize_on_categoricals(values: pd.Series, covariates: pd.DataFrame) -> tuple[np.ndarray, int]:
    ranks = values.rank(method="average").to_numpy(dtype=float)
    designs = [
        pd.get_dummies(covariates[col].astype(str), prefix=col, drop_first=True, dtype=float)
        for col in covariates.columns
    ]
    covariate_design = pd.concat(designs, axis=1) if designs else pd.DataFrame(index=covariates.index)
    design = np.column_stack([np.ones(len(covariates)), covariate_design.to_numpy(dtype=float)])
    coefficients = np.linalg.lstsq(design, ranks, rcond=None)[0]
    residuals = ranks - design @ coefficients
    return residuals, int(np.linalg.matrix_rank(design))


def partial_spearman(
    frame: pd.DataFrame,
    feature: str,
    pain: str,
    covariate_cols: list[str],
) -> tuple[float, float, int, np.ndarray, np.ndarray]:
    tmp = frame[[feature, pain, *covariate_cols]].dropna()
    feature_resid, design_rank = residualize_on_categoricals(tmp[feature], tmp[covariate_cols])
    pain_resid, _ = residualize_on_categoricals(tmp[pain], tmp[covariate_cols])
    rho = float(np.corrcoef(feature_resid, pain_resid)[0, 1])
    df = tmp.shape[0] - design_rank - 1
    if df <= 0 or not np.isfinite(rho) or abs(rho) >= 1:
        p_value = np.nan
    else:
        t_stat = rho * np.sqrt(df / max(1e-12, 1 - rho**2))
        p_value = float(2 * stats.t.sf(abs(t_stat), df=df))
    return rho, p_value, df, feature_resid, pain_resid


def add_trend(ax: plt.Axes, x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> None:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_arr.size < 3 or np.unique(x_arr).size < 2:
        return
    slope, intercept = np.polyfit(x_arr, y_arr, deg=1)
    x_line = np.linspace(x_arr.min(), x_arr.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, color="#202020", lw=1.0, zorder=1)


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "p = NA"
    if value < 0.0001:
        return "p < 0.0001"
    return f"p = {value:.4f}" if value < 0.1 else f"p = {value:.3f}"


def main() -> None:
    source = SOURCE
    if not io_path(source).exists():
        source = (
            ROOT
            / "data"
            / "MorphOMICs"
            / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
        )
    if not io_path(source).exists():
        raise FileNotFoundError(f"Missing source CSV: {source}")

    data = pd.read_csv(io_path(source))
    data = data[(data["side"] == "IPSI") & (data["lamina"] == 12)].copy()
    data = data[["animal_id", CONDITION, GROUP, FEATURE, PAIN]].dropna()
    data[GROUP] = data[GROUP].replace({"Metestrus": "Diestrus"})
    data = data.drop_duplicates(subset=["animal_id"]).reset_index(drop=True)

    raw = stats.spearmanr(data[FEATURE], data[PAIN])
    group_rho, group_p, group_df, feature_resid, pain_resid = partial_spearman(
        data, FEATURE, PAIN, [GROUP]
    )
    full_rho, full_p, full_df, full_feature_resid, full_pain_resid = partial_spearman(
        data, FEATURE, PAIN, [GROUP, CONDITION]
    )

    stats_out = pd.DataFrame(
        [
            {
                "comparison": "raw",
                "adjusted_for": "none",
                "n_animals": data.shape[0],
                "rho": float(raw.statistic),
                "p_value": float(raw.pvalue),
                "degrees_of_freedom": np.nan,
            },
            {
                "comparison": "partial_spearman",
                "adjusted_for": "sex_estrous_group",
                "n_animals": data.shape[0],
                "rho": group_rho,
                "p_value": group_p,
                "degrees_of_freedom": group_df,
            },
            {
                "comparison": "partial_spearman",
                "adjusted_for": "sex_estrous_group_and_condition",
                "n_animals": data.shape[0],
                "rho": full_rho,
                "p_value": full_p,
                "degrees_of_freedom": full_df,
            },
        ]
    )
    OUT_STATS.parent.mkdir(parents=True, exist_ok=True)
    stats_out.to_csv(io_path(OUT_STATS), index=False)

    ns.set_style()  # Arial, 5-7 pt, journal theme

    condition_colors = {"Control": "#999999", "CCI": "#CC79A7"}  # unified condition
    group_markers = {
        "Male": "s",
        "Proestrus": "o",
        "Estrus": "^",
        "Diestrus": "D",
    }
    group_labels = {
        "Male": "Male",
        "Proestrus": "Proestrus",
        "Estrus": "Estrus",
        "Diestrus": "Diestrus",
    }

    fig, ax = plt.subplots(1, 1, figsize=(ns.mm2inch(92), ns.mm2inch(60)),
                           constrained_layout=False)
    plot_data = data.copy()
    # The displayed analysis follows the prespecified figure: ranks are
    # adjusted for sex/estrous phase.  Adjustment for experimental condition is
    # reported separately as a sensitivity analysis in the companion table.
    plot_data["feature_resid"] = feature_resid
    plot_data["pain_resid"] = pain_resid

    for _, row in plot_data.iterrows():
        ax.scatter(
            row["feature_resid"],
            row["pain_resid"],
            s=22,
            marker=group_markers.get(row[GROUP], "o"),
            facecolor=condition_colors.get(row[CONDITION], "#666666"),
            edgecolor="black",
            linewidth=0.4,
            alpha=0.95,
            zorder=3,
        )
    add_trend(ax, plot_data["feature_resid"], plot_data["pain_resid"])
    ax.set_title(
        f"Sex/estrous phase-adjusted ranks\n"
        f"Partial Spearman rho = {group_rho:.2f}, {format_p(group_p)}; n = {data.shape[0]}",
        pad=4, fontsize=ns.FS_TITLE, fontweight="bold",
    )
    ax.set_xlabel("Ramification residual rank", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.set_ylabel("Pain residual rank", fontsize=ns.FS_AXIS, fontweight="bold")
    ax.grid(True, color="#E6E6E6", linewidth=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    condition_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="black",
            label=label,
            markersize=4,
        )
        for label, color in condition_colors.items()
    ]
    group_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=marker,
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="black",
            label=group_labels.get(group, group),
            markersize=4,
        )
        for group, marker in group_markers.items()
        if group in set(data[GROUP])
    ]
    leg1 = ax.legend(
        handles=condition_handles,
        title="Condition",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.72),
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=group_handles,
        title="Sex/estrous phase",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.36),
    )

    fig.subplots_adjust(left=0.14, right=0.70, bottom=0.18, top=0.82)
    OUT_TIFF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(io_path(OUT_TIFF), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(io_path(OUT_PDF))
    fig.savefig(io_path(OUT_SVG))
    plt.close(fig)

    print(f"Saved {OUT_TIFF}")
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_STATS}")


if __name__ == "__main__":
    main()
