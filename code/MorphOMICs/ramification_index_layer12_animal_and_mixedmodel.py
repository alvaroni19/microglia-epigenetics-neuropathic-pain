#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ramification index in laminae 1/2: cell-level versus animal-level inference.

Cells are nested within animals (n = 2-6 animals per group, 50-139 cells per
group), so a cell-level Mann-Whitney test treats non-independent cells as
replicates and inflates significance. Layer 12 is analysed in three ways and the
results are reported side by side:
  (1) cell-level Mann-Whitney        -> reference only (pseudoreplicated)
  (2) animal-level Mann-Whitney on per-animal means -> correct unit
  (3) linear mixed model  RI ~ Condition (+ Side) with (1 | AnimalID)
      -> uses all cells while modelling the nesting

A diagnostic figure shows cells (grey) and per-animal means (black), Control vs
CCI, per phase.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR  = Path(__file__).absolute().parent
BUNDLE_ROOT = SCRIPT_DIR.parent.parent
PERCELL = BUNDLE_ROOT / "data" / "MorphOMICs" / "reports_cam" / "cam_sholl_per_cell.csv"
OUT_DIR = BUNDLE_ROOT / "generated_outputs" / "MorphOMICs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RI = "ShollRamificationIndex"
LAYER = 12
PHASE_ORDER = ["Proestrus", "Estrus", "Diestrus", "Male"]

df = pd.read_csv(PERCELL)
df = df[(df["Layer"] == LAYER) & (df["Side"] == "IPSI")].copy()  # ipsilateral cells only
df = df.dropna(subset=[RI])
df["Condition"] = df["Condition"].astype(str)
df["Phase"] = df["Phase"].astype(str)
df["Side"] = df["Side"].astype(str)

phases = [c for c in PHASE_ORDER if c in set(df["Phase"])]

# Per-animal means (cells from both sides are averaged, as in the plotted panel).
animal = (
    df.groupby(["Phase", "Condition", "AnimalID"], as_index=False)[RI]
      .mean()
      .rename(columns={RI: "RI_animal_mean"})
)

rows = []
for phase in phases:
    sub = df[df["Phase"] == phase]
    cci_cells = sub.loc[sub["Condition"] == "CCI", RI].dropna()
    ctl_cells = sub.loc[sub["Condition"] == "Control", RI].dropna()

    # (1) cell-level Mann-Whitney (pseudoreplicado, referencia)
    p_cell = np.nan
    if len(cci_cells) >= 3 and len(ctl_cells) >= 3:
        _, p_cell = stats.mannwhitneyu(cci_cells, ctl_cells, alternative="two-sided")

    # (2) animal-level
    a = animal[animal["Phase"] == phase]
    cci_an = a.loc[a["Condition"] == "CCI", "RI_animal_mean"].dropna()
    ctl_an = a.loc[a["Condition"] == "Control", "RI_animal_mean"].dropna()
    p_animal_mwu = np.nan
    p_animal_t = np.nan
    if len(cci_an) >= 2 and len(ctl_an) >= 2:
        try:
            _, p_animal_mwu = stats.mannwhitneyu(cci_an, ctl_an, alternative="two-sided")
        except ValueError:
            pass
        _, p_animal_t = stats.ttest_ind(cci_an, ctl_an, equal_var=False)

    # (3) mixed model: RI ~ Condition + Side, random intercept per animal
    p_lmm = np.nan
    coef_lmm = np.nan
    lmm_note = ""
    try:
        sub_m = sub.copy()
        sub_m["Condition"] = pd.Categorical(sub_m["Condition"], categories=["Control", "CCI"])
        # include Side only when both sides are present
        formula = "ShollRamificationIndex ~ C(Condition)"
        if sub_m["Side"].nunique() > 1:
            formula += " + C(Side)"
        md = smf.mixedlm(formula, data=sub_m, groups=sub_m["AnimalID"])
        mf = md.fit(reml=True, method="lbfgs")
        term = [t for t in mf.params.index if "Condition" in t]
        if term:
            coef_lmm = float(mf.params[term[0]])
            p_lmm = float(mf.pvalues[term[0]])
    except Exception as e:  # noqa
        lmm_note = f"LMM error: {type(e).__name__}"

    rows.append({
        "Phase": phase,
        "n_animals_Control": int(ctl_an.shape[0]),
        "n_animals_CCI": int(cci_an.shape[0]),
        "n_cells_Control": int(len(ctl_cells)),
        "n_cells_CCI": int(len(cci_cells)),
        "median_RI_Control_animal": float(np.median(ctl_an)) if len(ctl_an) else np.nan,
        "median_RI_CCI_animal": float(np.median(cci_an)) if len(cci_an) else np.nan,
        "p_cell_MWU_pseudoreplicated": p_cell,
        "p_animal_MWU": p_animal_mwu,
        "p_animal_Welch_t": p_animal_t,
        "lmm_coef_CCI_vs_Control": coef_lmm,
        "p_lmm_animal_random_effect": p_lmm,
        "lmm_note": lmm_note,
    })

res = pd.DataFrame(rows)
res["_o"] = pd.Categorical(res["Phase"], categories=PHASE_ORDER, ordered=True)
res = res.sort_values("_o").drop(columns="_o")

out_csv = OUT_DIR / "cam_ramification_index_layer12_animal_and_mixedmodel.csv"
res.to_csv(out_csv, index=False)
animal.to_csv(OUT_DIR / "cam_ramification_index_layer12_animal_means.csv", index=False)

pd.set_option("display.width", 200, "display.max_columns", 30)
print("\n=== Ramification index (layer 12): cell vs animal vs modelo mixto ===\n")
print(res[["Phase", "n_animals_Control", "n_animals_CCI",
           "p_cell_MWU_pseudoreplicated", "p_animal_MWU",
           "p_animal_Welch_t", "p_lmm_animal_random_effect"]].to_string(index=False))

# --- Diagnostic figure: cells + per-animal mean ---
fig, axes = plt.subplots(1, len(phases), figsize=(3.0 * len(phases), 4.2), sharey=True)
if len(phases) == 1:
    axes = [axes]
order = ["Control", "CCI"]
for ax, phase in zip(axes, phases):
    sub = df[df["Phase"] == phase]
    a = animal[animal["Phase"] == phase]
    sns.stripplot(data=sub, x="Condition", y=RI, order=order, ax=ax,
                  color="0.7", size=2.5, jitter=0.18, alpha=0.5)
    sns.pointplot(data=a, x="Condition", y="RI_animal_mean", order=order, ax=ax,
                  color="black", linestyle="none", markers="D",
                  errorbar=("ci", 95), capsize=0.15)
    rrow = res[res["Phase"] == phase].iloc[0]
    ax.set_title(f"{phase}\nanimal n={int(rrow.n_animals_Control)}/{int(rrow.n_animals_CCI)}\n"
                 f"cell p={rrow.p_cell_MWU_pseudoreplicated:.4f}  LMM p={rrow.p_lmm_animal_random_effect:.3f}",
                 fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("Ramification index" if ax is axes[0] else "")
fig.suptitle("Ramification index (layer 12): grey = cells, black diamonds = per-animal mean (95% CI)",
             y=1.03, fontsize=11)
plt.tight_layout()
fig_path = OUT_DIR / "cam_ramification_index_layer12_animal_level.png"
plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.close("all")
print(f"\n[OK] Table: {out_csv}")
print(f"[OK] Figure: {fig_path}")
