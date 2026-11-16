"""Export t statistics, degrees of freedom and P values of the in vivo RT-qPCR
contrasts (the values reported in the figure legends).

Same tests as recompute_all_python_stats.py: males, Student's two-sided
unpaired t test; females, condition x phase two-way ANOVA with Bonferroni
cell-mean comparisons (15 comparisons, pooled residual variance).
Writes generated_outputs/PCR_260514/invivo_contrast_statistics_for_legends.csv
and checks that every P value matches the release bubble-statistics tables.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
from qpcr_raw import load_invivo_long  # noqa: E402

OUT = ROOT / "generated_outputs" / "PCR_260514"
PHASES = ["Proestrus", "Estrus", "Diestrus"]


def main() -> None:
    long = load_invivo_long()
    rows = []
    for gene, sub in long.groupby("gene"):
        m = sub[sub.phase == "Male"].dropna(subset=["relative_expression"])
        cci = m[m.condition == "CCI"].relative_expression.to_numpy()
        ctl = m[m.condition == "Control"].relative_expression.to_numpy()
        t, p = stats.ttest_ind(cci, ctl, equal_var=True)
        rows.append(dict(gene=gene, comparison="Male CCI vs control", n_cci=len(cci), n_control=len(ctl),
                         t=t, df=len(cci) + len(ctl) - 2, p_raw=p, p_reported=p, adjustment="none"))
        f = sub[sub.phase != "Male"].dropna(subset=["relative_expression"])
        fit = smf.ols("relative_expression ~ C(condition) * C(phase)", data=f).fit()
        mse, dfr = fit.mse_resid, int(fit.df_resid)
        cells = f.groupby(["condition", "phase"]).relative_expression.agg(["mean", "count"])
        ncomp = len(cells) * (len(cells) - 1) // 2
        for ph in PHASES:
            a, b = cells.loc[("CCI", ph)], cells.loc[("Control", ph)]
            se = np.sqrt(mse * (1 / a["count"] + 1 / b["count"]))
            tt = (a["mean"] - b["mean"]) / se
            pr = 2 * stats.t.sf(abs(tt), dfr)
            rows.append(dict(gene=gene, comparison=f"{ph} CCI vs control", n_cci=int(a["count"]),
                             n_control=int(b["count"]), t=tt, df=dfr, p_raw=pr,
                             p_reported=min(1.0, pr * ncomp), adjustment=f"Bonferroni x{ncomp}"))
    out = pd.DataFrame(rows)
    ref = pd.concat([pd.read_csv(OUT / f) for f in (
        "invivo_pain_microglia_markers_stats.csv", "invivo_methylation_regulators_stats.csv",
        "invivo_estrogen_genes_stats.csv")], ignore_index=True)
    alias = {"KCC2": "Slc12a5", "BDNF": "Bdnf", "CSF1R": "Csf1r", "TLR4": "Tlr4", "TNFa": "Tnf", "IL10": "Il10",
             "ERa": "Esr1", "ERb": "Esr2", "GPER": "Gper1", "ARO": "Cyp19a1"}
    ref["gene"] = ref["gene"].map(lambda g: alias.get(g, g))
    chk = out.merge(ref[["gene", "comparison", "p_value"]], on=["gene", "comparison"], how="inner")
    assert np.allclose(chk.p_reported, chk.p_value, rtol=1e-6, atol=1e-12), "P mismatch with release tables"
    out.to_csv(OUT / "invivo_contrast_statistics_for_legends.csv", index=False)
    print(f"[OK] {len(out)} contrasts; {len(chk)} checked against release tables")


if __name__ == "__main__":
    main()
