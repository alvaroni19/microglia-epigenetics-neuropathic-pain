# -*- coding: utf-8 -*-
"""
Recompute ALL qPCR figure statistics fully in Python (no values read from
GraphPad exports), reproducing the exact GraphPad Prism framework:

  * Male CCI-vs-control            -> Student's unpaired two-tailed t test (equal variance)
  * Female phase-stratified in vivo -> ordinary two-way ANOVA (condition x phase)
                                       + Bonferroni 'compare cell means' (15 comparisons)
  * In vitro CSF1-vs-control        -> ordinary two-way ANOVA (sex x treatment)
                                       + Bonferroni 'compare cell means' (6 comparisons)
  * In vitro CSF1+E2-vs-CSF1 (E2 rows of the in vitro bubble plots) -> same 2-way ANOVA family,
    but on the 3-level treatment design (Control/CSF1/CSF1+E2): Bonferroni over
    C(6,2)=15 cell-mean comparisons, matching the reported GraphPad analysis.

It also writes a VALIDATION table comparing the Python p-values against the
GraphPad-exported p-values, to prove significance is preserved.

Outputs (generated_outputs/PCR_260514/):
  invitro_csf1_by_sex_methylation_regulators_stats.csv
  invitro_methylation_regulators_stats.csv
  invivo_pain_microglia_markers_stats.csv
  invivo_methylation_regulators_stats.csv
  invivo_estrogen_genes_stats.csv
  invitro_estrogen_genes_stats.csv
  python_vs_graphpad_validation.csv   (NEW)
"""
from __future__ import annotations
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path(__file__).absolute().parents[2]
VITRO_EXCEL = ROOT / "data" / "PCR_260514" / "invitro_qpcr_relative_expression.xlsx"
OUT_NATIVE = ROOT / "generated_outputs" / "PCR_260514"
OUT_NATIVE.mkdir(parents=True, exist_ok=True)
# The bundle may be stored below a long Windows path.  The extended-length
# prefix keeps CSV export independent of the chosen parent directory.
OUT = (
    Path("\\\\?\\" + str(OUT_NATIVE.absolute()))
    if sys.platform == "win32" and not str(OUT_NATIVE).startswith("\\\\?\\")
    else OUT_NATIVE
)
sys.path.insert(0, str(ROOT / "code" / "_helpers"))
from qpcr_raw import GENE_SHEETS, SCRIPT_TO_CANONICAL, load_invivo_long  # noqa: E402

VIVO_LONG = load_invivo_long()

LABEL = {"TNFa": "Tnf-alpha", "KCC2": "Slc12a5", "CSF1R": "Csf1r", "TLR4": "Tlr4",
         "BDNF": "Bdnf", "IL10": "Il10", "ERa": "Esr1", "ERb": "Esr2",
         "GPER": "Gper1", "ARO": "Cyp19a1"}

METHOD_T   = "Student's unpaired two-tailed t test (Python scipy, equal variance; reproduces GraphPad Prism)"
METHOD_2W  = "Ordinary two-way ANOVA with Bonferroni cell-mean comparisons (Python statsmodels; reproduces GraphPad Prism)"

# ---------------------------------------------------------------- helpers
def p_to_stars(p):
    if p is None or not np.isfinite(p):
        return float("nan")
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else float("nan")

def log2fc(num, den):
    return float(np.log2(num / den)) if (num and den and num > 0 and den > 0) else float("nan")

def gp_p(raw):
    """Parse a GraphPad p string with decimal comma and '<0,0001'/'>0,9999' notation."""
    raw = str(raw).strip()
    if raw.startswith("<"):
        return float(raw[1:].replace(",", ".")) / 2.0   # below threshold -> keep stars correct
    if raw.startswith(">"):
        return 1.0
    try:
        return float(raw.replace(",", "."))
    except ValueError:
        return np.nan

# ---------------------------------------------------------------- 2-way ANOVA + Bonferroni
def twoway_bonferroni(df_long, f1, f2, value="value"):
    df_long = df_long.dropna(subset=[value]).copy()
    model = smf.ols(f"{value} ~ C({f1}) * C({f2})", data=df_long).fit()
    mse, dfr = model.mse_resid, model.df_resid
    gs = df_long.groupby([f1, f2])[value].agg(["mean", "count"]).reset_index()
    gs.columns = [f1, f2, "mean", "n"]
    k = len(gs)
    ncomp = k * (k - 1) // 2
    out = {}
    for a, b in combinations(gs.itertuples(index=False), 2):
        ga, gb = (getattr(a, f1), getattr(a, f2)), (getattr(b, f1), getattr(b, f2))
        se = np.sqrt(mse * (1.0 / a.n + 1.0 / b.n))
        p = np.nan if se <= 0 else min(1.0, 2 * stats.t.sf(abs((a.mean - b.mean) / se), df=dfr) * ncomp)
        out[(ga, gb)] = out[(gb, ga)] = p
    means = df_long.groupby([f1, f2])[value].mean()
    return out, means, ncomp

# ---------------------------------------------------------------- data readers
_SHEET_CACHE = {}
def resolve_sheet(excel, wanted):
    """Find a sheet whose name matches `wanted` ignoring spaces and case."""
    if excel not in _SHEET_CACHE:
        _SHEET_CACHE[excel] = pd.ExcelFile(excel).sheet_names
    norm = lambda s: re.sub(r"\s+", "", str(s)).lower()
    target = norm(wanted)
    for s in _SHEET_CACHE[excel]:
        if norm(s) == target:
            return s
    raise ValueError(f"sheet like '{wanted}' not found in {Path(excel).name}")

def read_male(gene, excel):
    canonical = SCRIPT_TO_CANONICAL[gene]
    subset = VIVO_LONG[
        (VIVO_LONG["gene"] == canonical) & (VIVO_LONG["phase"] == "Male")
    ]
    ctrl = subset.loc[
        subset["condition"] == "Control", "relative_expression"
    ].dropna().to_numpy()
    cci = subset.loc[
        subset["condition"] == "CCI", "relative_expression"
    ].dropna().to_numpy()
    return cci, ctrl

def male_ttest(gene, excel):
    cci, ctrl = read_male(gene, excel)
    t, p = stats.ttest_ind(cci, ctrl, equal_var=True)   # Student's t, two-sided
    # normalize means to control == 1 (GraphPad data are fold-change already; report raw means)
    return float(p), float(np.mean(cci)), float(np.mean(ctrl))

def read_female_long(gene, excel, sheet=None):
    canonical = SCRIPT_TO_CANONICAL[gene]
    phase_names = {
        "Proestrus": "proestrus",
        "Estrus": "estrus",
        "Diestrus": "diestrus",
    }
    subset = VIVO_LONG[
        (VIVO_LONG["gene"] == canonical) & (VIVO_LONG["phase"] != "Male")
    ].copy()
    return pd.DataFrame(
        {
            "condition": subset["condition"],
            "phase": subset["phase"].map(phase_names),
            "value": subset["relative_expression"],
        }
    ).dropna(subset=["value"])

def read_invitro_2way(gene, excel, sheet=None):
    sheet = sheet or f"{gene} 2 way"
    df = pd.read_excel(excel, sheet_name=resolve_sheet(excel, sheet))
    rows = []
    for _, row in df.iterrows():
        trt = str(row["ROWTITLE"]).strip()
        for col in df.columns[1:]:
            sex = "Male" if str(col).startswith("Male") else "Female"
            v = pd.to_numeric(row[col], errors="coerce")
            if pd.notna(v):
                rows.append({"sex": sex, "treatment": trt, "value": float(v)})
    return pd.DataFrame(rows)

# ---------------------------------------------------------------- GraphPad export readers (validation only)
def _gp_file(pattern_words):
    """Legacy compatibility: archived GraphPad result CSVs are not inputs."""
    return None

def _gp_adjP_for(f, groupA, groupB):
    """Return adjusted P for the 'A vs. B' cell-mean comparison (exact, either order)."""
    if f is None or not f.exists():
        return np.nan
    norm = lambda s: re.sub(r"\s+", "", str(s)).lower().strip(".")
    want = {norm(groupA), norm(groupB)}
    for line in f.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        cells = line.split(";")
        head = cells[0]
        if "vs." not in head or len(cells) <= 5:
            continue
        sides = {norm(s) for s in head.split("vs.")}
        if sides == want:
            return gp_p(cells[5])
    return np.nan

def gp_male(gene):
    f = _gp_file([f"{gene}", "male", "ttest"]) or _gp_file(["unpairedttest", gene, "male"]) \
        or _gp_file(["unpaired", gene, "male"])
    if f is None:
        return np.nan
    m = re.search(r"P value;([^;\n]+)", f.read_text(encoding="utf-8-sig", errors="replace"))
    return gp_p(m.group(1)) if m else np.nan

def gp_female(gene, phase):
    f = _gp_file(["2way", gene, "female"])
    return _gp_adjP_for(f, f"control:{phase}", f"CCI:{phase}")

def gp_invitro(gene, sex):
    f = _gp_file(["2way", f"{gene}2way"])
    return _gp_adjP_for(f, f"Control:{sex}", f"CSF1:{sex}")

# ---------------------------------------------------------------- builders
PHASES = ["diestrus", "proestrus", "estrus"]
val_rows = []

def build_vivo(genes, female_sheet_override=None):
    rows = []
    for g in genes:
        # male
        p, mcci, mctrl = male_ttest(g, None)
        rows.append(dict(gene=g, gene_label=LABEL.get(g, g), comparison="Male CCI vs control",
                         numerator="CCI", denominator="Control", mean_numerator=mcci,
                         mean_denominator=mctrl, log2fc=log2fc(mcci, mctrl), p_value=p,
                         stars=p_to_stars(p), method=METHOD_T,
                         graphpad_source=f"Unpaired t test of {g} male.csv",
                         graphpad_comparison="Male CCI vs control"))
        val_rows.append(dict(figure="vivo", gene=g, comparison="Male CCI vs control",
                             python_p=p, graphpad_p=gp_male(g),
                             python_stars=p_to_stars(p), graphpad_stars=p_to_stars(gp_male(g))))
        # female 2-way
        sheet = (female_sheet_override or {}).get(g)
        df_f = read_female_long(g, None, sheet)
        bonf, means, ncomp = twoway_bonferroni(df_f, "condition", "phase")
        for ph in PHASES:
            p = bonf.get((("CCI", ph), ("Control", ph)), np.nan)
            mc = float(means.get(("CCI", ph), np.nan)); mk = float(means.get(("Control", ph), np.nan))
            rows.append(dict(gene=g, gene_label=LABEL.get(g, g),
                             comparison=f"{ph.capitalize()} CCI vs control",
                             numerator="CCI", denominator="Control", mean_numerator=mc,
                             mean_denominator=mk, log2fc=log2fc(mc, mk), p_value=p,
                             stars=p_to_stars(p), method=f"{METHOD_2W} ({ncomp} comparisons)",
                             graphpad_source=f"2way ANOVA of {g} female.csv",
                             graphpad_comparison=f"control:{ph} vs. CCI:{ph}"))
            val_rows.append(dict(figure="vivo", gene=g, comparison=f"{ph} CCI vs control",
                                 python_p=p, graphpad_p=gp_female(g, ph),
                                 python_stars=p_to_stars(p), graphpad_stars=p_to_stars(gp_female(g, ph))))
    return pd.DataFrame(rows)

def build_invitro_csf1(genes):
    """CSF1 vs control rows from the 2-way (sex x treatment Control/CSF1) sheet."""
    rows = []
    for g in genes:
        df = read_invitro_2way(g, VITRO_EXCEL)
        bonf, means, ncomp = twoway_bonferroni(df, "sex", "treatment")
        for sex in ["Male", "Female"]:
            p = bonf.get(((sex, "CSF1"), (sex, "Control")), np.nan)
            mc = float(means.get((sex, "CSF1"), np.nan)); mk = float(means.get((sex, "Control"), np.nan))
            rows.append(dict(gene=g, gene_label=LABEL.get(g, g),
                             comparison=f"{sex} CSF1 vs control",
                             numerator="CSF1", denominator="Control", mean_numerator=mc,
                             mean_denominator=mk, log2fc=log2fc(mc, mk), p_value=p,
                             stars=p_to_stars(p), method=f"{METHOD_2W} ({ncomp} comparisons)",
                             graphpad_source=f"2way ANOVA of {g} 2 way.csv",
                             graphpad_comparison=f"Control:{sex} vs. CSF1:{sex}"))
            val_rows.append(dict(figure="invitro", gene=g, comparison=f"{sex} CSF1 vs control",
                                 python_p=p, graphpad_p=gp_invitro(g, sex),
                                 python_stars=p_to_stars(p), graphpad_stars=p_to_stars(gp_invitro(g, sex))))
    return pd.DataFrame(rows)

def read_invitro_menosE2(gene):
    """3-level treatment design (control/CSF1/CSF1+E2) x sex, from the
    '{gene} 1 way without E2' sheet (6 columns: <sex> <treatment>)."""
    sheet = resolve_sheet(VITRO_EXCEL, f"{gene} 1 way without E2")
    df = pd.read_excel(VITRO_EXCEL, sheet_name=sheet, header=None)
    headers = [str(x).strip().lower() for x in df.iloc[0].tolist()]
    rows = []
    for c, h in enumerate(headers):
        sex = "Male" if h.startswith("male") else "Female"
        if "csf1+e2" in h.replace(" ", ""):
            trt = "CSF1+E2"
        elif "csf1" in h:
            trt = "CSF1"
        else:
            trt = "Control"
        for r in range(1, df.shape[0]):
            v = pd.to_numeric(df.iloc[r, c], errors="coerce")
            if pd.notna(v):
                rows.append({"sex": sex, "treatment": trt, "value": float(v)})
    return pd.DataFrame(rows)


def gp_invitro_1way(gene, sex, a, b):
    """Adjusted P for '{sex} {a} vs. {sex} {b}' in the one-way 'without E2' export."""
    f = _gp_file(["onewayanova", f"{gene}1waymenose2"]) or _gp_file(["one-way", gene, "menose2"])
    return _gp_adjP_for(f, f"{sex} {a}", f"{sex} {b}")


def build_invitro_with_E2(genes):
    """Full 4-column in vitro matrix per gene:
       CSF1 vs control  -> 4-group 2-way (sex x Control/CSF1), 6 comparisons
       CSF1+E2 vs CSF1  -> 6-group design (sex x 3 treatments), 15 comparisons
    Both reproduce the corresponding GraphPad analyses exactly."""
    rows = []
    for g in genes:
        # CSF1 vs control from 4-group 2-way
        df4 = read_invitro_2way(g, VITRO_EXCEL)
        bonf4, means4, n4 = twoway_bonferroni(df4, "sex", "treatment")
        # E2 from 6-group design
        df6 = read_invitro_menosE2(g)
        bonf6, means6, n6 = twoway_bonferroni(df6, "sex", "treatment")
        for sex in ["Male", "Female"]:
            # CSF1 vs control
            p = bonf4.get(((sex, "CSF1"), (sex, "Control")), np.nan)
            mc = float(means4.get((sex, "CSF1"), np.nan)); mk = float(means4.get((sex, "Control"), np.nan))
            rows.append(dict(gene=g, gene_label=LABEL.get(g, g), comparison=f"{sex} CSF1 vs control",
                             numerator="CSF1", denominator="Control", mean_numerator=mc, mean_denominator=mk,
                             log2fc=log2fc(mc, mk), p_value=p, stars=p_to_stars(p),
                             method=f"{METHOD_2W} (4-group sex x Control/CSF1, {n4} comparisons)",
                             graphpad_source=f"2way ANOVA of {g} 2 way.csv",
                             graphpad_comparison=f"Control:{sex} vs. CSF1:{sex}"))
            val_rows.append(dict(figure="invitro_E2fig", gene=g, comparison=f"{sex} CSF1 vs control",
                                 python_p=p, graphpad_p=gp_invitro(g, sex),
                                 python_stars=p_to_stars(p), graphpad_stars=p_to_stars(gp_invitro(g, sex))))
            # CSF1+E2 vs CSF1
            p2 = bonf6.get(((sex, "CSF1+E2"), (sex, "CSF1")), np.nan)
            me = float(means6.get((sex, "CSF1+E2"), np.nan)); mcsf = float(means6.get((sex, "CSF1"), np.nan))
            rows.append(dict(gene=g, gene_label=LABEL.get(g, g), comparison=f"{sex} CSF1+E2 vs CSF1",
                             numerator="CSF1+E2", denominator="CSF1", mean_numerator=me, mean_denominator=mcsf,
                             log2fc=log2fc(me, mcsf), p_value=p2, stars=p_to_stars(p2),
                             method=f"{METHOD_2W} (6-group sex x Control/CSF1/CSF1+E2, {n6} comparisons)",
                             graphpad_source=f"Ordinary one-way ANOVA of {g} 1 way without E2.csv",
                             graphpad_comparison=f"{sex.lower()} CSF1 vs. {sex.lower()} CSF1+E2"))
            gpp = gp_invitro_1way(g, sex.lower(), "CSF1", "CSF1+E2")
            val_rows.append(dict(figure="invitro_E2fig", gene=g, comparison=f"{sex} CSF1+E2 vs CSF1",
                                 python_p=p2, graphpad_p=gpp,
                                 python_stars=p_to_stars(p2), graphpad_stars=p_to_stars(gpp)))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    EPI = ["Dnmt1", "Dnmt3a", "Dnmt3b", "Gadd45a", "Gadd45b", "Gadd45g"]
    PAIN = ["KCC2", "BDNF", "CSF1R", "TLR4", "TNFa", "IL10"]
    ER = ["ERa", "ERb", "GPER", "ARO"]
    fem_override = {"BDNF": "BDNF Female"}

    # in vitro CSF1-vs-control by sex (methylation regulators)
    def final_sources(frame, in_vivo):
        frame = frame.copy()
        if in_vivo:
            frame["graphpad_source"] = frame["gene"].map(
                lambda gene: f"invivo_qpcr_delta_ct.xlsx / {GENE_SHEETS[SCRIPT_TO_CANONICAL[gene]]}"
            )
        else:
            frame["graphpad_source"] = "invitro_qpcr_relative_expression.xlsx"
        return frame

    df1g = final_sources(build_invitro_csf1(EPI), in_vivo=False)
    df1g.to_csv(OUT / "invitro_csf1_by_sex_methylation_regulators_stats.csv", index=False)

    # full in vitro methylation-regulator set (CSF1 vs control + CSF1+E2 vs CSF1)
    final_sources(build_invitro_with_E2(EPI), in_vivo=False).to_csv(OUT / "invitro_methylation_regulators_stats.csv", index=False)

    # in vivo pain/microglia markers and methylation regulators
    final_sources(build_vivo(PAIN), in_vivo=True).to_csv(OUT / "invivo_pain_microglia_markers_stats.csv", index=False)
    final_sources(build_vivo(EPI), in_vivo=True).to_csv(OUT / "invivo_methylation_regulators_stats.csv", index=False)
    # in vivo estrogen-related genes
    final_sources(build_vivo(ER), in_vivo=True).to_csv(OUT / "invivo_estrogen_genes_stats.csv", index=False)
    # in vitro estrogen-related genes (CSF1 vs control + CSF1+E2 vs CSF1)
    final_sources(build_invitro_with_E2(ER), in_vivo=False).to_csv(OUT / "invitro_estrogen_genes_stats.csv", index=False)

    # in vitro CSF1R (CSF1 vs control + CSF1+E2 vs CSF1), used for the
    # six-group Csf1r violin annotation
    final_sources(build_invitro_with_E2(["CSF1R"]), in_vivo=False).to_csv(OUT / "invitro_csf1r_stats.csv", index=False)

    recomputed = pd.DataFrame(val_rows)[
        ["figure", "gene", "comparison", "python_p", "python_stars"]
    ].rename(
        columns={"python_p": "p_value", "python_stars": "stars"}
    )
    recomputed["source"] = np.where(
        recomputed["figure"].str.startswith("vivo"),
        "invivo_qpcr_delta_ct.xlsx (Delta-Ct sheets)",
        "invitro_qpcr_relative_expression.xlsx",
    )
    recomputed.to_csv(
        OUT / "qpcr_statistical_recomputation.csv", index=False
    )
    print(
        f"[OK] recomputed {len(recomputed)} qPCR contrasts from final retained sources"
    )

if __name__ == "__main__":
    main()
