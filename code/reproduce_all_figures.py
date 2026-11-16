#!/usr/bin/env python3
"""Regenerate every quantitative panel from the inputs that are present.

Domains: medip, morphology (distributed inputs) and qpcr, behaviour,
multivariate (behaviour and RT-qPCR inputs, not distributed; available from the
corresponding author on request). Steps whose inputs are absent are reported as
[SKIP] and the remaining analyses still run.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def native_path(path: Path) -> Path:
    text = str(path.absolute())
    if text.startswith("\\\\?\\UNC\\"):
        text = "\\\\" + text[8:]
    elif text.startswith("\\\\?\\"):
        text = text[4:]
    return Path(text)


ROOT = native_path(Path(__file__)).parents[1]
SCRIPTS = ROOT / "code"
DATA = ROOT / "data"

QPCR_VITRO = DATA / "PCR_260514" / "invitro_qpcr_relative_expression.xlsx"
QPCR_VIVO = DATA / "PCR_260514" / "invivo_qpcr_delta_ct.xlsx"
BEHAVIOUR = DATA / "Allodynia_0_16" / "allodynia_change_scores.xlsx"
VON_FREY = DATA / "Allodynia_0_16" / "allodynia_von_frey_scores.xlsx"
LINKAGE = DATA / "MorphOMICs" / "morphology_behaviour_linkage.xlsx"

DOMAIN_INPUTS = {
    "medip": [],
    "morphology": [],
    "qpcr": [QPCR_VITRO, QPCR_VIVO, VON_FREY],
    "behaviour": [BEHAVIOUR],
    "multivariate": [QPCR_VIVO, VON_FREY],
    "figures_v2": [],          # panel builders skip unavailable panels themselves
}
TASK_INPUTS = {
    "MorphOMICs/build_ramification_pain_source.py": [VON_FREY, LINKAGE],
    "MorphOMICs/plot_ramification_pain_correlation_layer12_english.py": [VON_FREY, LINKAGE],
}

PYTHON_TASKS = {
    "medip": [
        "_helpers/map_csf1r_esr_motifs_integrated.py",
        "MeDIP/figure_panels.py",
        "MeDIP/build_csf1r_locus_track.py",
    ],
    "qpcr": [
        "PCR_260514/recompute_all_python_stats.py",
        "PCR_260514/export_invivo_contrast_statistics.py",
        "PCR_260514/build_invitro_2way_violins.py",
        "PCR_260514/build_invitro_6group_violins.py",
        "PCR_260514/build_qpcr_male_violins.py",
        "PCR_260514/build_qpcr_female_violins.py",
        "PCR_260514/figure_panels.py",
    ],
    "behaviour": [
        "Behaviour/build_behaviour_violins.py",
    ],
    "multivariate": [
        "PCR_MFA/build_invivo_mfa_plot.py",
        "PCR_MFA/figure_panels.py",
    ],
    "morphology": [
        "MorphOMICs/build_sholl_report_from_swc.py",
        "MorphOMICs/plot_ramification_index_layer12_LMM_english.py",
        "MorphOMICs/plot_umap_layer12_full_and_zoom.py",
        "MorphOMICs/plot_persistence_diff_recolored.py",
        "MorphOMICs/build_ramification_pain_source.py",
        "MorphOMICs/plot_ramification_pain_correlation_layer12_english.py",
    ],
    "figures_v2": [
        "Figures_v2/build_manuscript_panels_v2.py",
        "Figures_v2/build_figure_blocks.py",
    ],
}

# R steps that must run before the Python steps of the same domain
R_BEFORE = {
    "qpcr": ["PCR_260514/build_invitro_estrogenic_score.R"],
    "behaviour": ["Behaviour/behaviour_reanalysis.R"],
    "multivariate": ["PCR_MFA/build_invivo_mfa.R"],
}
# R steps that run after all Python steps (they receive the repository root)
R_AFTER = {
    "medip": [
        "MeDIP/build_medip_robpca.R",
        "MeDIP/build_medip_upset.R",
        "MeDIP/build_medip_global_methylation.R",
    ],
    # runs after Figures_v2/build_manuscript_panels_v2.py, which recreates the folder
    "figures_v2": ["Figures_v2/build_medip_overview_panels.R"],
}
RESTRICTED = "behaviour and RT-qPCR inputs are not distributed (available from the corresponding author on request)"


def run(command: list[str]) -> None:
    print("+", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def find_rscript() -> str:
    for candidate in (os.environ.get("RSCRIPT"), shutil.which("Rscript")):
        if candidate and Path(candidate).exists():
            return str(candidate)
    raise FileNotFoundError("Rscript (R 4.4) was not found on PATH; set RSCRIPT to its location.")


def missing(paths: list[Path]) -> list[Path]:
    return [p for p in paths if not p.exists()]


def skip(label: str, absent: list[Path]) -> None:
    print(f"[SKIP] {label}: {absent[0].relative_to(ROOT).as_posix()} not found; {RESTRICTED}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain", action="append", choices=sorted(PYTHON_TASKS),
                        help="Regenerate only selected domains; may be repeated.")
    parser.add_argument("--skip-r", action="store_true", help="Skip R steps and run Python steps only.")
    args = parser.parse_args()
    requested = args.domain or list(PYTHON_TASKS)
    domains = []
    for domain in requested:
        absent = missing(DOMAIN_INPUTS[domain])
        if absent:
            skip(f"domain '{domain}'", absent)
        else:
            domains.append(domain)

    # qPCR panels and multivariate analyses share one animal-level table built
    # from the final Delta-Ct values and the von Frey scores.
    if {"qpcr", "multivariate"}.intersection(domains):
        run([sys.executable, str(SCRIPTS / "PCR_260514" / "build_invivo_sources.py")])

    if not args.skip_r:
        for domain in domains:
            for relative in R_BEFORE.get(domain, []):
                run([find_rscript(), str(SCRIPTS / relative)])

    for domain in domains:
        for relative in PYTHON_TASKS[domain]:
            absent = missing(TASK_INPUTS.get(relative, []))
            if absent:
                skip(relative, absent)
                continue
            run([sys.executable, str(SCRIPTS / relative)])

    if not args.skip_r:
        for domain in domains:
            for relative in R_AFTER.get(domain, []):
                run([find_rscript(), str(SCRIPTS / relative), str(ROOT)])

    print("[OK] Requested analyses completed.")


if __name__ == "__main__":
    main()
