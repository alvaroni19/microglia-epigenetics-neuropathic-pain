#!/usr/bin/env python3
"""Build the v2 Zenodo deposit, the GitHub repository content and the Zenodo ZIP.

Run from the repository root of the working master:
    python code/_maintenance/build_public_release_v2.py <output_dir>
"""
from __future__ import annotations

import hashlib
import shutil
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).absolute().parents[2]

CODE_FILES = {
    "MeDIP": ["medips_differential_methylation.Rmd", "rn6_annotation.R", "figure_panels.py",
              "build_csf1r_locus_track.py", "build_medip_robpca.R", "build_medip_upset.R",
              "build_medip_global_methylation.R", "validate_medip_outputs.py", "OUTPUTS.md",
              "REFERENCE_DATA.md"],
    "MorphOMICs": ["build_sholl_report_from_swc.py", "cam_morphomics_utils.py", "morphomics_plot_utils.py",
                   "plot_persistence_diff_recolored.py", "plot_ramification_index_layer12_LMM_english.py",
                   "ramification_index_layer12_animal_and_mixedmodel.py", "build_ramification_pain_source.py",
                   "plot_ramification_pain_correlation_layer12_english.py", "plot_umap_layer12_full_and_zoom.py",
                   "plot_umap_layer12_cycles_publication_recolored.py"],
    "PCR_260514": ["build_invivo_sources.py", "recompute_all_python_stats.py", "export_invivo_contrast_statistics.py",
                   "figure_panels.py", "build_invitro_2way_violins.py", "build_invitro_6group_violins.py",
                   "build_qpcr_male_violins.py", "build_qpcr_female_violins.py", "build_invitro_estrogenic_score.R"],
    "PCR_MFA": ["build_invivo_mfa.R", "build_invivo_mfa_plot.py", "figure_panels.py"],
    "Behaviour": ["behaviour_reanalysis.R", "build_behaviour_violins.py"],
    "Figures_v2": ["build_manuscript_panels_v2.py", "original_style.py", "build_panels_v2_mfa_morphology.py",
                   "build_medip_overview_panels.R", "figure_layout.py", "build_figure_blocks.py"],
    "_helpers": ["journal_style.py", "violin_style.py", "map_csf1r_esr_motifs_integrated.py", "qpcr_raw.py",
                 "pcr_python_stats.py", "qpcr_violin_style.py", "invitro_violin_style.py", "unified_palette.py",
                 "publication_style.R", "final_figure_controls.R"],
    "_maintenance": ["build_public_release_v2.py"],
}
CODE_GLOBS = {}
CODE_TREES = ["ramification_index", "imaris_swc_export"]
ROOT_CODE = ["reproduce_all_figures.py"]

MEDIP_TABLES = ["DMR_candidate_genes_ERE_scan_ws400.xlsx", "GENES_TOTAL_sig_contig_colap_ws400.xlsx",
                "GO_relevance_rank_distribution_ws400.xlsx", "csf1r_context_rn6.fasta",
                "csf1r_dmr_windows_rn6.fasta", "heatmap_candidate_genes_logFC_ws400_orders.xlsx",
                "medip_sample_metadata.csv", "methylation_map_candidate_genes_hits_ws400.xlsx",
                "significant_gene_rms_matrix_ws400.csv"]
DATA_TREES = ["MeDIP/annotations", "MeDIP/motifs", "MorphOMICs/reports_cam", "MorphOMICs/saved_instances_cam",
              "SWC_reconstructions"]
DATA_FILES = ["MeDIP/MEDIP_PROVIDER_QC_METRICS.tsv"]
OUTPUT_SUFFIXES = {".csv", ".xlsx", ".svg", ".pdf"}
OUTPUT_DIRS = ["MeDIP", "MorphOMICs", "Figures_v2"]
# Behaviour and RT-qPCR inputs and everything derived from them are not distributed.
EXCLUDED_OUTPUT_WORDS = ("pain", "allodynia")


def public_output(rel: str) -> bool:
    """Distributed outputs: MeDIP and morphology analyses without behavioural data."""
    name = rel.rsplit("/", 1)[-1].lower()
    if any(w in name for w in EXCLUDED_OUTPUT_WORDS):
        return False
    if rel.startswith("Figures_v2/"):
        return name.startswith(("medip_", "morphology_"))
    return True
DOCS = ["README.md", "CHANGELOG.md", "CITATION.cff", "LICENSE", "renv.lock", "requirements-python.txt",
        "python-version.txt", "ENVIRONMENT_v2.md"]
SKIP_NAMES = {"Thumbs.db", "desktop.ini", ".DS_Store"}


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_tree(src: Path, dst: Path, suffixes=None, keep=None) -> None:
    for item in src.rglob("*"):
        if item.is_file() and item.name not in SKIP_NAMES and "__pycache__" not in item.parts \
                and not item.name.startswith("~$") and "legacy_both_sides" not in item.name:
            if (suffixes is None or item.suffix.lower() in suffixes) and (keep is None or keep(item)):
                copy(item, dst / item.relative_to(src))


def build_code(target: Path) -> None:
    for folder, names in CODE_FILES.items():
        for name in names:
            copy(ROOT / "code" / folder / name, target / "code" / folder / name)
    for folder, patterns in CODE_GLOBS.items():
        for pattern in patterns:
            for f in (ROOT / "code" / folder).glob(pattern):
                copy(f, target / "code" / folder / f.name)
    for tree in CODE_TREES:
        copy_tree(ROOT / "code" / tree, target / "code" / tree)
    for name in ROOT_CODE:
        copy(ROOT / "code" / name, target / "code" / name)


def build_docs(target: Path, docs_dir: Path) -> None:
    for name in DOCS:
        src = docs_dir / name
        if src.exists():
            copy(src, target / name)


def build_data(target: Path) -> None:
    for name in MEDIP_TABLES:
        copy(ROOT / "data" / "MeDIP" / "workspace_and_tables" / name,
             target / "data" / "MeDIP" / "workspace_and_tables" / name)
    for tree in DATA_TREES:
        copy_tree(ROOT / "data" / tree, target / "data" / tree)
    for name in DATA_FILES:
        copy(ROOT / "data" / name, target / "data" / name)
    gen = ROOT / "generated_outputs"
    for folder in OUTPUT_DIRS:
        copy_tree(gen / folder, target / "generated_outputs" / folder, OUTPUT_SUFFIXES,
                  keep=lambda f: public_output(f.relative_to(gen).as_posix()))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 23), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(target: Path) -> None:
    rows = ["path\tbytes\tsha256"]
    for f in sorted(p for p in target.rglob("*") if p.is_file()):
        rel = f.relative_to(target).as_posix()
        if rel == "metadata/file_manifest_sha256.tsv":
            continue
        rows.append(f"{rel}\t{f.stat().st_size}\t{sha256(f)}")
    out = target / "metadata" / "file_manifest_sha256.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")


def make_zip(src: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for f in sorted(p for p in src.rglob("*") if p.is_file()):
            zf.write(f, f.relative_to(src).as_posix())


def main() -> None:
    out = Path(sys.argv[1]).absolute()
    docs = ROOT / "release_docs_v2"
    zen, gh = out / "zenodo_deposit", out / "github_repository"
    for d in (zen, gh):
        if d.exists():
            shutil.rmtree(d)
    build_code(zen); build_docs(zen, docs); build_data(zen)
    copy(docs / "zenodo_metadata.json", zen / "zenodo_metadata.json")
    write_manifest(zen)
    build_code(gh); build_docs(gh, docs)
    copy(docs / ".gitignore", gh / ".gitignore")
    copy(docs / ".gitattributes", gh / ".gitattributes")
    zip_path = out / "microglia-epigenetics-neuropathic-pain_v2.zip"
    make_zip(zen, zip_path)
    (out / (zip_path.name + ".sha256")).write_text(f"{sha256(zip_path)}  {zip_path.name}\n")
    print("[OK]", zen, gh, zip_path)


if __name__ == "__main__":
    main()
