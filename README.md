# CAM study analysis code

This repository contains the custom analysis code associated with the CAM study.
MeDIP-seq raw and processed sequencing files (FASTQ, WIG) are deposited in GEO
(Series accession [GSE338843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338843)). The MeDIP-seq analysis tables, MorphOMICs intermediate
data and SWC reconstructions are archived in the associated Zenodo record https://zenodo.org/records/21396206. The
aligned rn6 BAM files are **not** distributed because of their size; regenerate
them from the GEO FASTQ with HISAT2 v2.1.0 (see "Data layout") or request them
from the corresponding author.

## Included analyses

- `code/MeDIP/MEDIPS - Alvaro Bautista-Abad.Rmd`: MEDIPS/edgeR analysis using
  400-bp windows and export of the tables supporting the manuscript figures.
- `code/MeDIP/rn6_annotation.R`: local rn6 annotation functions used by the
  MEDIPS notebook.
- `code/_helpers/map_csf1r_esr_motifs_integrated.py`: ESR1/ESR2 PWM scanning
  of Csf1r-associated DMR windows using a relative score threshold of 0.92.
- `code/ramification_index/ramification_index.py`: Sholl profiles and
  ramification indices calculated from SWC reconstructions.
- `code/MorphOMICs/build_sholl_report_from_swc.py`: reconstruction and
  verification of the selected-cell Sholl table.
- `code/MorphOMICs/plot_ramification_index_layer12_LMM_english.py`: mixed
  models with animal as a random intercept and export of the statistical table.
- `code/imaris_swc_export/`: MATLAB scripts used for Imaris-to-SWC export and
  SWC standardization.

Scripts whose only function was publication-layout rendering are not included.
The numerical source tables and intermediate data required by the retained
analyses are included in the Zenodo archive.

## Environment

- R 4.4.2 and Bioconductor 3.20; exact packages are recorded in `renv.lock`.
- MEDIPS 1.58.0.
- Python 3.13.5; exact packages are recorded in `requirements-python.txt`.

Package names and versions should be verified against CRAN, Bioconductor or
PyPI before installation. The R environment can be restored with `renv` and
the Python environment with a package manager that respects the pinned
requirements file.

## Data layout

After downloading the Zenodo archive, place its `data/` and
`generated_outputs/` directories at the repository root. The expected layout
is:

```text
data/
  MeDIP/
    Aligned_Results_BAM/
    annotations/rn6/
    motifs/
    workspace_and_tables/
  MorphOMICs/
    reports_cam/
    saved_instances_cam/
  SWC_reconstructions/
generated_outputs/
  MeDIP/
  MorphOMICs/
```

The `Aligned_Results_BAM/` folder is **not** included in the Zenodo archive
because of its size. Recreate it from the GEO FASTQ files (HISAT2 v2.1.0, UCSC
rn6, sorted and indexed with samtools) before rendering the MEDIPS notebook, or
request the aligned BAM files from the corresponding author.

## Running the retained analyses

Render the MEDIPS notebook from the repository root:

```text
Rscript -e "rmarkdown::render('code/MeDIP/MEDIPS - Alvaro Bautista-Abad.Rmd')"
```

Rebuild and verify the Sholl table:

```text
python code/MorphOMICs/build_sholl_report_from_swc.py
```

Run the ramification-index mixed models and export their statistical table:

```text
python code/MorphOMICs/plot_ramification_index_layer12_LMM_english.py
```

Run the Csf1r ESR1/ESR2 PWM scan:

```text
python code/_helpers/map_csf1r_esr_motifs_integrated.py
```

The MeDIP-seq DMR workflow and exploratory GO over-representation analysis use
unadjusted p values.
