# microglia-epigenetics-neuropathic-pain: analysis code and data

Analysis code for a study of microglial DNA methylation, sex, estrous phase and trigeminal
neuropathic pain: MeDIP-seq of rat primary microglia treated with CSF1, RT-qPCR of mouse
brainstem after chronic constriction injury of the infraorbital nerve (CCI-IoN) and of primary
microglia treated with CSF1 and 17β-estradiol, von Frey behaviour, and microglial morphology
(Sholl analysis and MorphOMICs). Please cite this archive (`CITATION.cff`).

## Data availability

| Data | Location |
|---|---|
| MeDIP-seq FASTQ and WIG files | NCBI GEO, Series [GSE338843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338843) |
| Processed MeDIP-seq tables, rn6 annotation, ESR1/ESR2 motif matrices (JASPAR MA0112.3, MA0258.2) | Zenodo archive, `data/MeDIP/` |
| Microglial SWC reconstructions, Sholl reports and MorphOMICs objects | Zenodo archive, `data/SWC_reconstructions/`, `data/MorphOMICs/` |
| Behavioural (von Frey) and RT-qPCR data, and all results derived from them | Not distributed; available from the corresponding author on reasonable request |
| Aligned rn6 BAM files (about 32 GB) | Not distributed; regenerate from GEO (see below) or request them |

The Zenodo archive (concept DOI [10.5281/zenodo.21396205](https://doi.org/10.5281/zenodo.21396205),
which resolves to the latest version) contains the code, the data above and the regenerated
MeDIP and morphology outputs. This GitHub repository contains the code only; copy the `data/`
folder of the Zenodo archive to the repository root to run it.

The behaviour, RT-qPCR and multivariate code is included. When the corresponding inputs are
placed in `data/` (file names in "Input files"), the same commands regenerate those analyses;
without them, these steps are reported as `[SKIP]` and the remaining analyses run normally.

## Contents

| Folder | Analysis |
|---|---|
| `code/MeDIP/` | MEDIPS/edgeR analysis in 400-bp windows (`medips_differential_methylation.Rmd`), rn6 annotation, DMR directionality, GO rank distributions, candidate-gene heatmap and window map, robust PCA, global methylation heatmap, DMR gene-set UpSet, Csf1r locus map |
| `code/MorphOMICs/`, `code/ramification_index/` | Sholl ramification index with linear mixed models, MorphOMICs persistence images, persistence differences and UMAP, ramification index versus mechanical response |
| `code/imaris_swc_export/` | MATLAB scripts for Imaris-to-SWC export and SWC standardization |
| `code/PCR_260514/` | RT-qPCR in vivo and in vitro: animal-level reconstruction, two-way ANOVA with Bonferroni cell-mean contrasts, Student's t tests, violins, bubble plots, in vitro estrogenic score and PCA |
| `code/PCR_MFA/` | Multiple factor analysis with PERMANOVA and dispersion tests, marker-pain Spearman correlations, estrogenic score versus mechanical response |
| `code/Behaviour/` | Mechanical response (von Frey 0.16 g): change-score linear model and sensitivity analyses |
| `code/Figures_v2/` | Panels drawn from the result tables (`generated_outputs/Figures_v2/panels/`) and the same panels at final printed size (`generated_outputs/Figures_v2/blocks/`; sizes in `figure_layout.py`) |
| `code/_helpers/` | Shared loaders, statistics lookups and plotting style |

## Running the analyses

From the repository root:

```text
python code/reproduce_all_figures.py                       # every domain whose inputs are present
python code/reproduce_all_figures.py --domain morphology   # one domain (may be repeated)
python code/reproduce_all_figures.py --skip-r              # Python steps only
```

Domains: `medip`, `morphology`, `qpcr`, `behaviour`, `multivariate`, `figures_v2`. Outputs are
written to `generated_outputs/<domain folder>/`. `Rscript` must be on `PATH` (or set `RSCRIPT`).

The MEDIPS notebook needs the BAM files:

```text
Rscript -e "rmarkdown::render('code/MeDIP/medips_differential_methylation.Rmd')"
```

To recreate `data/MeDIP/Aligned_Results_BAM/`, download the FASTQ files from GEO GSE338843,
align them to UCSC rn6 with HISAT2 v2.1.0, and sort and index the alignments with samtools,
keeping the file names used in the notebook.

## Panels

Block identifiers (file stems in `generated_outputs/Figures_v2/blocks/`; the panel versions in
`panels/` use the same stems). Blocks marked † need the non-distributed behaviour or RT-qPCR inputs.

| Block | Content |
|---|---|
| `medip_dmr_directionality` | Proportion of hypo- and hypermethylated DMRs: CSF1 effect within sex and sex effect |
| `medip_go_rank_distribution` | GO over-representation ranks of DMR-associated genes, highlighted terms by biological theme |
| `medip_candidate_gene_heatmap` | Mean edgeR log2 fold change of candidate genes across contrasts |
| `medip_candidate_window_map` | Significant 400-bp windows of candidate genes by gene region and contrast |
| `medip_csf1r_locus_dmr_ere_map` | rn6 Csf1r locus: CpG density, DMR windows and ESR1/ESR2 PWM matches |
| `medip_robpca` | Robust PCA of the significant-gene methylation matrix |
| `medip_global_methylation_heatmap` | Methylation of ROBPCA-eligible significant genes |
| `medip_dmr_gene_set_upset` | Mirrored UpSet plot of hyper- and hypomethylated gene sets |
| `morphology_ramification_index` | Sholl ramification index, laminae I-II, ipsilateral, control vs CCI-IoN per group |
| `morphology_umap` | UMAP of MorphOMICs persistence images (laminae I-II, ipsilateral) |
| `morphology_persistence_difference` | CCI-IoN minus control persistence images with representative SWC traces |
| `morphology_ramification_pain_partial_spearman` † | Ramification index versus mechanical response (partial Spearman) |
| `invitro_qpcr_methylation_regulators_csf1_bubble` † | In vitro CSF1 vs control, DNA-methylation regulators |
| `invitro_qpcr_methylation_regulators_bubble` † | In vitro CSF1 and CSF1+E2 contrasts, DNA-methylation regulators |
| `invitro_qpcr_csf1r_violin` † | In vitro Csf1r expression, six groups |
| `invitro_qpcr_methylation_regulators_csf1_violins` † | In vitro DNA-methylation regulators, control vs CSF1 by sex |
| `invitro_qpcr_methylation_regulators_violins` † | In vitro DNA-methylation regulators, six groups |
| `invitro_qpcr_estrogen_genes_violins` † | In vitro Esr1, Esr2, Gper1 and Cyp19a1, six groups |
| `invitro_estrogen_genes_bubble` † | In vitro CSF1 and CSF1+E2 contrasts, estrogen-related genes |
| `invitro_estrogenic_score`, `invitro_estrogen_genes_pca` † | In vitro directed estrogenic score and PCA of estrogen-related genes |
| `behaviour_change_score_violins` † | Change in mechanical response after CCI-IoN by sex and estrous phase |
| `invivo_qpcr_pain_microglia_markers` † | In vivo pain-related and microglial markers |
| `invivo_qpcr_methylation_regulators`, `invivo_qpcr_gadd45ab_violins` † | In vivo DNA-methylation regulators |
| `invivo_qpcr_estrogen_genes_violins`, `invivo_estrogen_genes_bubble` † | In vivo estrogen-related genes |
| `invivo_marker_pain_spearman` † | Spearman correlation of markers with mechanical response |
| `estrogenic_score_vs_pain` † | In vivo estrogenic score versus mechanical response |
| `invivo_mfa_phase_condition` † | Multiple factor analysis by estrous phase and condition |

## Analysis notes

- **MeDIP-seq.** The DMR workflow and the exploratory GO over-representation analysis use
  unadjusted P values. Result tables store contrasts as CSF1 versus control and female versus
  male; panels show CSF1 contrasts CSF1 first (F-CSF1 vs F-C, M-CSF1 vs M-C) and sex contrasts
  male first (M-C vs F-C, M-CSF1 vs F-CSF1; log2 fold changes sign-inverted). "Hyper" denotes
  higher methylation in the first-named group.
- **Morphology.** Analyses use ipsilateral cells of laminae I-II. The ramification index is
  compared between control and CCI-IoN within each sex/estrous group with a linear mixed model
  (condition as fixed effect, animal as random intercept). The UMAP shows one global embedding,
  zoomed on the main cloud, with an inset showing all bootstrap bags.
- **Behaviour.** The outcome is the per-animal change score (test minus mean of three
  pre-surgery baseline days, 0.16 g filament), modelled with `lm(Delta ~ Condition * Group)`, HC3
  standard errors and type-III tests; injury effects are Holm-adjusted across groups and phase
  comparisons among CCI-IoN females within the female family. Baseline-adjusted ANCOVA and a
  Freedman-Lane permutation test are sensitivity analyses.
- **RT-qPCR.** In vivo values are reconstructed from animal-level Delta-Ct values with male
  controls as calibrator; two Cyp19a1 measurements are excluded for technical reasons
  (`TECHNICAL_EXCLUSIONS` in `code/_helpers/qpcr_raw.py`). Statistics are computed on Delta-Ct;
  fold changes are displayed relative to male controls (males) or diestrus controls (females).
- **Multivariate analysis.** Four active blocks (mechanical response, microglial markers,
  estrogen-related genes, DNA-methylation genes) after missMDA imputation; PERMANOVA and
  dispersion tests on the first four dimensions (9,999 permutations). Marker-pain Spearman
  correlations use the ten microglial and DNA-methylation markers with Benjamini-Hochberg
  correction.

## Input files

```text
data/
  MeDIP/                     distributed: workspace_and_tables/, annotations/rn6/, motifs/,
                             MEDIP_PROVIDER_QC_METRICS.tsv; Aligned_Results_BAM/ not distributed
  MorphOMICs/                distributed: reports_cam/, saved_instances_cam/
                             not distributed: morphology_behaviour_linkage.xlsx
  SWC_reconstructions/       distributed: <condition>/<group>/<sex>/<side>/<layer>/<animal>/*.swc
  Allodynia_0_16/            not distributed: allodynia_change_scores.xlsx, allodynia_von_frey_scores.xlsx
  PCR_260514/                not distributed: invivo_qpcr_delta_ct.xlsx, invitro_qpcr_relative_expression.xlsx
```

## Environment

- Python packages: `requirements-python.txt`; R packages for the MEDIPS notebook: `renv.lock`.
- Versions used to validate this release: `ENVIRONMENT_v2.md`.

## License

Code: MIT (see `LICENSE`).
