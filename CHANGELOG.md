# Changelog

## 2.0.0 (2026-10-20)

### Added
- Behaviour analysis code (`code/Behaviour/`): change-score linear model with HC3 standard
  errors, Holm-adjusted planned contrasts, baseline-adjusted ANCOVA and Freedman-Lane permutation
  test.
- RT-qPCR analysis code (`code/PCR_260514/`), including the in vitro directed estrogenic score and
  PCA (`build_invitro_estrogenic_score.R`), and the multiple factor analysis (`code/PCR_MFA/`).
- `code/Figures_v2/`: scripts that draw every data panel from the result tables, and the same
  panels at final printed size (`build_figure_blocks.py`; sizes in `figure_layout.py`) with one
  type scale (Arial; tick labels 6 pt, titles 7 pt bold, asterisks 8 pt).
- `code/ramification_index/` and `code/imaris_swc_export/`.
- `code/reproduce_all_figures.py` runs every domain whose inputs are present and reports the
  others as skipped.

### Changed
- Behavioural and RT-qPCR data, and all results derived from them, are no longer distributed;
  they are available from the corresponding author on request. The code is distributed.
- Code, comments, column names, group labels and file names are in English; outputs are named
  after the analysis they contain. MorphOMICs tables and objects use the columns `Condition`,
  `Phase`, `Sex`, `Side` and `Layer` with English group labels; the SWC tree uses English folder
  names. Numerical content is unchanged.
- In vivo RT-qPCR: two Cyp19a1 measurements excluded for technical reasons, as in the original
  GraphPad Prism analysis.
- Ramification-index mixed models restricted to ipsilateral cells, consistent with the UMAP and
  persistence analyses.
- MeDIP panels write CSF1 contrasts CSF1 first (F-CSF1 vs F-C, M-CSF1 vs M-C) and sex contrasts
  male first (M-C vs F-C, M-CSF1 vs F-CSF1); for display, sex contrasts are sign-inverted and
  Hyper/Hypo exchanged. "Hyper" means higher methylation in the first-named group.
- Csf1r locus map shows all ESR1/ESR2 PWM matches (relative score >= 0.92) along the rn6 locus
  and highlights the matches inside DMR windows.
- Violins use untrimmed Gaussian densities (R bw.nrd0 bandwidth, tails to min/max +/- 3 bw).
- Marker-pain Spearman correlations use the ten microglial and DNA-methylation markers with
  Benjamini-Hochberg correction; estrogen-related genes are analysed through the estrogenic score.
- Multiple factor analysis keeps four active blocks; phase ellipses are 68% normal data ellipses.
- Morphology UMAP: one global embedding of ipsilateral lamina I-II bags, zoomed on the main
  cloud (robust distance <= 6 MAD units) with a full-view inset.
- Candidate-gene heatmap: genes without a DMR window in either within-sex contrast are not drawn.
- MorphOMICs panels label the injured condition as CCI-IoN.

### Removed
- Paired Wilcoxon baseline-versus-test figures as the primary behavioural analysis.
- Exploratory notebooks superseded by the scripts above.

## 1.0.0

- Initial archive.
