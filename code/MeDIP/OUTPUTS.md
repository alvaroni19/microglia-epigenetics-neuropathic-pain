# MeDIP release outputs

The canonical MEDIPS notebook writes the processed tables consumed by the
final figure scripts directly to `data/MeDIP/workspace_and_tables/`:

- `GENES_TOTAL_sig_contig_colap_ws400.xlsx` — DMR directionality and DMR
  gene-set UpSet panels.
- `GO_relevance_rank_distribution_ws400.xlsx` — GO rank-distribution panel and
  UpSet GO themes.
- `methylation_map_candidate_genes_hits_ws400.xlsx` — candidate-gene heatmap.
- `heatmap_candidate_genes_logFC_ws400_orders.xlsx` — deterministic row and
  column clustering order used by the candidate-gene heatmap.
- `DMR_candidate_genes_ERE_scan_ws400.xlsx` — candidate-gene support table
  used by the DMR gene-set UpSet panel.
- `significant_gene_rms_matrix_ws400.csv` and `medip_sample_metadata.csv` —
  inputs for the robust PCA and global methylation heatmap.
- `csf1r_dmr_windows_rn6.fasta` — rn6 Csf1r DMR-window sequences used by the
  bundled ESR1/ESR2 PWM scan and the Csf1r locus map.
- `csf1r_context_rn6.fasta` — rn6 Csf1r locus and flanking sequence used for
  the Csf1r CpG-density track.

No frozen BioMart cache or serialized R workspace is required. The Csf1r locus map is
built from the rn6 GTF, the regenerated candidate-gene table, the regenerated
FASTA file and the bundled motif matrices.
