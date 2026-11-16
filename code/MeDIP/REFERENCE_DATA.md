# MeDIP reference data included in this archive

- `data/MeDIP/annotations/rn6/Rattus_norvegicus.Rnor_6.0.104.gtf.gz`:
  Ensembl Rnor_6.0 release 104 annotation used for every gene assignment and
  gene model. The workflow checks the BAM sequence lengths against rn6 before
  annotation.
- `data/MeDIP/motifs/MA0112.3.meme` and `MA0258.2.meme`:
  ESR1 and ESR2 position-weight matrices used for the Csf1r-window motif scan.
- `data/MeDIP/workspace_and_tables/csf1r_dmr_windows_rn6.fasta`:
  rn6 sequences exported by the MEDIPS notebook for DMR windows that overlap
  the rn6 Csf1r locus.
- `data/MeDIP/workspace_and_tables/csf1r_context_rn6.fasta`:
  rn6 Csf1r locus and flanking sequence used for the CpG-density track.
- `data/MeDIP/workspace_and_tables/significant_gene_rms_matrix_ws400.csv` and
  `medip_sample_metadata.csv`: plotting data exported by the MEDIPS notebook
  for the robust PCA and global methylation heatmap.

The retained workflow does not query Ensembl or BioMart. All genomic
annotation is derived locally from the bundled rn6 GTF.
