assert_rn6_bam_headers <- function(bam_files, genome) {
  expected <- GenomeInfoDb::seqlengths(genome)
  required_chromosomes <- paste0("chr", c(1:20, "X", "Y", "M"))
  required_chromosomes <- intersect(required_chromosomes, names(expected))

  for (bam_file in bam_files) {
    observed <- Rsamtools::scanBamHeader(bam_file)[[1]]$targets
    comparable <- intersect(required_chromosomes, names(observed))
    if (length(comparable) != length(required_chromosomes)) {
      stop("BAM header is missing rn6 reference sequences: ", basename(bam_file))
    }
    mismatched <- comparable[observed[comparable] != expected[comparable]]
    if (length(mismatched) > 0) {
      stop(
        "BAM reference lengths do not match rn6 in ", basename(bam_file),
        ": ", paste(mismatched, collapse = ", ")
      )
    }
  }
  invisible(TRUE)
}

normalize_rn6_seqnames <- function(values) {
  values <- as.character(values)
  values[values == "MT"] <- "M"
  ifelse(grepl("^chr", values), values, paste0("chr", values))
}

build_rn6_medips_annotations <- function(gtf_path, orgdb) {
  if (!file.exists(gtf_path)) stop("Missing bundled rn6 GTF: ", gtf_path)
  gtf <- rtracklayer::import(gtf_path)
  metadata <- as.data.frame(S4Vectors::mcols(gtf), stringsAsFactors = FALSE)
  feature <- as.character(gtf$type)
  chromosome <- normalize_rn6_seqnames(GenomicRanges::seqnames(gtf))
  start_position <- BiocGenerics::start(gtf)
  end_position <- BiocGenerics::end(gtf)
  strand_value <- as.character(BiocGenerics::strand(gtf))

  gene_symbols <- unique(stats::na.omit(as.character(metadata$gene_name)))
  descriptions <- suppressMessages(AnnotationDbi::mapIds(
    orgdb,
    keys = gene_symbols,
    keytype = "SYMBOL",
    column = "GENENAME",
    multiVals = "first"
  ))

  make_information <- function(feature_name, id_column) {
    keep <- feature == feature_name & !is.na(metadata[[id_column]])
    symbol <- as.character(metadata$gene_name[keep])
    description <- unname(descriptions[symbol])
    description[is.na(description) | description == ""] <- symbol[
      is.na(description) | description == ""
    ]
    result <- data.frame(
      ID = as.character(metadata[[id_column]][keep]),
      external_gene_name = symbol,
      gene_biotype = as.character(metadata$gene_biotype[keep]),
      description = description,
      stringsAsFactors = FALSE
    )
    result[!duplicated(result$ID), , drop = FALSE]
  }

  make_annotation <- function(feature_name, id_column) {
    keep <- feature == feature_name & !is.na(metadata[[id_column]])
    result <- data.frame(
      id = as.character(metadata[[id_column]][keep]),
      chr = chromosome[keep],
      start = as.integer(start_position[keep]),
      end = as.integer(end_position[keep]),
      stringsAsFactors = FALSE
    )
    result[!duplicated(result), , drop = FALSE]
  }

  gene_annotation <- make_annotation("gene", "gene_id")
  exon_annotation <- make_annotation("exon", "exon_id")

  transcript_rows <- feature == "transcript" & !is.na(metadata$transcript_id)
  tss <- ifelse(
    strand_value[transcript_rows] == "-",
    end_position[transcript_rows],
    start_position[transcript_rows]
  )
  tss_annotation <- data.frame(
    id = as.character(metadata$transcript_id[transcript_rows]),
    chr = chromosome[transcript_rows],
    start = as.integer(ifelse(
      strand_value[transcript_rows] == "-", tss - 500, tss - 1000
    )),
    end = as.integer(ifelse(
      strand_value[transcript_rows] == "-", tss + 1000, tss + 500
    )),
    stringsAsFactors = FALSE
  )
  tss_annotation$start <- pmax(1L, tss_annotation$start)
  tss_annotation <- tss_annotation[!duplicated(tss_annotation), , drop = FALSE]

  csf1r <- gene_annotation[gene_annotation$id == "ENSRNOG00000018414", ]
  if (nrow(csf1r) != 1 || csf1r$chr != "chr18" ||
      csf1r$start != 56414488L || csf1r$end != 56458300L) {
    stop("The bundled GTF does not contain the expected rn6 Csf1r coordinates.")
  }

  list(
    annotations = list(
      TSS = list(TSS = tss_annotation),
      EXON = list(EXON = exon_annotation),
      GENE = list(Gene = gene_annotation)
    ),
    information = list(
      TSS = make_information("transcript", "transcript_id"),
      EXON = make_information("exon", "exon_id"),
      GENE = make_information("gene", "gene_id")
    )
  )
}
