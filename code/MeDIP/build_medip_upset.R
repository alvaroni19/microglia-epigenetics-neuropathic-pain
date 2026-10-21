# Mirrored DMR gene-set UpSet plot.
suppressMessages({
  library(dplyr); library(tidyr); library(readxl); library(ggplot2)
  library(ComplexUpset); library(patchwork)
})

current_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                         error = function(e) NA_character_)
if (is.na(current_file) || !nzchar(current_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  current_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
}
SCRIPT_DIR <- dirname(normalizePath(current_file, winslash = "/", mustWork = FALSE))
BUNDLE_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/")
DATA <- file.path(BUNDLE_ROOT, "data", "MeDIP", "workspace_and_tables")
OUT <- file.path(BUNDLE_ROOT, "generated_outputs", "MeDIP")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.pt2mm <- function(pt) pt / 2.845276

UPSET_COMPS <- c("F_CSF_vs_F_C", "M_CSF_vs_M_C", "F_C_vs_M_C", "F_CSF_vs_M_CSF")
UPSET_LABELS2 <- c("F_CSF_vs_F_C" = "F-CSF1 vs F-C", "M_CSF_vs_M_C" = "M-CSF1 vs M-C",
                   "F_C_vs_M_C" = "F-C vs M-C", "F_CSF_vs_M_CSF" = "F-CSF1 vs M-CSF1")
# Colour-blind-safe region palette.
REGION_COLORS <- c("TSS" = "#E69F00", "EXON" = "#0072B2", "GENEBODY" = "#009E73")
# Unified Okabe-Ito GO-theme palette (matches the GO rank-distribution panel).
THEME_COLORS2 <- c("Morphology" = "#009E73", "Inflammation / neuroinflammation" = "#D55E00",
                   "Sex hormones / estrogen" = "#CC79A7", "Mechanical stimulus / pain" = "#E69F00",
                   "None" = "grey82")

# Font sizes (pt); nothing below 7 pt.
FS <- list(title = 8, subtitle = 7, panel_title = 7, bar_label = 7,
           axis_title = 7, axis_text = 7, set_axis = 7, matrix_text = 7,
           legend_title = 7, legend_text = 7)

sig <- read_excel(file.path(DATA, "GENES_TOTAL_sig_contig_colap_ws400.xlsx"), sheet = 1)

# candidate symbols
cand <- read_excel(file.path(DATA, "DMR_candidate_genes_ERE_scan_ws400.xlsx"), sheet = "candidate_genes")
cand_lower <- tolower(cand$symbol)

# gene -> GO theme map from the GO rank table
go <- read_excel(file.path(DATA, "GO_relevance_rank_distribution_ws400.xlsx"))
go$biological_theme <- dplyr::recode(
  as.character(go$biological_theme),
  "Steroid / sex hormones" = "Sex hormones / estrogen"
)
gene_theme_map <- go %>%
  dplyr::filter(!is.na(geneID), geneID != "", biological_theme %in% names(THEME_COLORS2)) %>%
  tidyr::separate_rows(geneID, sep = "[/,;|]+") %>%
  dplyr::mutate(gene = tolower(trimws(geneID))) %>%
  dplyr::filter(gene != "", biological_theme != "Other") %>%
  dplyr::distinct(gene, biological_theme)

build_upset_df <- function(direction) {
  base <- sig %>%
    dplyr::filter(diffexpressed == direction, comparison %in% UPSET_COMPS,
                  !is.na(external_transcript_gene_name)) %>%
    dplyr::mutate(gene = tolower(as.character(external_transcript_gene_name)),
                  region = dplyr::recode(as.character(region), "GENE" = "GENEBODY")) %>%
    dplyr::distinct(gene, comparison, region)
  if (nrow(base) == 0) return(tibble::tibble())
  rp <- c("TSS" = 1L, "EXON" = 2L, "GENEBODY" = 3L)
  primary_region <- base %>%
    dplyr::mutate(rpv = dplyr::coalesce(rp[region], 4L)) %>%
    dplyr::group_by(gene) %>% dplyr::slice_min(rpv, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
    dplyr::transmute(gene, primary_region = factor(region, levels = names(REGION_COLORS)))
  memb <- base %>% dplyr::distinct(gene, comparison) %>% dplyr::mutate(present = TRUE) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = present, values_fill = FALSE)
  for (cp in UPSET_COMPS) if (!cp %in% names(memb)) memb[[cp]] <- FALSE
  for (old in names(UPSET_LABELS2)) if (old %in% names(memb)) names(memb)[names(memb) == old] <- UPSET_LABELS2[old]
  prim_theme <- gene_theme_map %>% dplyr::group_by(gene) %>% dplyr::slice(1) %>% dplyr::ungroup() %>%
    dplyr::transmute(gene, go_primary_theme = factor(biological_theme, levels = names(THEME_COLORS2)))
  memb %>%
    dplyr::left_join(primary_region, by = "gene") %>%
    dplyr::left_join(prim_theme, by = "gene") %>%
    dplyr::mutate(
      primary_region = dplyr::coalesce(primary_region, factor("GENEBODY", levels = names(REGION_COLORS))),
      go_primary_theme = dplyr::coalesce(go_primary_theme, factor("None", levels = names(THEME_COLORS2))),
      is_candidate = gene %in% cand_lower
    )
}

make_upset <- function(df, sets, bar_fill, direction_title) {
  sets_ok <- intersect(sets, names(df))
  p1 <- ComplexUpset::intersection_size(
      mapping = aes(fill = primary_region), text = list(size = .pt2mm(FS$bar_label)),
      bar_number_threshold = 1) +
    scale_fill_manual(values = REGION_COLORS, name = "Region", drop = FALSE) +
    labs(title = direction_title, y = "n genes") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = FS$panel_title, color = bar_fill),
          axis.title = element_text(size = FS$axis_title), axis.text = element_text(size = FS$axis_text),
          legend.title = element_text(size = FS$legend_title), legend.text = element_text(size = FS$legend_text),
          legend.position = "right", legend.key.size = unit(2.4, "mm"))
  p2 <- ggplot(mapping = aes(x = intersection)) +
    geom_bar(data = function(d) dplyr::filter(d, is_candidate), fill = "#8B0000", alpha = 0.85, width = 0.55) +
    ylab("Candidate") +
    scale_y_continuous(breaks = function(x) unique(floor(pretty(x))), expand = expansion(mult = c(0, 0.2))) +
    theme(axis.title = element_text(size = FS$axis_title), axis.text = element_text(size = FS$axis_text))
  p3 <- ggplot(mapping = aes(x = intersection)) +
    geom_bar(data = function(d) dplyr::mutate(d, go_primary_theme = factor(go_primary_theme, levels = names(THEME_COLORS2))),
             mapping = aes(fill = go_primary_theme), position = "fill", width = 0.55) +
    scale_fill_manual(values = THEME_COLORS2, name = "GO theme", drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
    ylab("GO theme") +
    theme(axis.title = element_text(size = FS$axis_title), axis.text = element_text(size = FS$axis_text),
          legend.title = element_text(size = FS$legend_title), legend.text = element_text(size = FS$legend_text),
          legend.position = "right", legend.key.size = unit(2.4, "mm"))
  ComplexUpset::upset(
    data = df, intersect = sets_ok, name = NULL, min_size = 1, width_ratio = 0.20, sort_sets = FALSE,
    matrix = ComplexUpset::intersection_matrix(
      geom = geom_point(size = 0.70),
      segment = geom_segment(linewidth = 0.35)
    ),
    set_sizes = FALSE,
    base_annotations = list("n genes" = p1),
    themes = ComplexUpset::upset_modify_themes(list(
      "intersections_matrix" = theme(text = element_text(size = FS$matrix_text)),
      "overall_sizes" = theme(text = element_text(size = FS$axis_text))))
  )
}

hyper <- build_upset_df("Hyper"); hypo <- build_upset_df("Hypo")
sets_display <- intersect(unname(UPSET_LABELS2), names(hyper))
p_hi <- make_upset(hyper, sets_display, "#C0392B", "Hypermethylated")
p_lo <- make_upset(hypo, sets_display, "#2471A3", "Hypomethylated")
p_mirror <- patchwork::wrap_plots(p_hi, p_lo, ncol = 1, heights = c(1, 1)) +
  patchwork::plot_annotation(
    title = "DMR gene-set membership — Hyper (top) / Hypo (bottom)",
    theme = theme(plot.title = element_text(face = "bold", size = FS$title)))

W <- 180; H <- 150
base <- file.path(OUT, "medip_dmr_gene_set_upset")
ragg::agg_tiff(paste0(base, ".tiff"), width = W, height = H, units = "mm",
               res = 600, background = "white")
print(p_mirror); dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = W / 25.4, height = H / 25.4, bg = "white")
print(p_mirror); dev.off()
svglite::svglite(paste0(base, ".svg"), width = W / 25.4, height = H / 25.4, bg = "white")
print(p_mirror); dev.off()
cat("[OK] medip_dmr_gene_set_upset (tiff + pdf + svg)\n")
