# MeDIP overview panels: robust PCA, global methylation heatmap and DMR gene-set UpSet.
#   a  ROBPCA of the significant-gene RMS matrix (same computation as
#      code/MeDIP/build_medip_robpca.R)
#   b  heatmap of the 80 most variable ROBPCA-eligible genes (same selection
#      as code/MeDIP/build_medip_global_methylation.R)
#   c  mirrored UpSet of DMR gene sets (same membership rules as
#      code/MeDIP/build_medip_upset.R)
# MeDIP contrasts are written CSF1 first (F-CSF1 vs F-C, M-CSF1 vs M-C) and
# male first (M-C vs F-C, M-CSF1 vs F-CSF1). The release tables store the CSF1
# contrasts in the same order and the sex contrasts female first, so
# 'Hypermethylated' here (higher in the first group) is the release direction
# 'Hyper' for the CSF1 contrasts and 'Hypo' for the sex contrasts.
# Usage: Rscript build_medip_overview_panels.R [BUNDLE_ROOT]
suppressMessages({
  library(dplyr); library(tidyr); library(readxl); library(ggplot2); library(tibble)
  library(ComplexUpset); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  ROOT <- normalizePath(args[[1]], winslash = "/")
} else {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1]])), "..", ".."), winslash = "/")
}
DATA <- file.path(ROOT, "data", "MeDIP", "workspace_and_tables")
OUT <- file.path(ROOT, "generated_outputs", "Figures_v2", "panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.pt2mm <- function(pt) pt / 2.845276
# Register Arial (regular/bold/italic) when the TTF files are available, so
# ragg/svglite render italic gene symbols; otherwise the system sans is used.
font_dirs <- c(Sys.getenv("ARIAL_FONT_DIR"), file.path(Sys.getenv("CONDA_PREFIX"), "fonts"),
               "/usr/share/fonts/truetype/msttcorefonts", "C:/Windows/Fonts")
for (d in font_dirs[nzchar(font_dirs)]) {
  f <- file.path(d, c("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"))
  if (all(file.exists(f))) {
    systemfonts::register_font("Arial", plain = f[1], bold = f[2], italic = f[3], bolditalic = f[4])
    break
  }
}

save_all <- function(p, stem, w, h) {
  ragg::agg_png(file.path(OUT, paste0(stem, ".png")), width = w, height = h, units = "mm",
                res = 600, background = "white"); print(p); invisible(dev.off())
  grDevices::cairo_pdf(file.path(OUT, paste0(stem, ".pdf")), width = w / 25.4, height = h / 25.4,
                       bg = "white", family = "Arial"); print(p); invisible(dev.off())
  svglite::svglite(file.path(OUT, paste0(stem, ".svg")), width = w / 25.4, height = h / 25.4,
                   bg = "white"); print(p); invisible(dev.off())
}

theme_ms <- function(base = 7) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(text = element_text(colour = "black"),
          axis.title = element_text(size = 7, face = "bold"),
          axis.text = element_text(size = 6.5, colour = "black"),
          axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5),
          legend.key.size = unit(2.6, "mm"), plot.margin = margin(2, 2, 1, 1))
}

# ----------------------------------------------------------------------------- matrix
mt <- utils::read.csv(file.path(DATA, "significant_gene_rms_matrix_ws400.csv"), check.names = FALSE)
mat <- as.matrix(mt[, -1, drop = FALSE]); rownames(mat) <- mt$gene
keep <- apply(mat, 1, function(v) all(is.finite(v))) & apply(mat, 1, function(v) stats::mad(v) > 0) &
  !grepl("^(AA|AC)", rownames(mat), ignore.case = TRUE)
mat_use <- mat[keep, , drop = FALSE]
meta <- utils::read.csv(file.path(DATA, "medip_sample_metadata.csv"), stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------- a: ROBPCA
set.seed(1)
rob <- rrcov::PcaHubert(t(mat_use), k = min(5, ncol(mat_use) - 1, nrow(mat_use) - 1), scale = TRUE)
sc <- as.data.frame(rob@scores[, 1:2, drop = FALSE]); names(sc) <- c("RPC1", "RPC2")
sc <- sc %>% rownames_to_column("sample_label") %>% left_join(meta, by = "sample_label")
grp_lab <- c("MC" = "M-C", "MCSF1" = "M-CSF1", "FC" = "F-C", "FCSF1" = "F-CSF1")
sc$Group <- factor(grp_lab[as.character(sc$group)], levels = unname(grp_lab))
sc$sex_col <- ifelse(grepl("^M", sc$group), "#F2A16A", "#6666FF")
write.csv(sc, file.path(OUT, "medip_robpca_scores.csv"), row.names = FALSE)
grp_col <- c("M-C" = "#F2A16A", "M-CSF1" = "#F2A16A", "F-C" = "#6666FF", "F-CSF1" = "#6666FF")
grp_shp <- c("M-C" = 16, "M-CSF1" = 17, "F-C" = 16, "F-CSF1" = 17)
pa <- ggplot(sc, aes(RPC1, RPC2, colour = Group, shape = Group)) +
  geom_point(size = 1.8) +
  ggrepel::geom_text_repel(aes(label = sub("^([MF])C", "\\1C", sample_label)), size = .pt2mm(6.5),
                           family = "Arial", show.legend = FALSE, seed = 1, box.padding = 0.2,
                           min.segment.length = Inf, max.overlaps = Inf) +
  scale_colour_manual(values = grp_col, name = "Group") +
  scale_shape_manual(values = grp_shp, name = "Group") +
  labs(x = "Robust PC1", y = "Robust PC2") +
  theme_ms() + theme(legend.position = c(0.86, 0.55), legend.background = element_blank(),
                     legend.key.spacing.y = unit(0, "mm"))
save_all(pa, "medip_robpca", 82, 66)

# ----------------------------------------------------------------------------- b: heatmap
top <- mat_use
if (nrow(top) > 80) top <- top[order(apply(top, 1, stats::var), decreasing = TRUE)[1:80], , drop = FALSE]
z <- t(scale(t(top))); z[!is.finite(z)] <- 0
ord <- function(x) {
  cc <- suppressWarnings(stats::cor(t(x))); cc[!is.finite(cc)] <- 0; diag(cc) <- 1
  stats::hclust(stats::as.dist(1 - cc), method = "average")$order
}
z <- z[ord(z), ord(t(z)), drop = FALSE]
relab <- function(s) { s <- sub("^MCSF1", "M-CSF1", s); s <- sub("^FCSF1", "F-CSF1", s)
  s <- sub("^MC", "M-C", s); sub("^FC", "F-C", s) }
gene_levels <- rev(rownames(z)); sample_levels <- colnames(z)
dfz <- as.data.frame(z) %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "zscore") %>%
  mutate(gene = factor(gene, levels = gene_levels), sample = factor(sample, levels = sample_levels))
pb <- ggplot(dfz, aes(sample, gene, fill = zscore)) + geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#EF7660", midpoint = 0, name = "Row z-score") +
  scale_x_discrete(expand = c(0, 0), labels = relab) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_ms(6) +
  theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 5.2, face = "italic"), axis.line = element_blank(),
        legend.title = element_text(angle = 270, size = 7, face = "bold"),
        legend.key.height = unit(14, "mm"), legend.key.width = unit(2.4, "mm")) +
  guides(fill = guide_colourbar(title.position = "right", title.hjust = 0.5))
save_all(pb, "medip_global_methylation_heatmap", 74, 150)

# ----------------------------------------------------------------------------- c: UpSet
sig <- read_excel(file.path(DATA, "GENES_TOTAL_sig_contig_colap_ws400.xlsx"), sheet = 1)
LAB <- c("F_CSF_vs_F_C" = "F-CSF1 vs F-C", "M_CSF_vs_M_C" = "M-CSF1 vs M-C",
         "F_C_vs_M_C" = "M-C vs F-C", "F_CSF_vs_M_CSF" = "M-CSF1 vs F-CSF1")
SETS <- c("M-CSF1 vs F-CSF1", "M-C vs F-C", "M-CSF1 vs M-C", "F-CSF1 vs F-C")
# release direction per contrast for a displayed direction: CSF1 contrasts are displayed as
# stored (CSF1 first); sex contrasts are stored female first and displayed male first
FLIP <- c("F_CSF_vs_F_C" = FALSE, "M_CSF_vs_M_C" = FALSE, "F_C_vs_M_C" = TRUE, "F_CSF_vs_M_CSF" = TRUE)
REGION_COLORS <- c("TSS" = "#E69F00", "EXON" = "#0072B2", "GENEBODY" = "#009E73")
build <- function(display_direction) {
  other <- c(Hyper = "Hypo", Hypo = "Hyper")[[display_direction]]
  base <- sig %>% filter(comparison %in% names(LAB), !is.na(external_transcript_gene_name)) %>%
    filter(diffexpressed == ifelse(FLIP[as.character(comparison)], other, display_direction)) %>%
    mutate(gene = tolower(as.character(external_transcript_gene_name)),
           region = recode(as.character(region), "GENE" = "GENEBODY")) %>%
    distinct(gene, comparison, region)
  rp <- c("TSS" = 1L, "EXON" = 2L, "GENEBODY" = 3L)
  prim <- base %>% mutate(r = coalesce(rp[region], 4L)) %>% group_by(gene) %>%
    slice_min(r, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(gene, Region = factor(region, levels = names(REGION_COLORS)))
  memb <- base %>% distinct(gene, comparison) %>% mutate(v = TRUE) %>%
    pivot_wider(names_from = comparison, values_from = v, values_fill = FALSE)
  for (k in names(LAB)) { if (!k %in% names(memb)) memb[[k]] <- FALSE; names(memb)[names(memb) == k] <- LAB[[k]] }
  memb %>% left_join(prim, by = "gene") %>% mutate(Region = coalesce(Region, factor("GENEBODY", levels = names(REGION_COLORS))))
}
up <- function(df, title, col, show_legend) {
  upset(df, intersect = SETS, name = NULL, min_size = 1, width_ratio = 0.001, sort_sets = FALSE,
        set_sizes = FALSE, n_intersections = 16,
        matrix = intersection_matrix(geom = geom_point(size = 0.8), segment = geom_segment(linewidth = 0.35)),
        base_annotations = list("n genes" = intersection_size(
          mapping = aes(fill = Region), text = list(size = .pt2mm(6), vjust = -0.1), bar_number_threshold = 1) +
            scale_fill_manual(values = REGION_COLORS, drop = FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
            labs(title = title, y = "n genes") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5, colour = col),
                  axis.title.y = element_text(size = 6.5), axis.text = element_text(size = 6),
                  legend.position = if (show_legend) c(0.86, 0.78) else "none",
                  legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
                  legend.key.size = unit(2.2, "mm"), legend.background = element_blank())),
        themes = upset_modify_themes(list(intersections_matrix = theme(text = element_text(size = 6.5)))))
}
hyper_df <- build("Hyper")   # higher methylation in the first group of each label
hypo_df <- build("Hypo")
counts <- bind_rows(
  hyper_df %>% summarise(across(all_of(SETS), sum)) %>% mutate(direction = "Hypermethylated"),
  hypo_df %>% summarise(across(all_of(SETS), sum)) %>% mutate(direction = "Hypomethylated"))
write.csv(counts, file.path(OUT, "medip_dmr_gene_set_upset_set_sizes.csv"), row.names = FALSE)
pc <- wrap_plots(up(hyper_df, "Hypermethylated", "#C0392B", TRUE),
                 up(hypo_df, "Hypomethylated", "#2471A3", FALSE), ncol = 1)
save_all(pc, "medip_dmr_gene_set_upset", 100, 98)

# ----------------------------------------------------------------------------- slot-size blocks
# Same panels redrawn at their final printed size
# (sizes in mm; see code/Figures_v2/figure_layout.py, SUPP[1]) with the shared type sizes:
# tick labels 6 pt, axis/legend titles 7 pt bold, legend text 6 pt.
BLK <- file.path(ROOT, "generated_outputs", "Figures_v2", "blocks")
dir.create(BLK, recursive = TRUE, showWarnings = FALSE)
save_block <- function(p, id, w, h) { OUT_old <- OUT; OUT <<- BLK; on.exit(OUT <<- OUT_old); save_all(p, id, w, h) }
lw <- 0.6 / .pt
blk_theme <- theme(axis.title = element_text(size = 7, face = "bold"), axis.text = element_text(size = 6),
                   legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6),
                   axis.line = element_line(linewidth = lw), axis.ticks = element_line(linewidth = lw))
pa_b <- pa + blk_theme + theme(legend.position = "right", legend.key.size = unit(3, "mm"))
pa_b$layers[[2]]$aes_params$size <- .pt2mm(5.5)
save_block(pa_b, "medip_robpca", 75, 57)
pb_b <- pb + theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1, face = "bold"),
                   legend.title = element_text(angle = 270, size = 7, face = "bold"),
                   legend.text = element_text(size = 6), legend.key.height = unit(18, "mm"))
save_block(pb_b, "medip_global_methylation_heatmap", 62, 170)
up_b <- function(df, title, col, show_legend) {
  upset(df, intersect = SETS, name = NULL, min_size = 1, width_ratio = 0.001, height_ratio = 1.1,
        sort_sets = FALSE, set_sizes = FALSE, n_intersections = 12,
        matrix = intersection_matrix(geom = geom_point(size = 0.7), segment = geom_segment(linewidth = 0.3)),
        base_annotations = list("n genes" = intersection_size(
          mapping = aes(fill = Region), text = list(size = .pt2mm(5), angle = 90, hjust = -0.15, vjust = 0.5, family = "Arial"),
          bar_number_threshold = 1e9) +
            scale_fill_manual(values = REGION_COLORS, drop = FALSE) +
            scale_y_continuous(expand = expansion(mult = c(0, 1.3)), breaks = scales::breaks_extended(3)) +
            coord_cartesian(clip = "off") +
            labs(title = title, y = "n genes", x = NULL) +
            theme(axis.title.y = element_text(size = 6, face = "bold"), legend.position = if (show_legend) "right" else "none",
                  legend.key.size = unit(2.2, "mm"), legend.background = element_blank())),
        themes = upset_modify_themes(list(intersections_matrix = theme(axis.title = element_blank(),
                                                                         text = element_text(size = 6))))) &
    theme(text = element_text(family = "Arial"), axis.text = element_text(size = 6, colour = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 7, colour = col),
          legend.title = element_text(size = 6, face = "bold"), legend.text = element_text(size = 6))
}
pc_b <- wrap_plots(up_b(hyper_df, "Hypermethylated", "#C0392B", TRUE),
                   up_b(hypo_df, "Hypomethylated", "#2471A3", FALSE), ncol = 1)
save_block(pc_b, "medip_dmr_gene_set_upset", 77, 60)
cat("[OK] MeDIP overview panels ->", OUT, "and", BLK, "\n")
