# Global methylation heatmap of ROBPCA-eligible significant genes.
suppressMessages({ library(ggplot2); library(dplyr); library(tibble); library(tidyr) })

current_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                         error = function(e) NA_character_)
if (is.na(current_file) || !nzchar(current_file)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  current_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
}
SCRIPT_DIR <- dirname(normalizePath(current_file, winslash = "/", mustWork = FALSE))
BUNDLE_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/")
OUT <- file.path(BUNDLE_ROOT, "generated_outputs", "MeDIP")

table_dir <- file.path(BUNDLE_ROOT, "data", "MeDIP", "workspace_and_tables")
matrix_table <- utils::read.csv(
  file.path(table_dir, "significant_gene_rms_matrix_ws400.csv"),
  check.names = FALSE
)
mat <- as.matrix(matrix_table[, -1, drop = FALSE])
rownames(mat) <- matrix_table$gene
robpca_eligible <- apply(mat, 1, function(values) all(is.finite(values))) &
  apply(mat, 1, function(values) stats::mad(values, na.rm = TRUE) > 0)
informative_symbol <- !grepl("^(AA|AC)", rownames(mat), ignore.case = TRUE)
mat <- mat[robpca_eligible & informative_symbol, , drop = FALSE]
if (nrow(mat) < 2L) {
  stop("The heatmap requires at least two complete ROBPCA-eligible genes.")
}
select_top_variable_rows <- function(values, n_top = 80L) {
  if (nrow(values) <= n_top) return(values)
  row_variance <- apply(values, 1, stats::var, na.rm = TRUE)
  values[order(row_variance, decreasing = TRUE)[seq_len(n_top)], , drop = FALSE]
}
row_zscore <- function(values) {
  scaled <- t(scale(t(values)))
  scaled[!is.finite(scaled)] <- 0
  scaled
}
mat <- row_zscore(select_top_variable_rows(mat, 80L))
message(
  "Heatmap selected 80 top-variable genes from ", nrow(matrix_table),
  " significant genes after ROBPCA eligibility and AA/AC-prefix filtering."
)

cluster_order <- function(m, rows = TRUE) {
  x <- if (rows) m else t(m)
  if (nrow(x) < 2) return(seq_len(nrow(x)))
  cc <- suppressWarnings(stats::cor(t(x), use = "pairwise.complete.obs"))
  cc[!is.finite(cc)] <- 0; diag(cc) <- 1
  stats::hclust(stats::as.dist(1 - cc), method = "average")$order
}
mat <- mat[cluster_order(mat, TRUE), cluster_order(mat, FALSE), drop = FALSE]
nr <- nrow(mat); nc <- ncol(mat)

dfp <- as.data.frame(mat) %>% tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
  mutate(gene = factor(gene, levels = rev(rownames(mat))),
         sample = factor(sample, levels = colnames(mat)))

# Relabel x-axis groups as M-C / F-C / M-CSF1 / F-CSF1 (keeping the replicate number).
relab_sample <- function(s) {
  s <- sub("^MCSF1", "M-CSF1", s); s <- sub("^FCSF1", "F-CSF1", s)
  s <- sub("^MC", "M-C", s);       s <- sub("^FC", "F-C", s)
  s
}

p <- ggplot(dfp, aes(sample, gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                       na.value = "grey90", name = "Row z-score") +
  scale_x_discrete(expand = c(0, 0), labels = relab_sample) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Global methylation — significant genes (RMS)", x = NULL, y = NULL) +
  theme_classic(base_size = 6) +
  theme(text = element_text(family = "sans", colour = "black"),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 2)),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 5, face = "italic", colour = "black"),
        axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_text(size = 6.5, face = "bold"), legend.text = element_text(size = 6),
        legend.key.height = unit(6, "mm"), legend.key.width = unit(2.5, "mm"),
        plot.margin = margin(2, 3, 1, 1))

W <- 92; H <- max(120, nr * 1.7 + 18)
base <- file.path(OUT, "medip_global_methylation_heatmap")
ragg::agg_tiff(paste0(base, ".tiff"), width = W, height = H, units = "mm",
               res = 600, background = "white")
print(p); dev.off()
# Vector outputs: PDF (text-as-text, feeds the editable-PPTX pipeline) + SVG.
grDevices::cairo_pdf(paste0(base, ".pdf"), width = W / 25.4, height = H / 25.4, bg = "white")
print(p); dev.off()
svglite::svglite(paste0(base, ".svg"), width = W / 25.4, height = H / 25.4, bg = "white")
print(p); dev.off()
cat(sprintf("[OK] medip_global_methylation_heatmap (%d genes x %d samples)\n", nr, nc))
