# Robust PCA (ROBPCA) of the significant-gene RMS matrix.
suppressMessages({ library(ggplot2); library(dplyr); library(tibble) })

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
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

.pt2mm <- function(pt) pt / 2.845276
MM <- 1 / 25.4

# Unified sex x treatment palette (Male=vermillion tints, Female=blue tints).
GROUP_COLORS <- c(
  "FC" = "#0072B2", "FCSF1" = "#7FB8D8", "MC" = "#D55E00", "MCSF1" = "#EAAE7F",
  "Female control" = "#0072B2", "Female CSF" = "#7FB8D8", "Female CSF1" = "#7FB8D8",
  "Male control" = "#D55E00", "Male CSF" = "#EAAE7F", "Male CSF1" = "#EAAE7F"
)

table_dir <- file.path(BUNDLE_ROOT, "data", "MeDIP", "workspace_and_tables")
matrix_table <- utils::read.csv(
  file.path(table_dir, "significant_gene_rms_matrix_ws400.csv"),
  check.names = FALSE
)
pca_sig_mat <- as.matrix(matrix_table[, -1, drop = FALSE])
rownames(pca_sig_mat) <- matrix_table$gene
sample_metadata <- utils::read.csv(
  file.path(table_dir, "medip_sample_metadata.csv"),
  stringsAsFactors = FALSE
)

finite_rows <- apply(pca_sig_mat, 1, function(values) all(is.finite(values)))
positive_robust_scale <- apply(
  pca_sig_mat,
  1,
  function(values) stats::mad(values, na.rm = TRUE) > 0
)
informative_symbol <- !grepl("^(AA|AC)", rownames(pca_sig_mat), ignore.case = TRUE)
mat_use <- pca_sig_mat[
  finite_rows & positive_robust_scale & informative_symbol,
  ,
  drop = FALSE
]
if (nrow(mat_use) < 2L) {
  stop("ROBPCA requires at least two complete genes with positive robust scale.")
}
message(
  "ROBPCA retained ", nrow(mat_use), " of ", nrow(pca_sig_mat),
  " genes after complete-case, positive-MAD and AA/AC-prefix filtering."
)
meta <- sample_metadata %>%
  dplyr::filter(as.character(sample_label) %in% colnames(mat_use)) %>%
  dplyr::mutate(sample_label = as.character(sample_label))

stopifnot(requireNamespace("rrcov", quietly = TRUE))
stopifnot(requireNamespace("ggrepel", quietly = TRUE))
set.seed(1)
rob <- rrcov::PcaHubert(t(mat_use), k = min(5, ncol(mat_use) - 1, nrow(mat_use) - 1), scale = TRUE)
sc <- as.data.frame(rob@scores[, 1:2, drop = FALSE]); names(sc) <- c("RPC1", "RPC2")
sc <- sc %>% tibble::rownames_to_column("sample_label") %>% dplyr::left_join(meta, by = "sample_label")

# normalise group labels to colour keys
sc$group <- as.character(sc$group)

# Minimum text size for axis / legend / in-plot labels is 7 pt.
base <- 7
p <- ggplot(sc, aes(RPC1, RPC2, color = group, shape = condition, label = sample_label)) +
  geom_point(size = 1.7, stroke = 0.5) +
  scale_color_manual(values = GROUP_COLORS, name = "Group", na.value = "grey50") +
  scale_shape_manual(values = c(Control = 16, CSF1 = 17), name = "Treatment", na.value = 15) +
  labs(title = "ROBPCA — significant genes", x = "Robust PC1", y = "Robust PC2") +
  theme_classic(base_size = base) +
  theme(
    text = element_text(family = "sans", colour = "black"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    axis.title = element_text(size = 7, face = "bold", colour = "black"),
    axis.text = element_text(size = 7, colour = "black"),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.position = "bottom", legend.box = "horizontal",
    legend.key.size = unit(2.4, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(2, 3, 1, 1)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.7), nrow = 2, order = 1),
         shape = guide_legend(nrow = 2, order = 2))

p <- p + ggrepel::geom_text_repel(size = .pt2mm(7), fontface = "plain",
                                  family = "sans", show.legend = FALSE,
                                  max.overlaps = Inf, min.segment.length = 0,
                                  segment.size = 0.3, seed = 1, box.padding = 0.25)

W <- 82; H <- 84
ragg::agg_tiff(file.path(OUT, "medip_robpca_samples.tiff"),
               width = W, height = H, units = "mm", res = 600, background = "white")
print(p); dev.off()
ggplot2::ggsave(file.path(OUT, "medip_robpca_samples.pdf"), p,
                width = W, height = H, units = "mm", device = cairo_pdf)
ggplot2::ggsave(file.path(OUT, "medip_robpca_samples.svg"), p,
                width = W, height = H, units = "mm", device = svglite::svglite)
cat("[OK] medip_robpca_samples (tiff+pdf+svg)\n")
