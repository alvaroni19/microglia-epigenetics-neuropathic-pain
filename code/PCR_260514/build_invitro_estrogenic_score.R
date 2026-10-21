# =============================================================================
# In vitro estrogen-related genes (Esr1, Esr2, Gper1, Cyp19a1) in primary
# microglia: directed estrogenic score per replicate and exploratory PCA.
#
# Input : data/PCR_260514/invitro_qpcr_relative_expression.xlsx
#         ("<gene> 1 way without E2" sheets; six columns = male control, male
#         CSF1, male CSF1+E2, female control, female CSF1, female CSF1+E2;
#         rows = biological replicates)
# Output: generated_outputs/PCR_260514/
#         invitro_estrogenic_score_per_replicate.csv
#         invitro_estrogen_genes_pca_{scores,loadings,variance}.csv
#
# Score: each gene is z-scored across all replicates; the score of a replicate
# is the mean z of the four genes (full), of the three receptors (reduced),
# of Esr2 + Cyp19a1, or of Esr2 + Gper1. "delta_*" columns subtract the median
# score of the control replicates of the same sex.
# PCA: stats::prcomp(center = TRUE, scale. = TRUE) on replicates with all four
# genes measured.
# =============================================================================
suppressMessages({ library(readxl); library(dplyr); library(tidyr) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), "..", ".."), mustWork = TRUE)
WB <- file.path(ROOT, "data", "PCR_260514", "invitro_qpcr_relative_expression.xlsx")
OUT <- file.path(ROOT, "generated_outputs", "PCR_260514")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GENES <- c("ERa", "ERb", "GPER", "ARO")            # sheet keys: Esr1, Esr2, Gper1, Cyp19a1
GENE_LABELS <- c(ERa = "ER\u03b1", ERb = "ER\u03b2", GPER = "GPER", ARO = "ARO")
GROUP_IDS <- c("male_ctrl", "male_csf1", "male_csf1e2", "fem_ctrl", "fem_csf1", "fem_csf1e2")
GROUP_ORDER <- c("Male Control", "Male CSF1", "Male CSF1+E2",
                 "Female Control", "Female CSF1", "Female CSF1+E2")
GROUP_LABELS <- setNames(GROUP_ORDER, GROUP_IDS)

long <- bind_rows(lapply(GENES, function(gene) {
  df <- read_excel(WB, sheet = paste(gene, "1 way without E2"), col_names = FALSE, skip = 1,
                   .name_repair = "minimal")
  colnames(df) <- GROUP_IDS[seq_len(ncol(df))]
  df %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.))), gene = gene, rep = row_number()) %>%
    pivot_longer(-c(gene, rep), names_to = "group_id", values_to = "value") %>%
    filter(!is.na(value))
})) %>%
  mutate(gene = factor(gene, levels = GENES),
         group_label = factor(GROUP_LABELS[group_id], levels = GROUP_ORDER))

# ---------------------------------------------------------------- score
score <- long %>%
  group_by(gene) %>% mutate(z = as.numeric(scale(value))) %>% ungroup() %>%
  group_by(group_label, rep) %>%
  summarise(ESI_full = mean(z, na.rm = TRUE),
            ESI_reduced = mean(z[as.character(gene) != "ARO"], na.rm = TRUE),
            ESI_erb_aro = mean(z[as.character(gene) %in% c("ERb", "ARO")], na.rm = TRUE),
            ESI_erb_gper = mean(z[as.character(gene) %in% c("ERb", "GPER")], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(group_chr = as.character(group_label),
         sex = if_else(grepl("^Male", group_chr), "Male", "Female"),
         treatment = case_when(grepl("Control", group_chr) ~ "Control",
                               grepl("CSF1\\+E2", group_chr) ~ "CSF1+E2",
                               grepl("CSF1", group_chr) ~ "CSF1"),
         treatment = factor(treatment, levels = c("Control", "CSF1", "CSF1+E2"))) %>%
  group_by(sex) %>%
  mutate(delta_full = ESI_full - median(ESI_full[treatment == "Control"], na.rm = TRUE),
         delta_red = ESI_reduced - median(ESI_reduced[treatment == "Control"], na.rm = TRUE),
         delta_ea = ESI_erb_aro - median(ESI_erb_aro[treatment == "Control"], na.rm = TRUE),
         delta_eg = ESI_erb_gper - median(ESI_erb_gper[treatment == "Control"], na.rm = TRUE)) %>%
  ungroup()
write.csv(score, file.path(OUT, "invitro_estrogenic_score_per_replicate.csv"), row.names = FALSE)

# ---------------------------------------------------------------- PCA
wide <- long %>%
  mutate(group_chr = as.character(group_label), sample_id = paste(group_id, rep, sep = "_"),
         sex = if_else(grepl("^Male", group_chr), "Male", "Female"),
         treatment = case_when(grepl("CSF1\\+E2$", group_chr) ~ "CSF1+E2",
                               grepl("CSF1$", group_chr) ~ "CSF1", TRUE ~ "Control"),
         gene = as.character(gene)) %>%
  select(sample_id, group_chr, sex, treatment, rep, gene, value) %>%
  pivot_wider(names_from = gene, values_from = value) %>%
  mutate(group_label = factor(group_chr, levels = GROUP_ORDER),
         sex = factor(sex, levels = c("Male", "Female")),
         treatment = factor(treatment, levels = c("Control", "CSF1", "CSF1+E2"))) %>%
  select(sample_id, group_label, sex, treatment, rep, all_of(GENES)) %>%
  filter(if_all(all_of(GENES), ~ is.finite(.)))
mat <- as.matrix(wide[, GENES]); rownames(mat) <- wide$sample_id
pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
var_expl <- pca$sdev^2 / sum(pca$sdev^2)
k <- seq_len(min(4, ncol(pca$x)))
scores <- as.data.frame(pca$x[, k, drop = FALSE]) %>%
  mutate(sample_id = rownames(pca$x)) %>%
  left_join(wide %>% select(sample_id, group_label, sex, treatment, rep), by = "sample_id")
loadings <- as.data.frame(pca$rotation[, k, drop = FALSE]) %>%
  mutate(gene = rownames(pca$rotation), gene_label = GENE_LABELS[gene]) %>%
  select(gene, gene_label, everything())
variance <- data.frame(PC = paste0("PC", seq_along(var_expl)), variance_explained = var_expl,
                       variance_explained_pct = 100 * var_expl)
write.csv(scores, file.path(OUT, "invitro_estrogen_genes_pca_scores.csv"), row.names = FALSE)
write.csv(loadings, file.path(OUT, "invitro_estrogen_genes_pca_loadings.csv"), row.names = FALSE)
write.csv(variance, file.path(OUT, "invitro_estrogen_genes_pca_variance.csv"), row.names = FALSE)
cat(sprintf("[OK] in vitro estrogenic score (%d replicates) and PCA (PC1 %.1f%%, PC2 %.1f%%)\n",
            nrow(score), 100 * var_expl[1], 100 * var_expl[2]))
