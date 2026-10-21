# =============================================================================
# Multiple Factor Analysis (MFA) of brainstem pain + molecular markers,
# rebuilt at true mm size (5-7 pt Arial). Faithful to
# the original exploratory MFA: same variable blocks, missMDA::imputeMFA(seed=123)
# and FactoMineR::MFA. Phase colours match the shared palette.
# =============================================================================
suppressMessages({ library(dplyr); library(ggplot2); library(FactoMineR); library(vegan) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}
BUNDLE_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/")
OUT <- file.path(BUNDLE_ROOT, "generated_outputs", "PCR_MFA")
.pt2mm <- function(pt) pt / 2.845276

correl <- read.csv(file.path(OUT, "invivo_qpcr_mfa_input.csv"), check.names = FALSE)
names(correl)[names(correl) == "phase"] <- "phase"
names(correl)[names(correl) == "condition"] <- "condition"
names(correl)[names(correl) == "pain_016_delta"] <- "pain016"
phase_order <- c("Male","Proestrus","Estrus","Diestrus"); condition_order <- c("Control","CCI")
correl$phase <- factor(correl$phase, levels = phase_order)
correl$condition <- factor(correl$condition, levels = condition_order)
correl$group <- interaction(correl$phase, correl$condition, drop = TRUE)

pain_var <- "pain016"
microglia_vars <- c("Slc12a5","Bdnf","Csf1r","Tlr4")
hormone_vars <- c("Esr1","Esr2","Gper1","Cyp19a1")
methylation_vars <- c("Dnmt1","Dnmt3a","Dnmt3b","Gadd45a","Gadd45b","Gadd45g")
# Active blocks as in the original exploratory MFA: pain, microglial markers,
# estrogen-related genes and DNA-methylation machinery.
mfa_vars <- c(pain_var, microglia_vars, hormone_vars, methylation_vars)
mfa_group <- c(1, length(microglia_vars), length(hormone_vars), length(methylation_vars))
mfa_type <- rep("s", length(mfa_group))

correlmfa <- as.data.frame(correl[, mfa_vars])
rownames(correlmfa) <- make.unique(paste(correl$group, correl$id, sep = "_"))
if (anyNA(correlmfa) && requireNamespace("missMDA", quietly = TRUE)) {
  set.seed(123)
  correlmfa_complete <- as.data.frame(missMDA::imputeMFA(correlmfa, group = mfa_group, type = mfa_type, ncp = 2, seed = 123)$completeObs)
} else {
  correlmfa_complete <- correlmfa
  for (v in mfa_vars) if (anyNA(correlmfa_complete[[v]])) correlmfa_complete[[v]][is.na(correlmfa_complete[[v]])] <- mean(correlmfa_complete[[v]], na.rm = TRUE)
}
res1 <- FactoMineR::MFA(correlmfa_complete, group = mfa_group, type = mfa_type,
                        name.group = c("pain","microglia","hormones","methylation"), ncp = 5, graph = FALSE)
eig <- res1$eig
var1 <- eig[1, 2]; var2 <- eig[2, 2]

df <- as.data.frame(res1$ind$coord); colnames(df)[1:2] <- c("Dim1","Dim2")
df$animal_id <- correl$animal_id
df$phase <- factor(correl$phase, levels = phase_order)
df$condition <- factor(correl$condition, levels = condition_order)

# Export the individual MFA coordinates + variance so the panel can be re-plotted
# in Python (editable SVG/PPTX) with proestrus drawn on top (F4E z-order).
write.csv(data.frame(animal_id = df$animal_id, Dim1 = df$Dim1, Dim2 = df$Dim2, phase = df$phase, condition = df$condition),
          file.path(OUT, "invivo_mfa_coords.csv"), row.names = FALSE)
write.csv(data.frame(var1 = var1, var2 = var2),
          file.path(OUT, "invivo_mfa_eig.csv"), row.names = FALSE)

# Quantify group separation in the minimum number of MFA dimensions explaining
# at least 70% cumulative variance, matching the documented analysis.
cumulative_variance <- eig[, 3]
k_dim <- which(cumulative_variance >= 70)[1]
if (is.na(k_dim)) k_dim <- min(5, ncol(res1$ind$coord))
k_dim <- max(2, min(k_dim, ncol(res1$ind$coord)))
mfa_scores <- as.data.frame(res1$ind$coord[, seq_len(k_dim), drop = FALSE])
metadata <- data.frame(
  phase = correl$phase,
  condition = correl$condition,
  group = correl$group
)
dist_mfa <- dist(mfa_scores, method = "euclidean")
set.seed(123)
permanova <- vegan::adonis2(
  dist_mfa ~ phase * condition,
  data = metadata,
  permutations = 9999,
  by = "terms"
)
permanova_df <- as.data.frame(permanova)
permanova_df$term <- rownames(permanova_df)
permanova_df <- permanova_df[, c("term", setdiff(names(permanova_df), "term"))]
write.csv(permanova_df, file.path(OUT, "mfa_permanova_terms.csv"), row.names = FALSE)

extract_betadisper <- function(group_variable) {
  bd <- vegan::betadisper(dist_mfa, metadata[[group_variable]], bias.adjust = TRUE)
  set.seed(123)
  test <- vegan::permutest(bd, permutations = 9999)
  tab <- as.data.frame(test$tab)
  data.frame(
    variable = group_variable,
    F_value = tab[1, "F"],
    p_value = tab[1, "Pr(>F)"]
  )
}
dispersion <- do.call(
  rbind,
  lapply(c("phase", "condition", "group"), extract_betadisper)
)
write.csv(dispersion, file.path(OUT, "mfa_betadisper_results.csv"), row.names = FALSE)

phase_colors <- c("Male" = "#4575B4", "Proestrus" = "#E7298A", "Estrus" = "#F4C95D", "Diestrus" = "#E8826B")
condition_shapes <- c("Control" = 21, "CCI" = 24)
condition_labels <- c("Control" = "Control", "CCI" = "CCI-IoN")

p <- ggplot(df, aes(Dim1, Dim2)) +
  stat_ellipse(aes(group = phase, fill = phase), type = "t", level = 0.95, geom = "polygon", alpha = 0.18, color = NA) +
  stat_ellipse(aes(group = phase, color = phase), type = "t", level = 0.95, linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_point(aes(shape = condition, fill = phase), size = 1.9, color = "black", stroke = 0.4) +
  scale_color_manual(values = phase_colors, name = "Biological phase") +
  scale_fill_manual(values = phase_colors, name = "Biological phase") +
  scale_shape_manual(values = condition_shapes, labels = condition_labels, name = "Condition") +
  labs(title = "Multiple Factor Analysis (MFA)",
       x = sprintf("Dimension 1 (%.1f%%)", var1), y = sprintf("Dimension 2 (%.1f%%)", var2)) +
  theme_classic(base_size = 6) +
  theme(text = element_text(family = "sans", colour = "black"),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 2)),
        axis.title = element_text(size = 6.5, face = "bold"), axis.text = element_text(size = 5.5, colour = "black"),
        axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        legend.title = element_text(size = 6, face = "bold"), legend.text = element_text(size = 5.5),
        legend.position = "right", legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
        plot.margin = margin(2, 2, 1, 1)) +
  guides(fill = guide_legend(override.aes = list(size = 2, shape = 21)),
         color = "none")

W <- 96; H <- 66
ragg::agg_tiff(file.path(OUT, "invivo_mfa.tiff"), width = W, height = H, units = "mm", res = 600, background = "white")
print(p); dev.off()
ggplot2::ggsave(file.path(OUT, "invivo_mfa.pdf"), p, width = W, height = H, units = "mm", device = cairo_pdf)
cat(sprintf(
  "[OK] invivo_mfa  Dim1=%.1f%% Dim2=%.1f%%; %d retained markers; k=%d\n",
  var1, var2, length(microglia_vars) + length(hormone_vars) + length(methylation_vars), k_dim
))
