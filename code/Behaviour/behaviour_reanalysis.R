# Reanalysis of mechanical allodynia (von Frey 0.16 g)
# Outcome: Pain_norm = test score (post-CCI) - mean baseline (3 days pre-CCI)
# Input : data/Allodynia_0_16/allodynia_change_scores.xlsx (one row per animal)
# Output: generated_outputs/Behaviour_reanalysis/ (statistics workbook, plotting tables)
# Run from the repository root:  Rscript code/Behaviour/behaviour_reanalysis.R
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(emmeans)
  library(sandwich); library(lmtest); library(car); library(writexl); library(ggplot2); library(patchwork)
})
set.seed(20261014)
options(contrasts = c("contr.sum", "contr.poly"))  # required for valid type-III tests
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "."
OUT <- file.path(ROOT, "generated_outputs", "Behaviour_reanalysis"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
d <- read_excel(file.path(ROOT, "data", "Allodynia_0_16", "allodynia_change_scores.xlsx"))
stopifnot(nrow(d) == 105, all(abs(d$Test016 - d$BL016 - d$Pain_norm) < 1e-8))

d <- d %>% mutate(
  Group = factor(unname(c(M = "Male", Proestrus = "Proestrus", Estrus = "Estrus", Diestrus = "Diestrus")[phase_surgery]),
                 levels = c("Male", "Proestrus", "Estrus", "Diestrus")),
  Condition = factor(Surgery, levels = c("Control", "CCI")),
  Delta = Pain_norm)

# ---------------------------------------------------------------- descriptives
desc <- d %>% group_by(Condition, Group) %>% summarise(
  n = n(), BL_median = median(BL016), Test_median = median(Test016),
  Delta_mean = mean(Delta), Delta_sd = sd(Delta), Delta_median = median(Delta),
  Delta_IQR = IQR(Delta), .groups = "drop")

# paired pre/post: DESCRIPTIVE only (kept for transparency)
paired <- d %>% group_by(Condition, Group) %>% summarise(
  n = n(),
  median_BL = median(BL016), median_Test = median(Test016),
  median_paired_diff = median(Delta),
  wilcoxon_p_descriptive = suppressWarnings(wilcox.test(Test016, BL016, paired = TRUE, exact = FALSE)$p.value),
  .groups = "drop")

# ---------------------------------------------------------------- primary model
# Linear model on change score, Condition x Group (Group = Male + 3 estrous phases at surgery)
# Heteroscedasticity-consistent (HC3) standard errors for all inference.
m <- lm(Delta ~ Condition * Group, data = d)
V <- vcovHC(m, type = "HC3")
anova_hc3 <- Anova(m, type = 3, vcov. = V)
# Females only: Condition x estrous phase interaction (key hypothesis)
df_f <- droplevels(filter(d, Group != "Male"))
mf <- lm(Delta ~ Condition * Group, data = df_f, contrasts = list(Condition = contr.sum, Group = contr.sum))
anova_f <- Anova(mf, type = 3, vcov. = vcovHC(mf, type = "HC3"))

emm <- emmeans(m, ~ Condition * Group, vcov. = V)
# (a) injury effect (CCI - Control) within each group  = difference-in-differences vs baseline
inj <- contrast(emm, method = "revpairwise", by = "Group", adjust = "none")
inj_df <- summary(inj, infer = TRUE) %>% as.data.frame()
inj_df$p_holm <- p.adjust(inj_df$p.value, "holm")
# (b) direct comparison of phases within CCI animals (change score)
cci_ph <- contrast(emmeans(m, ~ Group | Condition, vcov. = V, at = list()), method = "pairwise", adjust = "none")
cci_ph_df <- summary(cci_ph, infer = TRUE) %>% as.data.frame() %>% filter(Condition == "CCI")
fem_rows <- !grepl("Male", cci_ph_df$contrast)
cci_ph_df$p_holm_females <- NA; cci_ph_df$p_holm_females[fem_rows] <- p.adjust(cci_ph_df$p.value[fem_rows], "holm")
# (c) interaction contrasts: does the injury effect differ between phases? (DiD of DiD)
int_c <- contrast(emm, interaction = c(Condition = "revpairwise", Group = "pairwise"), adjust = "none")
int_df <- summary(int_c, infer = TRUE) %>% as.data.frame()
fem_int <- !grepl("Male", int_df$Group_pairwise)
int_df$p_holm_females <- NA; int_df$p_holm_females[fem_int] <- p.adjust(int_df$p.value[fem_int], "holm")

# ---------------------------------------------------------------- sensitivity analyses
# Sensitivity analysis 1: ANCOVA, Test ~ Baseline + Condition*Group
m_anc <- lm(Test016 ~ BL016 + Condition * Group, data = d)
emm_anc <- emmeans(m_anc, ~ Condition * Group, vcov. = vcovHC(m_anc, type = "HC3"))
sens1_inj  <- summary(contrast(emm_anc, "revpairwise", by = "Group", adjust = "holm"), infer = TRUE) %>% as.data.frame()
sens1_ph   <- summary(contrast(emmeans(m_anc, ~ Group | Condition, vcov. = vcovHC(m_anc, type = "HC3")), "pairwise", adjust = "none"),
                   infer = TRUE) %>% as.data.frame() %>% filter(Condition == "CCI", !grepl("Male", contrast))
sens1_ph$p_holm <- p.adjust(sens1_ph$p.value, "holm")
sens1_int  <- Anova(lm(Test016 ~ BL016 + Condition * Group, data = df_f,
                    contrasts = list(Condition = contr.sum, Group = contr.sum)), type = 3,
                 vcov. = vcovHC(lm(Test016 ~ BL016 + Condition * Group, data = df_f,
                                   contrasts = list(Condition = contr.sum, Group = contr.sum)), type = "HC3"))

# Sensitivity analysis 2: distribution-free permutation tests (10,000 permutations)
B <- 10000
perm_diff <- function(x, g, a, b) {         # difference in means, two groups, label permutation
  obs <- mean(x[g == a]) - mean(x[g == b]); xx <- x[g %in% c(a, b)]; gg <- g[g %in% c(a, b)]
  null <- replicate(B, { p <- sample(gg); mean(xx[p == a]) - mean(xx[p == b]) })
  c(diff = obs, p_perm = (sum(abs(null) >= abs(obs) - 1e-12) + 1) / (B + 1))
}
cci_f <- filter(d, Condition == "CCI", Group != "Male")
pairs <- list(c("Proestrus", "Estrus"), c("Proestrus", "Diestrus"), c("Estrus", "Diestrus"))
sens2 <- bind_rows(lapply(pairs, function(p) {
  r <- perm_diff(cci_f$Delta, as.character(cci_f$Group), p[1], p[2])
  data.frame(contrast = paste(p[1], "-", p[2]), diff_mean_delta = r["diff"], p_perm = r["p_perm"])
}))
sens2$p_holm <- p.adjust(sens2$p_perm, "holm")
# permutation test of Condition x Phase interaction in females (F statistic, permuting residuals of additive model: Freedman-Lane)
m0 <- lm(Delta ~ Condition + Group, data = df_f); F_obs <- anova(m0, mf)$F[2]
fit0 <- fitted(m0); res0 <- resid(m0)
F_null <- replicate(B, { y <- fit0 + sample(res0); anova(lm(y ~ Condition + Group, data = df_f),
                                                          lm(y ~ Condition * Group, data = df_f))$F[2] })
sens2_int <- data.frame(test = "Condition x estrous phase (females), Freedman-Lane permutation",
                     F_obs = F_obs, p_perm = (sum(F_null >= F_obs) + 1) / (B + 1))
# Sensitivity analysis 3: Kruskal-Wallis on change score among CCI females + pairwise Wilcoxon (Holm)
sens3_kw <- kruskal.test(Delta ~ Group, data = droplevels(cci_f))
sens3_pw <- pairwise.wilcox.test(cci_f$Delta, droplevels(cci_f$Group), p.adjust.method = "holm", exact = FALSE)

# ---------------------------------------------------------------- export
tidy_anova <- function(a) data.frame(term = rownames(a), as.data.frame(a), row.names = NULL, check.names = FALSE)
write_xlsx(list(
  README = data.frame(item = c("Outcome", "Primary model", "Inference", "Groups", "Paired tests", "Data source", "Date"),
    value = c("Delta = von Frey 0.16 g test score minus 3-day mean baseline (Pain_norm)",
              "lm(Delta ~ Condition * Group); Group = Male / estrous phase at surgery",
              "HC3 robust SE; Holm within each contrast family (female phase contrasts)",
              "Control and CCI animals; estrous phase recorded on surgery day",
              "Pre/post Wilcoxon kept as DESCRIPTIVE only",
              "Copy of allodynia_change_scores.xlsx (bundle unchanged)", as.character(Sys.Date()))),
  descriptives = desc, paired_descriptive = paired,
  ANOVA_all_HC3 = tidy_anova(anova_hc3), ANOVA_females_HC3 = tidy_anova(anova_f),
  injury_effect_by_group = inj_df, CCI_phase_comparisons = cci_ph_df, interaction_contrasts = int_df,
  Sens1_ANCOVA_injury = sens1_inj, Sens1_ANCOVA_CCI_phase = sens1_ph, Sens1_ANCOVA_int_females = tidy_anova(sens1_int),
  Sens2_permutation_CCI_phase = sens2, Sens2_permutation_interaction = sens2_int,
  Sens3_KW_CCI_females = data.frame(H = unname(sens3_kw$statistic), df = unname(sens3_kw$parameter), p = sens3_kw$p.value),
  Sens3_pairwise_wilcox_holm = as.data.frame(as.table(sens3_pw$p.value)) %>% filter(!is.na(Freq)) %>% rename(p_holm = Freq)
), file.path(OUT, "behaviour_reanalysis_stats.xlsx"))

# ---------------------------------------------------------------- figure (exploratory preview)
cols <- c(Control = "#999999", CCI = "#CC79A7")
emm_df <- as.data.frame(summary(emm))
p1 <- ggplot(d, aes(Group, Delta, colour = Condition)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .6, seed = 1), alpha = .6, size = 1.6) +
  geom_pointrange(data = emm_df, aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = Condition),
                  position = position_dodge(.6), colour = "black", shape = 21, size = .4) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(x = "Group (estrous phase at surgery)", y = expression(Delta~"pain score (test - baseline)"),
       title = "Change from baseline, von Frey 0.16 g") + theme_classic(base_size = 10)
inj_plot <- inj_df %>% mutate(Group = factor(Group, levels = levels(d$Group)))
p2 <- ggplot(inj_plot, aes(Group, estimate)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), colour = "#CC79A7") +
  labs(x = NULL, y = "CCI - Control (difference in Delta, 95% CI)", title = "Injury effect by group") +
  theme_classic(base_size = 10)
fig <- p1 + p2 + plot_layout(widths = c(1.6, 1))
ggsave(file.path(OUT, "behaviour_reanalysis_preview.png"), fig, width = 9, height = 3.8, dpi = 200)
saveRDS(list(desc = desc, paired = paired, anova_hc3 = anova_hc3, anova_f = anova_f, inj = inj_df,
             cci_ph = cci_ph_df, int = int_df, sens1_inj = sens1_inj, sens1_ph = sens1_ph, sens1_int = sens1_int, sens2 = sens2,
             sens2_int = sens2_int, sens3_kw = sens3_kw, sens3_pw = sens3_pw), file.path(OUT, "results.rds"))
cat("DONE\n")

# ---------------------------------------------------------------- tables for the change-score violins
write.csv(d %>% transmute(n, Condition, Group, BL016, Test016, Delta),
          file.path(OUT, "plot_data.csv"), row.names = FALSE)
emm_hc3 <- emmeans(m, ~ Condition * Group, vcov. = V)
write.csv(as.data.frame(summary(emm_hc3)), file.path(OUT, "emm_HC3.csv"), row.names = FALSE)
