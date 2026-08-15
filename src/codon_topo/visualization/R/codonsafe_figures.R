#!/usr/bin/env Rscript
# ============================================================
# CodonSafe Meta-Analysis Figures (CS1, CS2, CS-S1)
# Clayworth & Kornilov 2026
# ggplot2 + ggpubr, 300 DPI, colorblind-friendly
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

# Source shared theme (relative path from output/codonsafe/)
source("src/codon_topo/visualization/R/theme_codon.R")

outdir <- "output/codonsafe"

# Verify required data
required <- c("syn57_s9_genuine_deviations.csv", "ostrov_segments_for_ggplot.csv")
for (f in required) {
  fp <- file.path(outdir, f)
  if (!file.exists(fp)) stop(paste("Required file missing:", fp))
}


# ═══════════════════════════════════════════════════════════════
# Figure CS1: Syn57 S9 Deviation Directionality (3 panels)
# ═══════════════════════════════════════════════════════════════
dev_data <- read.csv(file.path(outdir, "syn57_s9_genuine_deviations.csv"))

# Compute statistics dynamically
dev_nonzero <- dev_data[dev_data$delta_local != 0, ]
n_genuine   <- nrow(dev_data)
n_better    <- sum(dev_data$delta_local < 0, na.rm = TRUE)
n_worse     <- sum(dev_data$delta_local > 0, na.rm = TRUE)
n_same      <- sum(dev_data$delta_local == 0, na.rm = TRUE)
n_test      <- n_better + n_worse
binom_res   <- binom.test(n_better, n_test, p = 0.5, alternative = "greater")
binom_p     <- binom_res$p.value
ci_lo       <- binom_res$conf.int[1]
ci_hi       <- binom_res$conf.int[2]
pct_h1      <- round(100 * sum(dev_data$hamming == 1) / n_genuine, 0)
mean_ham    <- round(mean(dev_data$hamming), 2)

# Panel A: Distribution of delta_local_mismatch
fig_cs1a <- ggplot(dev_nonzero, aes(x = delta_local)) +
  geom_histogram(binwidth = 25, fill = PAL_BLUE_LT, color = "white",
                 linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = PAL_RED,
             linewidth = 0.7) +
  annotate("label", x = -300, y = Inf, vjust = 1.5, hjust = 0,
           size = ANNOT_SIZE, fill = "white", linewidth = 0.3,
           label = sprintf("n = %d genuine deviations\n%d better / %d worse / %d same\np(binomial) = %.2f\n95%% CI [%.2f, %.2f]",
                           n_genuine, n_better, n_worse, n_same,
                           binom_p, ci_lo, ci_hi)) +
  labs(
    x = expression(Delta ~ "local mismatch (Grantham)"),
    y = "Count",
    title = "Local physicochemical change"
  ) +
  theme_codon_pub()

# Panel B: Hamming distance distribution
fig_cs1b <- ggplot(dev_data, aes(x = factor(hamming))) +
  geom_bar(fill = PAL_BLUE_LT, color = "grey30", linewidth = 0.3) +
  annotate("label", x = 1, y = Inf, vjust = 1.5, hjust = 0.5,
           size = ANNOT_SIZE, fill = "white", linewidth = 0.3,
           label = sprintf("%d%% single-bit\nmean = %.2f", pct_h1, mean_ham)) +
  labs(
    x = "Hamming distance",
    y = "Count",
    title = expression("Step size in GF(2)"^6)
  ) +
  theme_codon_pub()

# Panel C: Synonymous vs nonsynonymous delta_local
dev_data$swap_type <- ifelse(dev_data$is_synonymous, "Synonymous", "Nonsynonymous")
dev_data$swap_type <- factor(dev_data$swap_type, levels = c("Synonymous", "Nonsynonymous"))

fig_cs1c <- ggplot(dev_data, aes(x = swap_type, y = delta_local)) +
  geom_boxplot(aes(fill = swap_type), alpha = 0.8, outlier.size = 1,
               color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("Synonymous" = PAL_BLUE, "Nonsynonymous" = PAL_RED_LT)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  labs(
    x = "", y = expression(Delta ~ "local mismatch"),
    title = "By swap type"
  ) +
  theme_codon_pub() +
  theme(legend.position = "none")

fig_cs1 <- ggarrange(fig_cs1a, fig_cs1b, fig_cs1c,
                     ncol = 3, nrow = 1, labels = "AUTO",
                     font.label = list(size = 14, face = "bold"))

save_figure(fig_cs1, file.path(outdir, "fig_cs1_syn57_deviations"),
            width = 15, height = 4.5)
cat("Figure CS1 saved\n")


# ═══════════════════════════════════════════════════════════════
# Figure CS2: Ostrov Segment Case-Control (3 panels)
# ═══════════════════════════════════════════════════════════════
seg_data <- read.csv(file.path(outdir, "ostrov_segments_for_ggplot.csv"))
seg_data <- seg_data[!is.na(seg_data$segment), ]

# Panel A: Recoding burden vs fitness
seg_fitness <- seg_data[!is.na(seg_data$doubling_time_deletion) &
                        !is.na(seg_data$n_essential_recoded), ]
n_fit <- nrow(seg_fitness)

cor_res <- cor.test(seg_fitness$n_essential_recoded,
                    seg_fitness$doubling_time_deletion,
                    method = "spearman", exact = FALSE)
rho_val <- round(cor_res$estimate, 2)
p_raw   <- cor_res$p.value
p_bonf  <- min(p_raw * 3, 1.0)

seg_fitness$issue_label <- ifelse(seg_fitness$has_issue, "Lethal exception", "Viable")
seg_fitness$rev_label   <- ifelse(seg_fitness$has_reversions, "Yes", "No")

fig_cs2a <- ggplot(seg_fitness,
                   aes(x = n_essential_recoded, y = doubling_time_deletion)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey50", fill = "grey90",
              linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(color = issue_label, shape = rev_label), size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Viable" = PAL_BLUE, "Lethal exception" = PAL_RED),
                     name = "Status") +
  scale_shape_manual(values = c("No" = 16, "Yes" = 17), name = "Reversions") +
  labs(
    x = "Recoded codons in essential genes",
    y = "Relative doubling time",
    title = bquote("Recoding burden vs fitness (" * italic(n) *
                   " = " * .(n_fit) * ", " * rho * " = " * .(rho_val) *
                   ", " * italic(p)[Bonf] * " = " * .(sprintf("%.3f", p_bonf)) * ")")
  ) +
  theme_codon_pub() +
  theme(legend.position = "right")

# Panel B: Problem vs normal segments
seg_data$status <- ifelse(seg_data$has_issue, "Problem", "Normal")
seg_plot_b <- seg_data[!is.na(seg_data$n_essential_recoded), ]

fig_cs2b <- ggplot(seg_plot_b, aes(x = status, y = n_essential_recoded)) +
  geom_boxplot(aes(fill = status), alpha = 0.8, outlier.size = 1,
               color = "grey30", linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, color = "grey30") +
  scale_fill_manual(values = c("Normal" = PAL_BLUE, "Problem" = PAL_RED)) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center", vjust = -0.5, size = ANNOT_SIZE) +
  labs(
    x = "", y = "Recoded codons in essential genes",
    title = "Essential gene recoding load"
  ) +
  theme_codon_pub() +
  theme(legend.position = "none")

# Panel C: Reversions per segment
seg_rev <- seg_data[!is.na(seg_data$n_reversions), ]
seg_rev$issue_label <- ifelse(seg_rev$has_issue, "Lethal exception", "Viable")

fig_cs2c <- ggplot(seg_rev, aes(x = n_recoded, y = n_reversions)) +
  geom_point(aes(color = issue_label), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Viable" = PAL_BLUE, "Lethal exception" = PAL_RED),
                     name = "Status") +
  labs(
    x = "Total recoded codons",
    y = "Codon reversions",
    title = "Design instability: reversions"
  ) +
  theme_codon_pub()

fig_cs2 <- ggarrange(fig_cs2a, fig_cs2b, fig_cs2c,
                     ncol = 3, nrow = 1, labels = "AUTO",
                     font.label = list(size = 14, face = "bold"),
                     common.legend = FALSE)

save_figure(fig_cs2, file.path(outdir, "fig_cs2_ostrov_segments"),
            width = 15, height = 5)
cat("Figure CS2 saved\n")


# ═══════════════════════════════════════════════════════════════
# Supplementary Figure CS-S1: Sensitivity Analysis
# ═══════════════════════════════════════════════════════════════
sens_file <- file.path(outdir, "syn57_s9_sensitivity_table.csv")
if (file.exists(sens_file)) {
  sens <- read.csv(sens_file)
  sens <- sens[sens$cutoff != "all", ]
  sens$cutoff_num <- as.numeric(gsub("[^0-9]", "", sens$cutoff))

  fig_css1 <- ggplot(sens, aes(x = cutoff_num, y = binomial_p)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60",
               linewidth = 0.4) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = PAL_RED,
               linewidth = 0.5) +
    geom_line(color = PAL_BLUE, linewidth = 0.9) +
    geom_point(color = PAL_BLUE, size = 3, shape = 16) +
    annotate("text", x = max(sens$cutoff_num), y = 0.08,
             hjust = 1, size = ANNOT_SIZE, color = PAL_RED,
             label = expression(italic(p) ~ "= 0.05 threshold")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Max deviations per gene (filter cutoff)",
      y = expression("Binomial " * italic(p) * "-value (better > 0.5)"),
      title = "Sensitivity Analysis: Filtering Threshold vs Null Result"
    ) +
    theme_codon_pub()

  save_figure(fig_css1, file.path(outdir, "fig_css1_sensitivity"))
  cat("Figure CS-S1 (sensitivity) saved\n")
}

cat("All CodonSafe figures generated.\n")
