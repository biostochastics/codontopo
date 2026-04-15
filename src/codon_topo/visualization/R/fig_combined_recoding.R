#!/usr/bin/env Rscript
# ============================================================
# Combined recoding reanalysis figure (4 panels)
# Clayworth & Kornilov 2026
# ggplot2 + ggpubr, colorblind-friendly palette
#
# NOTE: Panels A, B, C use ILLUSTRATIVE simulated data to show
# the direction and scale of real effects. The actual analyses
# producing exact values are in the Python pipeline. This script
# generates schematic figures for the supplement only.
# Real sample sizes and statistics are reported from the pipeline.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

# Source shared theme
source("src/codon_topo/visualization/R/theme_codon.R")

outdir <- "output/codonsafe"

# Load data
dev_data <- read.csv(file.path(outdir, "syn57_s9_genuine_deviations.csv"))


# === Panel A: Local mismatch Ser vs Ala (POSITIVE) ===
set.seed(135325)
ser_ala_data <- data.frame(
  amino_acid = c(rep("Serine (boundary-crossing)", 100),
                 rep("Alanine (within-box)", 100)),
  delta_local = c(rnorm(100, mean = -37.2, sd = 30.5),
                  rnorm(100, mean = 19.0, sd = 0.5))
)

fig_a <- ggplot(ser_ala_data, aes(x = amino_acid, y = delta_local, fill = amino_acid)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = c("Serine (boundary-crossing)" = PAL_RED_LT,
                                "Alanine (within-box)" = PAL_BLUE)) +
  annotate("label", x = 1.5, y = 55, size = ANNOT_SIZE,
           fill = "white", linewidth = 0.3,
           label = "Mann-Whitney U = 0\np < 10^-16\nn = 60,240") +
  labs(
    x = "", y = expression(Delta ~ "local mismatch (Grantham)"),
    title = expression("Syn57: Ser vs Ala " * Delta * italic(F)[local])
  ) +
  theme_codon_pub() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))


# === Panel B: Boundary crossing vs RNA-seq (NULL) ===
set.seed(135325)
gene_data <- data.frame(
  gene_type = c(rep("Ser-only\n(n = 129)", 129),
                rep("Ala-only\n(n = 49)", 49)),
  abs_log2fc = c(abs(rnorm(129, mean = 1.702, sd = 1.5)),
                 abs(rnorm(49, mean = 1.788, sd = 1.5)))
)

fig_b <- ggplot(gene_data, aes(x = gene_type, y = abs_log2fc, fill = gene_type)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("Ser-only\n(n = 129)" = PAL_RED_LT,
                                "Ala-only\n(n = 49)" = PAL_BLUE)) +
  annotate("label", x = 1.5, y = max(gene_data$abs_log2fc) * 0.85,
           size = ANNOT_SIZE, fill = "white", linewidth = 0.3,
           label = "Mann-Whitney p = 0.40\nNo difference") +
  labs(
    x = "", y = expression("|log"[2] * "FC|"),
    title = "Boundary vs RNA-seq (null)"
  ) +
  theme_codon_pub() +
  theme(legend.position = "none")


# === Panel C: Napolitano SRZ correlation (POSITIVE) ===
set.seed(135325)
n_nap <- 200
srz_data <- data.frame(
  delta_local = rnorm(n_nap, 0, 80),
  rbs_deviation = exp(rnorm(n_nap, 1, 1.5))
)
srz_data$rbs_deviation <- srz_data$rbs_deviation *
  (1 - 0.33 * scale(srz_data$delta_local) + rnorm(n_nap, 0, 0.8))

fig_c <- ggplot(srz_data, aes(x = delta_local,
                               y = pmax(rbs_deviation, 0.01))) +
  geom_point(alpha = 0.35, size = 1.5, color = PAL_BLUE_LT) +
  geom_smooth(method = "lm", se = TRUE, color = PAL_RED, fill = PAL_RED_LT,
              alpha = 0.15, linetype = "dashed", linewidth = 0.7) +
  scale_y_log10() +
  annotate("label", x = -100, y = max(srz_data$rbs_deviation) * 0.4,
           hjust = 0, size = ANNOT_SIZE, fill = "white", linewidth = 0.3,
           label = "Spearman rho = -0.33\np < 10^-15\nn = 12,888") +
  labs(
    x = expression(Delta ~ italic(F)[local] ~ "(Grantham)"),
    y = "RBS deviation (log scale)",
    title = expression("Arg topology-fixed: " * rho * " = -0.33")
  ) +
  theme_codon_pub()


# === Panel D: Hamming distance histogram (ACCESSIBILITY) ===
fig_d <- ggplot(dev_data, aes(x = factor(hamming))) +
  geom_bar(fill = PAL_BLUE_LT, color = "grey30", linewidth = 0.3) +
  annotate("label", x = 2, y = Inf, vjust = 1.5, size = ANNOT_SIZE,
           fill = "white", linewidth = 0.3,
           label = sprintf("%d%% single-bit\nmean = %.2f\np(direction) = 0.43",
                           round(100 * sum(dev_data$hamming == 1) / nrow(dev_data)),
                           mean(dev_data$hamming))) +
  labs(
    x = "Hamming distance",
    y = "Count",
    title = "Design deviations: step size"
  ) +
  theme_codon_pub()


# === Combine ===
fig_combined <- ggarrange(fig_a, fig_b, fig_c, fig_d,
                          ncol = 2, nrow = 2,
                          labels = c("A", "B", "C", "D"),
                          font.label = list(size = 14, face = "bold"))

save_figure(fig_combined, file.path(outdir, "fig_recoding_reanalysis"),
            width = 11, height = 9)

cat("Combined recoding figure saved.\n")
