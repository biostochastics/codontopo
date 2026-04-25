#!/usr/bin/env Rscript
# ============================================================
# CODON-TOPO: Main manuscript figures (1-7)
# Clayworth & Kornilov 2026
# ggplot2 + ggpubr, 300 DPI, colorblind-friendly
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(patchwork)
  library(scales)
})

# Source shared theme
theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
if (!file.exists(theme_path)) {
  # Try relative to script location
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  theme_path <- file.path(script_dir, "theme_codon.R")
}
source(theme_path)

args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "output"
output_dir <- if (length(args) >= 2) args[2] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating CODON-TOPO main manuscript figures...\n\n")


# ═══════════════════════════════════════════════════════════════
# Figure 1: Persistent Homology Barcode Plot (WS1)
# ═══════════════════════════════════════════════════════════════
cat("  [1/7] Persistent homology barcode plot\n")
ph <- read.csv(file.path(input_dir, "persistent_homology.csv"))

# Focus on AAs that disconnect, plus connected AAs for contrast
disconnected <- ph %>%
  filter(epsilon == 1, beta_0 > 1) %>%
  pull(aa) %>%
  unique()

show_aas <- unique(c(disconnected, "Leu", "Arg", "Val", "Pro"))
ph_sub <- ph %>% filter(aa %in% show_aas)

# Highlight Serine vs others
ph_sub$highlight <- ifelse(ph_sub$aa == "Ser", "Serine", "Other")

p1 <- ggplot(ph_sub, aes(x = epsilon, y = beta_0, color = aa, group = aa)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(aes(linewidth = highlight, alpha = highlight)) +
  geom_point(aes(size = highlight)) +
  scale_linewidth_manual(values = c("Serine" = 1.4, "Other" = 0.7), guide = "none") +
  scale_alpha_manual(values = c("Serine" = 1.0, "Other" = 0.6), guide = "none") +
  scale_size_manual(values = c("Serine" = 3, "Other" = 1.8), guide = "none") +
  scale_color_viridis_d(option = "D", end = 0.9) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(breaks = 1:4) +
  labs(
    title = "Persistent Homology of Amino Acid Codon Graphs",
    subtitle = expression("Connected components (" * beta[0] * ") vs Hamming distance threshold " * epsilon),
    x = expression("Hamming distance threshold " * epsilon),
    y = expression(beta[0] ~ "(connected components)"),
    color = "Amino acid"
  ) +
  annotate("label", x = 2.5, y = 2.5,
           label = "Serine: disconnected at\nepsilon = 1, reconnects at epsilon = 4",
           size = ANNOT_SIZE, fill = "grey97", linewidth = 0.3,
           color = "grey25", label.padding = unit(0.3, "lines")) +
  annotate("segment", x = 2.2, xend = 1.3, y = 2.35, yend = 2.05,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = "grey50", linewidth = 0.4) +
  theme_codon_pub() +
  theme(legend.position = "right")

save_figure(p1, file.path(output_dir, "fig1_persistent_homology"), width = 7, height = 5)


# ═══════════════════════════════════════════════════════════════
# Figure 2: Disconnection Catalogue Heatmap (WS1)
# ═══════════════════════════════════════════════════════════════
cat("  [2/7] Disconnection catalogue heatmap\n")
disc <- read.csv(file.path(input_dir, "disconnection_catalogue.csv"))

# Precompute text color for heatmap readability
disc$text_color <- ifelse(disc$reconnect_eps >= 3.5, "white", "grey15")

p2 <- ggplot(disc, aes(x = reorder(aa, -reconnect_eps),
                        y = reorder(table_name, table_id),
                        fill = reconnect_eps)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = reconnect_eps, color = text_color),
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_viridis_c(option = "D", direction = -1,
                        name = expression(epsilon[reconnect])) +
  labs(
    title = "Disconnection Catalogue Across NCBI Translation Tables",
    subtitle = expression("Each tile: amino acid disconnected at " * epsilon *
                          " = 1; number = reconnection " * epsilon),
    x = "Amino acid",
    y = "Translation table"
  ) +
  theme_codon_pub() +
  theme(
    axis.text.y = element_text(size = 8.5),
    panel.grid.major.y = element_blank()
  )

save_figure(p2, file.path(output_dir, "fig2_disconnection_catalogue"), width = 7, height = 7)


# ═══════════════════════════════════════════════════════════════
# Figure 3: Bit-Position Bias (WS2)
# ═══════════════════════════════════════════════════════════════
cat("  [3/7] Bit-position bias bar chart\n")

bit_data <- data.frame(
  bit_position = factor(0:5),
  label = c("Bit 0", "Bit 1", "Bit 2", "Bit 3", "Bit 4", "Bit 5"),
  region = factor(c("1st base", "1st base", "2nd base", "2nd base",
                     "3rd base (wobble)", "3rd base (wobble)"),
                  levels = c("1st base", "2nd base", "3rd base (wobble)")),
  count = c(5, 4, 0, 8, 13, 5)
)
total_n <- sum(bit_data$count)
expected <- total_n / 6

p3 <- ggplot(bit_data, aes(x = label, y = count, fill = region)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = expected, linetype = "dashed", color = PAL_RED,
             linewidth = 0.6) +
  geom_text(aes(label = count), vjust = -0.5, fontface = "bold", size = 4) +
  annotate("text", x = 6, y = expected + 0.6,
           label = paste0("Expected: ", round(expected, 1)),
           color = PAL_RED, size = ANNOT_SIZE, hjust = 1) +
  scale_fill_manual(values = c("1st base" = PAL_PURPLE,
                                "2nd base" = PAL_BLUE,
                                "3rd base (wobble)" = PAL_TEAL),
                    name = "Nucleotide\nposition") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Bit-Position Bias in Codon Reassignment Events",
    subtitle = bquote(chi^2 ~ "= 16.26," ~ italic(p) ~ "= 0.006 (" *
               italic(n) ~ "=" ~ .(total_n) ~ "bit-changes, 27 NCBI tables)"),
    x = expression("Bit position in GF(2)"^6 ~ "vector"),
    y = "Number of bit changes"
  ) +
  theme_codon_pub()

save_figure(p3, file.path(output_dir, "fig3_bit_position_bias"))


# ═══════════════════════════════════════════════════════════════
# Figure 4: Evolutionary Depth Calibration (WS3)
# ═══════════════════════════════════════════════════════════════
cat("  [4/7] Depth calibration scatter\n")
depth <- read.csv(file.path(input_dir, "depth_calibration.csv"))

p4 <- ggplot(depth, aes(x = age_midpoint_mya, y = reconnect_eps,
                          color = aa, shape = aa)) +
  geom_errorbar(aes(xmin = age_mya_low, xmax = age_mya_high),
                width = 0.12, linewidth = 0.7, alpha = 0.7,
                orientation = "y") +
  geom_point(size = 4, stroke = 0.8) +
  scale_color_manual(values = c("Ala" = PAL_PURPLE, "Leu" = PAL_BLUE,
                                 "Ser" = PAL_TEAL, "Thr" = PAL_ORANGE)) +
  scale_shape_manual(values = c("Ala" = 16, "Leu" = 17, "Ser" = 15, "Thr" = 3)) +
  scale_y_continuous(breaks = 1:6) +
  scale_x_log10(labels = comma) +
  labs(
    title = "Evolutionary Depth Calibration",
    subtitle = expression("Spearman " * rho * " = 0.0 -- ordinal correspondence, not monotonic rank"),
    x = "Estimated divergence age (Mya, log scale)",
    y = expression("Reconnection " * epsilon),
    color = "Amino acid", shape = "Amino acid"
  ) +
  annotate("label", x = 3750, y = 4.4, label = "Serine\n(pre-LUCA)",
           size = ANNOT_SIZE, fill = "grey97", linewidth = 0.3,
           color = "grey25") +
  annotate("label", x = 135, y = 3.4, label = "CUG clade\n(~150 Mya)",
           size = ANNOT_SIZE, fill = "grey97", linewidth = 0.3,
           color = "grey25") +
  annotate("label", x = 550, y = 1.6, label = "Chlorophyceae\n(~600 Mya)",
           size = ANNOT_SIZE, fill = "grey97", linewidth = 0.3,
           color = "grey25") +
  theme_codon_pub()

save_figure(p4, file.path(output_dir, "fig4_depth_calibration"), width = 7, height = 5)


# ═══════════════════════════════════════════════════════════════
# Figure 5: KRAS Fano-Line Predictions (WS4)
# ═══════════════════════════════════════════════════════════════
cat("  [5/7] KRAS Fano predictions\n")
fano <- read.csv(file.path(input_dir, "fano_predictions.csv"))

fano$variant_label <- fano$variant
fano$partner_label <- paste0(fano$fano_partner_codon, " (", fano$fano_partner_aa, ")")
fano$mutation_label <- paste0(fano$wt_codon, " \u2192 ", fano$mutant_codon)

p5 <- ggplot(fano, aes(x = reorder(variant, variant), y = 1)) +
  geom_tile(aes(fill = fano_partner_aa), width = 0.85, height = 0.55,
            color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = partner_label), size = 3.8, fontface = "bold",
            color = "grey15") +
  geom_text(aes(label = mutation_label, y = 0.58),
            size = ANNOT_SIZE_S, color = "grey45") +
  scale_fill_manual(
    values = c("Ala" = PAL_BLUE_LT, "Arg" = PAL_BLUE, "His" = PAL_TEAL,
               "Leu" = PAL_GOLD, "Ser" = PAL_ORANGE, "Thr" = PAL_RED_LT),
    name = "Fano partner\namino acid"
  ) +
  scale_y_continuous(limits = c(0.3, 1.35), expand = c(0, 0)) +
  labs(
    title = expression("KRAS G12 Fano-Line Predictions in GF(2)"^6),
    subtitle = "For each G12X mutation (GGU \u2192 X), the Fano partner completes XOR = 0",
    x = "KRAS G12 variant"
  ) +
  theme_codon_pub() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y  = element_blank(),
    panel.grid.major.y = element_blank()
  )

save_figure(p5, file.path(output_dir, "fig5_kras_fano"), width = 7, height = 3.5)


# ═══════════════════════════════════════════════════════════════
# Figure 6: Synbio Feasibility Landscape (WS6)
# ═══════════════════════════════════════════════════════════════
cat("  [6/7] Synbio feasibility landscape\n")
synbio <- read.csv(file.path(input_dir, "synbio_landscape.csv"))

# Grouped bar chart — score has only 4 discrete values from boolean components
synbio$filtration <- ifelse(synbio$twofold_intact == "True" & synbio$fourfold_intact == "True",
                             "Both preserved", "Filtration broken")
n_total <- nrow(synbio)
n_high  <- sum(synbio$feasibility_score >= 0.8)
pct_high <- round(100 * n_high / n_total, 0)

# Aggregate counts per score x filtration
synbio_agg <- as.data.frame(table(
  score = factor(synbio$feasibility_score),
  filtration = synbio$filtration
))
names(synbio_agg)[3] <- "count"
synbio_agg$pct <- round(100 * synbio_agg$count / n_total, 1)

# Compute threshold position dynamically from factor levels
score_levels <- levels(synbio_agg$score)
cut_idx <- which(as.numeric(score_levels) >= 0.8)[1]
vline_x <- cut_idx - 0.5  # boundary between last <0.8 and first >=0.8 bin

p6 <- ggplot(synbio_agg, aes(x = score, y = count, fill = filtration)) +
  annotate("rect", xmin = vline_x, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = PAL_BLUE, alpha = 0.04) +
  geom_col(position = "dodge", width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            position = position_dodge(width = 0.7), vjust = -0.4,
            size = 3.5, fontface = "bold", color = "grey20") +
  geom_vline(xintercept = vline_x, linetype = "dashed", color = PAL_RED,
             linewidth = 0.5) +
  annotate("text", x = vline_x - 0.1, y = max(synbio_agg$count) * 0.95,
           label = paste0(pct_high, "% score\n>= 0.8"),
           hjust = 1, size = ANNOT_SIZE, color = PAL_RED) +
  scale_fill_manual(
    values = c("Both preserved" = PAL_BLUE, "Filtration broken" = PAL_ORANGE),
    name = "Degeneracy\nfiltration"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Feasibility Landscape of Single-Codon Reassignments",
    subtitle = paste0("All ", format(n_total, big.mark = ","),
                      " variants; filtration-preserving variants concentrated at highest scores"),
    x = "Feasibility score (composite: filtration + disconnection criteria)",
    y = "Number of variants"
  ) +
  theme_codon_pub()

save_figure(p6, file.path(output_dir, "fig6_synbio_landscape"))


# ═══════════════════════════════════════════════════════════════
# Figure 7: Prediction Catalogue Summary (WS5)
# ═══════════════════════════════════════════════════════════════
cat("  [7/7] Prediction catalogue summary\n")
cat_data <- read.csv(file.path(input_dir, "catalogue.csv"))

cat_data$evidence_strength <- factor(
  gsub("_", " ", tools::toTitleCase(gsub("_", " ", cat_data$evidence_strength))),
  levels = c("Very Strong", "Strong", "Moderate", "Weak", "Untested")
)
cat_data$status <- factor(cat_data$status,
  levels = c("verified", "tested", "pending", "null"))
cat_data$workstream <- factor(cat_data$workstream,
  levels = c("WS1", "WS2", "WS3", "WS4", "WS6"))

p7a <- ggplot(cat_data, aes(x = workstream, fill = status)) +
  geom_bar(width = 0.6, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = PAL_STATUS,
                    labels = c("Verified", "Tested", "Pending", "Null")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = seq(0, 10, 2)) +
  labs(title = "A. Predictions by workstream", x = NULL, y = "Count",
       fill = "Status") +
  theme_codon_pub() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

p7b <- ggplot(cat_data, aes(x = evidence_strength, fill = evidence_strength)) +
  geom_bar(width = 0.6, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Very Strong" = PAL_BLUE, "Strong" = PAL_TEAL,
    "Moderate" = PAL_ORANGE, "Weak" = PAL_GREY_LT, "Untested" = "grey90"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = seq(0, 10, 2)) +
  labs(title = "B. Evidence strength distribution", x = NULL, y = "Count") +
  theme_codon_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

p7 <- p7a + p7b +
  plot_annotation(
    title = "CODON-TOPO Prediction Catalogue",
    subtitle = "15 predictions: 9 verified, 5 tested, 1 null",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "grey35", hjust = 0),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

save_figure(p7, file.path(output_dir, "fig7_catalogue_summary"), width = 7, height = 5)


cat("\nAll main manuscript figures saved to:", output_dir, "\n")
