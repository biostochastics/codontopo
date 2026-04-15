#!/usr/bin/env Rscript
# ============================================================
# CODON-TOPO: Supplementary figures (A-G + combined panel)
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
  library(jsonlite)
  library(scales)
})

# Source shared theme
theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
if (!file.exists(theme_path)) {
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  theme_path <- file.path(script_dir, "theme_codon.R")
}
source(theme_path)

args <- commandArgs(trailingOnly = TRUE)
table_dir  <- if (length(args) >= 1) args[1] else "output/tables"
output_dir <- if (length(args) >= 2) args[2] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating supplementary figures...\n\n")


# ═══════════════════════════════════════════════════════════════
# Fig A: Coloring optimality null distribution
# ═══════════════════════════════════════════════════════════════
cat("  [A] Coloring optimality null distribution\n")
mc <- fromJSON(file.path(table_dir, "T3_coloring_optimality.json"))
obs   <- mc$observed_score
mu    <- mc$null_mean
sigma <- mc$null_std

x_range <- seq(mu - 4 * sigma, mu + 4 * sigma, length.out = 400)
y_dens  <- dnorm(x_range, mu, sigma)
df_null <- data.frame(score = x_range, density = y_dens)

# Shade area below observed
df_shade <- df_null %>% filter(score <= obs)

pA <- ggplot(df_null, aes(x = score, y = density)) +
  geom_area(fill = PAL_BLUE_PALE, alpha = 0.6) +
  geom_area(data = df_shade, fill = PAL_BLUE, alpha = 0.25) +
  geom_line(color = PAL_BLUE, linewidth = 0.9) +
  geom_vline(xintercept = obs, color = PAL_RED, linewidth = 1.1) +
  annotate("label", x = obs - 80, y = max(y_dens) * 0.88,
           label = paste0("Observed\n", format(round(obs), big.mark = ",")),
           color = PAL_RED, fill = "white", linewidth = 0.3,
           size = ANNOT_SIZE, fontface = "bold", hjust = 1) +
  annotate("label", x = mu + 0.6 * sigma, y = max(y_dens) * 0.65,
           label = paste0("Null: ", format(round(mu), big.mark = ","),
                          " \u00b1 ", round(sigma),
                          "\np = ", format(mc$p_value_conservative, digits = 2)),
           color = PAL_BLUE, fill = "white", linewidth = 0.3,
           size = ANNOT_SIZE, hjust = 0) +
  labs(
    title = "Hypercube Coloring Optimality",
    subtitle = paste0("Standard code vs ", format(mc$n_samples, big.mark = ","),
                      " Freeland-Hurst block-preserving random colorings"),
    x = expression("Grantham edge-mismatch score " * italic(F) * "(code)"),
    y = "Density (normal approximation)"
  ) +
  theme_codon_pub()

save_figure(pA, file.path(output_dir, "FigA_coloring_null"))


# ═══════════════════════════════════════════════════════════════
# Fig B: Per-table optimality
# ═══════════════════════════════════════════════════════════════
cat("  [B] Per-table optimality\n")
pt <- read.csv(file.path(table_dir, "T4_per_table_optimality.csv"))
pt$significant <- ifelse(pt$p_value < 0.05, "p < 0.05", "n.s.")

pB <- ggplot(pt, aes(x = reorder(factor(table_id), quantile_pct),
                      y = quantile_pct, fill = significant)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = 5, linetype = "dashed", color = PAL_RED, linewidth = 0.5) +
  annotate("text", x = 2, y = 5.8, label = "5% threshold",
           color = PAL_RED, hjust = 0, size = ANNOT_SIZE) +
  scale_fill_manual(values = c("p < 0.05" = PAL_BLUE, "n.s." = PAL_GREY_LT)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Coloring Optimality Preserved Across Variant Codes",
    subtitle = expression("Each NCBI table tested against its own block-preserving null (" *
                          italic(n) * " = 1,000)"),
    x = "NCBI translation table",
    y = "Quantile of observed score (%)",
    fill = NULL
  ) +
  coord_flip() +
  theme_codon_pub() +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3))

save_figure(pB, file.path(output_dir, "FigB_per_table_optimality"), width = 7, height = 6.5)


# ═══════════════════════════════════════════════════════════════
# Fig C: Rho robustness sweep
# ═══════════════════════════════════════════════════════════════
cat("  [C] Rho robustness sweep\n")
rho <- read.csv(file.path(table_dir, "T5_rho_robustness.csv"))

# Shade region below p = 0.01 for emphasis
pC <- ggplot(rho, aes(x = rho, y = p_value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.01,
           fill = PAL_TEAL, alpha = 0.08) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = PAL_RED, linewidth = 0.5) +
  geom_hline(yintercept = 0.01, linetype = "dotted", color = PAL_RED, linewidth = 0.5) +
  geom_line(color = PAL_BLUE, linewidth = 1) +
  geom_point(color = PAL_BLUE, size = 2.5, shape = 16) +
  annotate("text", x = 0.88, y = 0.054, label = expression(italic(p) ~ "= 0.05"),
           color = PAL_RED, size = ANNOT_SIZE) +
  annotate("text", x = 0.88, y = 0.014, label = expression(italic(p) ~ "= 0.01"),
           color = PAL_RED, size = ANNOT_SIZE) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, max(0.065, max(rho$p_value) * 1.2)),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = expression("Optimality Robust Across Transversion Weight " * rho),
    subtitle = expression(italic(F)[rho] * " = " * italic(F)[H1] * " + " * rho *
               " \u00b7 " * italic(F)[diag] *
               ";   " * rho * " = 0: " * Q[6] * " only;   " *
               rho * " = 1: full " * K[4]^3),
    x = expression("Transversion weight " * rho),
    y = expression(italic(p) * "-value (conservative)")
  ) +
  theme_codon_pub()

save_figure(pC, file.path(output_dir, "FigC_rho_robustness"))


# ═══════════════════════════════════════════════════════════════
# Fig D: Score decomposition by nucleotide position
# ═══════════════════════════════════════════════════════════════
cat("  [D] Score decomposition\n")
dec <- fromJSON(file.path(table_dir, "T6_score_decomposition.json"))
df_dec <- data.frame(
  position = c("Position 1\n(bits 0-1)", "Position 2\n(bits 2-3)",
                "Position 3\n(wobble, bits 4-5)"),
  score = c(dec$by_nucleotide_position$pos1,
            dec$by_nucleotide_position$pos2,
            dec$by_nucleotide_position$pos3_wobble),
  fraction = c(dec$position_fractions$pos1,
               dec$position_fractions$pos2,
               dec$position_fractions$pos3_wobble)
)
df_dec$position <- factor(df_dec$position,
                          levels = df_dec$position[order(df_dec$score, decreasing = TRUE)])

pD <- ggplot(df_dec, aes(x = position, y = score, fill = position)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(fraction * 100, 1), "%")),
            vjust = -0.5, size = 4.2, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c(PAL_PURPLE, PAL_BLUE, PAL_TEAL)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Mismatch Score Decomposition by Codon Position",
    subtitle = paste0("Total score: ", format(round(dec$total_score), big.mark = ","),
                      " -- position 2 dominates (amino acid identity changes)"),
    x = NULL,
    y = "Grantham mismatch contribution"
  ) +
  theme_codon_pub() +
  theme(legend.position = "none")

save_figure(pD, file.path(output_dir, "FigD_score_decomposition"), width = 6, height = 5)


# ═══════════════════════════════════════════════════════════════
# Fig E: tRNA enrichment — AA rank per pairing
# ═══════════════════════════════════════════════════════════════
cat("  [E] tRNA enrichment AA rank\n")
trna <- read.csv(file.path(table_dir, "T7_trna_per_pairing.csv"))
trna$label <- paste0(trna$reassigned_aa,
                     " (", gsub(" vs .*", "", trna$pairing), ")")

# Clamp exact_p for color scale
trna$exact_p_clamped <- pmin(trna$exact_p, 0.3)

pE <- ggplot(trna, aes(x = reorder(label, aa_rank), y = aa_rank,
                         fill = exact_p_clamped)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = paste0("p = ", sprintf("%.3f", exact_p))),
            hjust = -0.1, size = ANNOT_SIZE_S, color = "grey35") +
  geom_hline(yintercept = 1, color = "grey80", linewidth = 0.3) +
  scale_fill_gradient(low = PAL_BLUE, high = PAL_GOLD,
                      limits = c(0, 0.3), name = expression(italic(p)[exact])) +
  scale_y_reverse(breaks = 1:20,
                  expand = expansion(mult = c(0.15, 0.02))) +
  labs(
    title = "Rank of Reassigned AA Among All tRNA Gene Counts",
    subtitle = "Rank 1 = most enriched amino acid in disconnection vs control comparison",
    x = NULL,
    y = "Rank among 20 amino acids (1 = most enriched)"
  ) +
  coord_flip() +
  theme_codon_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3))

save_figure(pE, file.path(output_dir, "FigE_trna_aa_rank"), width = 7.5, height = 5.5)


# ═══════════════════════════════════════════════════════════════
# Fig F: Topology avoidance
# ═══════════════════════════════════════════════════════════════
cat("  [F] Topology avoidance\n")
ta <- read.csv(file.path(table_dir, "T9_topology_avoidance.csv"))
df_ta <- data.frame(
  category = factor(c("Observed\n(natural)", "Possible\n(all variants)"),
                    levels = c("Observed\n(natural)", "Possible\n(all variants)")),
  rate = c(ta$rate_observed_pct[1], ta$rate_possible_pct[1]),
  group = c("Observed", "Possible")
)

pF <- ggplot(df_ta, aes(x = category, y = rate, fill = group)) +
  geom_col(width = 0.5, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(rate, 1), "%")),
            vjust = -0.5, size = 5, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("Observed" = PAL_BLUE, "Possible" = PAL_GREY_LT)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "Natural Reassignments Avoid\nTopology-Breaking Changes",
    subtitle = paste0("Permutation p = ",
                      format(ta$permutation_p[1], digits = 2, scientific = TRUE),
                      " (n = 10,000)  |  Hypergeometric p = ",
                      format(ta$hypergeom_p[1], digits = 2, scientific = TRUE)),
    x = NULL,
    y = "Fraction creating new disconnections (%)"
  ) +
  theme_codon_pub() +
  theme(legend.position = "none")

save_figure(pF, file.path(output_dir, "FigF_topology_avoidance"), width = 5.5, height = 5)


# ═══════════════════════════════════════════════════════════════
# Fig G: Bit-position bias (6-bin observed vs expected)
# ═══════════════════════════════════════════════════════════════
cat("  [G] Bit-position bias\n")
bb <- read.csv(file.path(table_dir, "T8_bit_bias.csv"))
obs_str  <- as.character(bb$observed[bb$test == "6-bin uniform"])
obs_vals <- as.numeric(strsplit(gsub("\\[|\\]| ", "", obs_str), ",")[[1]])

df_bb <- data.frame(
  bit = paste0("Bit ", 0:5),
  position = factor(c("Pos 1", "Pos 1", "Pos 2", "Pos 2",
                       "Pos 3 (wobble)", "Pos 3 (wobble)"),
                    levels = c("Pos 1", "Pos 2", "Pos 3 (wobble)")),
  observed = obs_vals,
  expected = rep(sum(obs_vals) / 6, 6)
)
df_bb$bit <- factor(df_bb$bit, levels = paste0("Bit ", 0:5))

pG <- ggplot(df_bb, aes(x = bit, y = observed, fill = position)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = df_bb$expected[1], linetype = "dashed",
             color = PAL_RED, linewidth = 0.5) +
  geom_text(aes(label = observed), vjust = -0.4, size = 3.8,
            fontface = "bold", color = "grey20") +
  annotate("text", x = 6, y = df_bb$expected[1] + 0.6,
           label = paste0("Expected: ", round(df_bb$expected[1], 1)),
           color = PAL_RED, size = ANNOT_SIZE, hjust = 1) +
  scale_fill_manual(values = c("Pos 1" = PAL_PURPLE, "Pos 2" = PAL_BLUE,
                                "Pos 3 (wobble)" = PAL_TEAL),
                    name = "Nucleotide\nposition") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Bit-Position Bias in Codon Reassignment Events",
    subtitle = paste0("n = ", sum(obs_vals), " bit-flips across ",
                      bb$n_events[bb$test == "6-bin uniform"],
                      " events  |  Bit 4 (wobble) dominates"),
    x = expression("GF(2)"^6 ~ "bit position"),
    y = "Number of bit-flips"
  ) +
  theme_codon_pub()

save_figure(pG, file.path(output_dir, "FigG_bit_position_bias"))


# ═══════════════════════════════════════════════════════════════
# Combined panel figure
# ═══════════════════════════════════════════════════════════════
cat("  [Panel] Combined supplementary figure\n")

panel <- (pA | pC) / (pB | pF) / (pE | pD) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E", "F")),
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

save_figure(panel, file.path(output_dir, "panel_strengthened"), width = 14, height = 16)

cat("\nAll supplementary figures saved to:", output_dir, "\n")
