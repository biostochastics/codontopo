#!/usr/bin/env Rscript
# CODON-TOPO: Strengthened analyses — publication figures
# ggplot2 + ggpubr, 300 DPI, colorblind-friendly (viridis)
#
# Figures:
#   Fig A: Coloring optimality null distribution with observed score
#   Fig B: Per-table optimality across all 25 NCBI tables
#   Fig C: Rho robustness sweep
#   Fig D: Score decomposition by nucleotide position
#   Fig E: tRNA enrichment — AA rank per pairing
#   Fig F: Topology avoidance — observed vs possible rates
#   Fig G: Bit-position bias — 6-bin and 3-bin histograms

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(patchwork)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
table_dir  <- if (length(args) >= 1) args[1] else "output/tables"
output_dir <- if (length(args) >= 2) args[2] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

theme_codon <- theme_pubr(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right"
  )

DPI <- 300
W <- 7
H <- 5

cat("Generating strengthened analysis figures...\n")


# ═══════════════════════════════════════════════════════════════════
# Fig A: Coloring optimality null distribution
# ═══════════════════════════════════════════════════════════════════
cat("  [A] Coloring optimality null distribution...\n")
mc <- fromJSON(file.path(table_dir, "T3_coloring_optimality.json"))
# mc$null_scores is not exported; use summary stats to draw a normal approx
obs <- mc$observed_score
mu <- mc$null_mean
sigma <- mc$null_std
x_range <- seq(mu - 4*sigma, mu + 4*sigma, length.out = 300)
y_dens <- dnorm(x_range, mu, sigma)
df_null <- data.frame(score = x_range, density = y_dens)

pA <- ggplot(df_null, aes(x = score, y = density)) +
  geom_area(fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = obs, color = "firebrick", linewidth = 1.2, linetype = "solid") +
  annotate("text", x = obs - 100, y = max(y_dens) * 0.9,
           label = paste0("Observed\n", round(obs, 0)),
           color = "firebrick", hjust = 1, size = 3.5, fontface = "bold") +
  annotate("text", x = mu + 0.5 * sigma, y = max(y_dens) * 0.7,
           label = paste0("Null: ", round(mu, 0), " +/- ", round(sigma, 0), "\n",
                          "p = ", format(mc$p_value_conservative, digits = 3)),
           color = "steelblue", hjust = 0, size = 3.5) +
  labs(
    title = "Hypercube Coloring Optimality",
    subtitle = paste0("Standard code vs ", mc$n_samples,
                      " Freeland-Hurst block-preserving random colorings"),
    x = "Grantham edge-mismatch score F(code)",
    y = "Density (normal approximation)"
  ) +
  theme_codon

ggsave(file.path(output_dir, "FigA_coloring_null.png"), pA,
       width = W, height = H, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Fig B: Per-table optimality
# ═══════════════════════════════════════════════════════════════════
cat("  [B] Per-table optimality...\n")
pt <- read.csv(file.path(table_dir, "T4_per_table_optimality.csv"))
pt$significant <- ifelse(pt$p_value < 0.05, "p < 0.05", "n.s.")

pB <- ggplot(pt, aes(x = reorder(factor(table_id), quantile_pct),
                      y = quantile_pct, fill = significant)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "firebrick", linewidth = 0.6) +
  annotate("text", x = 1, y = 6, label = "5% threshold",
           color = "firebrick", hjust = 0, size = 3) +
  scale_fill_manual(values = c("p < 0.05" = "steelblue", "n.s." = "grey60")) +
  labs(
    title = "Coloring Optimality Preserved Across Variant Codes",
    subtitle = "Each NCBI table tested against its own block-preserving null (n=1000)",
    x = "NCBI Translation Table",
    y = "Quantile of observed score (%)",
    fill = ""
  ) +
  coord_flip() +
  theme_codon +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "FigB_per_table_optimality.png"), pB,
       width = W, height = 6.5, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Fig C: Rho robustness sweep
# ═══════════════════════════════════════════════════════════════════
cat("  [C] Rho robustness sweep...\n")
rho <- read.csv(file.path(table_dir, "T5_rho_robustness.csv"))

pC <- ggplot(rho, aes(x = rho, y = p_value)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "firebrick") +
  geom_hline(yintercept = 0.01, linetype = "dotted", color = "firebrick") +
  annotate("text", x = 0.85, y = 0.055, label = "p = 0.05",
           color = "firebrick", size = 3) +
  annotate("text", x = 0.85, y = 0.015, label = "p = 0.01",
           color = "firebrick", size = 3) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, max(0.06, max(rho$p_value) * 1.2))) +
  labs(
    title = expression(paste("Optimality Robust Across Transversion Weight ", rho)),
    subtitle = expression(paste("F"[rho], " = F"["H1"], " + ", rho, " * F"["diag"],
                                ";  ", rho, "=0: Q"[6], " only;  ", rho, "=1: full K"[4]^3)),
    x = expression(paste("Transversion weight ", rho)),
    y = "P-value (conservative)"
  ) +
  theme_codon

ggsave(file.path(output_dir, "FigC_rho_robustness.png"), pC,
       width = W, height = H, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Fig D: Score decomposition by nucleotide position
# ═══════════════════════════════════════════════════════════════════
cat("  [D] Score decomposition...\n")
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
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(fraction * 100, 1), "%")),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_viridis_d(option = "D", end = 0.8) +
  labs(
    title = "Mismatch Score by Codon Position",
    subtitle = paste0("Total score: ", round(dec$total_score, 0),
                      " — Position 2 dominates (AA identity changes)"),
    x = "",
    y = "Grantham mismatch contribution"
  ) +
  theme_codon +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "FigD_score_decomposition.png"), pD,
       width = 6, height = H, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Fig E: tRNA enrichment — AA rank per pairing
# ═══════════════════════════════════════════════════════════════════
cat("  [E] tRNA enrichment AA rank...\n")
trna <- read.csv(file.path(table_dir, "T7_trna_per_pairing.csv"))
trna$label <- paste0(trna$reassigned_aa, " (", sub(" vs .*", "", trna$pairing), ")")

pE <- ggplot(trna, aes(x = reorder(label, aa_rank), y = aa_rank, fill = exact_p)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey80") +
  scale_fill_viridis_c(option = "C", direction = -1, limits = c(0, 0.3),
                        name = "Exact p") +
  scale_y_reverse(breaks = 1:20) +
  labs(
    title = "Rank of Reassigned AA Among All tRNA Gene Counts",
    subtitle = "Rank 1 = most enriched AA in disconnection vs control comparison",
    x = "",
    y = "Rank among 20 amino acids (1 = most enriched)"
  ) +
  coord_flip() +
  theme_codon

ggsave(file.path(output_dir, "FigE_trna_aa_rank.png"), pE,
       width = W, height = 5.5, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Fig F: Topology avoidance
# ═══════════════════════════════════════════════════════════════════
cat("  [F] Topology avoidance...\n")
ta <- read.csv(file.path(table_dir, "T9_topology_avoidance.csv"))
df_ta <- data.frame(
  category = c("Observed\n(natural)", "Possible\n(all variants)"),
  rate = c(ta$rate_observed_pct[1], ta$rate_possible_pct[1]),
  group = c("Observed", "Possible")
)

pF <- ggplot(df_ta, aes(x = category, y = rate, fill = group)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(round(rate, 1), "%")),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Observed" = "steelblue", "Possible" = "grey60")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Natural Reassignments Avoid Topology-Breaking Changes",
    subtitle = paste0("Permutation p = ",
                      format(ta$permutation_p[1], digits = 2, scientific = TRUE),
                      " (n=10,000) | Hypergeometric p = ",
                      format(ta$hypergeom_p[1], digits = 2, scientific = TRUE)),
    x = "",
    y = "Fraction creating new disconnections (%)"
  ) +
  theme_codon +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "FigF_topology_avoidance.png"), pF,
       width = 5, height = H, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Fig G: Bit-position bias (6-bin observed vs expected)
# ═══════════════════════════════════════════════════════════════════
cat("  [G] Bit-position bias...\n")
# Parse from the bit_bias CSV
bb <- read.csv(file.path(table_dir, "T8_bit_bias.csv"))
# Manually build the 6-bin histogram from the "6-bin de-duplicated" row
# Using the full dataset for the figure
obs_str <- as.character(bb$observed[bb$test == "6-bin uniform"])
obs_vals <- as.numeric(strsplit(gsub("\\[|\\]| ", "", obs_str), ",")[[1]])

df_bb <- data.frame(
  bit = paste0("Bit ", 0:5),
  position = c("Pos 1", "Pos 1", "Pos 2", "Pos 2", "Pos 3\n(wobble)", "Pos 3\n(wobble)"),
  observed = obs_vals,
  expected = rep(sum(obs_vals) / 6, 6)
)
df_bb$bit <- factor(df_bb$bit, levels = paste0("Bit ", 0:5))

pG <- ggplot(df_bb, aes(x = bit, y = observed, fill = position)) +
  geom_col(width = 0.7) +
  geom_hline(aes(yintercept = expected[1]), linetype = "dashed", color = "firebrick") +
  annotate("text", x = 6, y = df_bb$expected[1] + 0.8,
           label = paste0("Expected (uniform): ", round(df_bb$expected[1], 1)),
           color = "firebrick", size = 3, hjust = 1) +
  scale_fill_viridis_d(option = "D", end = 0.8) +
  labs(
    title = "Bit-Position Bias in Codon Reassignment Events",
    subtitle = paste0("n=", sum(obs_vals), " bit-flips across ",
                      bb$n_events[bb$test == "6-bin uniform"],
                      " events | Bit 4 (wobble) dominates"),
    x = "GF(2)^6 bit position",
    y = "Number of bit-flips"
  ) +
  theme_codon

ggsave(file.path(output_dir, "FigG_bit_position_bias.png"), pG,
       width = W, height = H, dpi = DPI)


# ═══════════════════════════════════════════════════════════════════
# Combined panel figure
# ═══════════════════════════════════════════════════════════════════
cat("  [Panel] Combined figure...\n")
panel <- (pA | pC) / (pB | pF) / (pE | pD)
ggsave(file.path(output_dir, "panel_strengthened.png"), panel,
       width = 14, height = 16, dpi = DPI)

cat("Done. Figures written to", output_dir, "\n")
