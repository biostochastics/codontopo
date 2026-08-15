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

# Plot the ACTUAL empirical histogram of the 10,000 null draws (not a
# Gaussian approximation). null_samples was added to monte_carlo_null in
# T15; the caption in supplement.typ describes this as an empirical
# histogram rather than a smooth density.
stopifnot(!is.null(mc$null_samples),
          length(mc$null_samples) == mc$n_samples)
df_null <- data.frame(score = as.numeric(mc$null_samples))

pA <- ggplot(df_null, aes(x = score)) +
  geom_histogram(bins = 60, fill = PAL_BLUE, alpha = 0.65,
                 colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = obs, color = PAL_RED, linewidth = 1.1) +
  annotate("label",
           x = obs, y = Inf,
           label = paste0("Observed\n", format(round(obs), big.mark = ",")),
           color = PAL_RED, fill = "white", linewidth = 0.3,
           size = ANNOT_SIZE, fontface = "bold",
           hjust = 1.05, vjust = 1.4) +
  annotate("label",
           x = mu + 0.6 * sigma, y = Inf,
           label = paste0("Null: ", format(round(mu), big.mark = ","),
                          " \u00b1 ", round(sigma),
                          "\np = ", format(mc$p_value_conservative, digits = 2),
                          " (", mc$n_beaten_observed, "/",
                          format(mc$n_samples, big.mark = ","), ")"),
           color = PAL_BLUE, fill = "white", linewidth = 0.3,
           size = ANNOT_SIZE, hjust = 0, vjust = 1.4) +
  labs(
    title = "Hypercube Coloring Optimality",
    subtitle = paste0("Standard code vs ", format(mc$n_samples, big.mark = ","),
                      " Freeland-Hurst block-preserving random colorings"),
    x = expression("Grantham edge-mismatch score " * italic(F) * "(code)"),
    y = "Number of null codes"
  ) +
  theme_codon_pub()

save_figure(pA, file.path(output_dir, "FigA_coloring_null"))


# ═══════════════════════════════════════════════════════════════
# Fig B: Per-table optimality
# ═══════════════════════════════════════════════════════════════
cat("  [B] Per-table optimality\n")
pt <- read.csv(file.path(table_dir, "T4_per_table_optimality.csv"))
pt$significant <- ifelse(pt$p_value < 0.05, "p < 0.05", "n.s.")

# Round-2 audit fix (T17): the per-table Monte Carlo runs with the same
# n as the headline monte_carlo block (both default to n = 10,000 in
# per_table_optimality). Read n from coloring_optimality.json so this
# subtitle can never silently drift past the actual per-table sample size.
# The previous literal "(n = 1,000)" was already stale after T3.
co_json <- fromJSON(file.path(dirname(table_dir), "coloring_optimality.json"))
n_per_table <- co_json$monte_carlo$n_samples
if (is.null(n_per_table)) n_per_table <- co_json$n_samples
if (is.null(n_per_table)) {
  # Fall back to the p-value floor: p_conservative = (rank+1)/(n+1),
  # so at rank = 0 (best) we recover n = round(1/min(p) - 1).
  n_per_table <- as.integer(round(1 / min(pt$p_value) - 1))
}

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
    subtitle = bquote(
      "Each NCBI table tested against its own block-preserving null (" *
      italic(n) == .(format(n_per_table, big.mark = ",")) * ")"
    ),
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
cat("  [C] Rho robustness sweep (z-score)\n")
rho <- read.csv(file.path(table_dir, "T5_rho_robustness.csv"))
rho$z <- (rho$null_mean - rho$observed_score) / rho$null_std

# n_MC drives the resolution floor. rho_sweep in coloring_optimality.json
# does not currently emit its own n_samples, so we fall back to
# monte_carlo$n_samples (same value under the T3 rerun) and, as a last
# resort, infer from the CSV's minimum p-value under
#   p_conservative = (rank + 1) / (n + 1)  =>  n = round(2 / min(p) - 1)
# using the fact that the smallest observed rank is rank = 1.
co_json <- fromJSON(file.path(dirname(table_dir), "coloring_optimality.json"))
n_mc <- co_json$rho_sweep$n_samples
if (is.null(n_mc)) n_mc <- co_json$monte_carlo$n_samples
if (is.null(n_mc)) n_mc <- co_json$n_samples
if (is.null(n_mc)) n_mc <- as.integer(round(2 / min(rho$p_value) - 1))

z_at <- function(target_rho) rho$z[which.min(abs(rho$rho - target_rho))]
z_min <- min(rho$z)
z_max <- max(rho$z)
y_pad <- 0.12 * (z_max - z_min)

# z = qnorm(1 - 0.01) ~= 2.33  (one-sided p = 0.01, the paper's threshold)
Z_P01 <- qnorm(1 - 0.01)

# Round-2 revision (T14): plot effect-size z (not p) so the monotonic
# strengthening across rho is visible. Under the MC resolution floor
# (p >= 1/(n_MC + 1)), p-values pin near the floor across most of the
# rho range and appear artificially flat, hiding the claim that the
# signal strengthens toward rho = 1. Effect-size z shows the continuous
# strengthening from Q_6 (rho = 0) to H(3,4) (rho = 1). Subtitle notation
# harmonised with the supplement caption: F_(H_1), F_(H_2), H(3,4).
pC <- ggplot(rho, aes(x = rho, y = z)) +
  geom_hline(yintercept = Z_P01, linetype = "dashed",
             color = PAL_RED, linewidth = 0.5) +
  annotate("text", x = 0.02, y = Z_P01 + 0.04,
           label = expression(italic(z) == 2.33 ~ "(" * italic(p) == 0.01 * ")"),
           color = PAL_RED, size = ANNOT_SIZE, hjust = 0) +
  geom_line(colour = PAL_BLUE, linewidth = 1) +
  geom_point(colour = PAL_BLUE, size = 2.5, shape = 16) +
  # Endpoint annotations (inside plot area) — match supplement caption notation
  annotate("text", x = 0.02, y = z_at(0),
           label = "Q[6]", parse = TRUE,
           hjust = 0, vjust = 1.9, size = 4.0, fontface = "bold",
           colour = PAL_BLUE) +
  annotate("text", x = 0.98, y = z_at(1),
           label = "H(3,4)", parse = FALSE,
           hjust = 1, vjust = -1.1, size = 4.0, fontface = "bold",
           colour = PAL_BLUE) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(min(Z_P01, z_min) - y_pad, z_max + y_pad),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = expression("Optimality Robust Across Transversion Weight " * rho),
    subtitle = bquote(
      italic(F)[rho] == italic(F)[H[1]] + rho %.% italic(F)[H[2]] * ";" ~~
      rho == 0 * ":" ~ Q[6] * ";" ~~ rho == 1 * ":" ~ H(3,4) * ";" ~~
      italic(n) == .(format(n_mc, big.mark = ","))
    ),
    x = expression("Transversion weight " * rho),
    y = expression("Effect size " * italic(z) *
                   " (standardised distance below null mean)")
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
# Include both disconnection and control organism so each of the 24 pairings
# is a unique factor level. Stripping the control (old behavior) collapsed
# 6 pairings into duplicate labels, and geom_col's default position="stack"
# then silently summed their ranks — falsely inflating the tallest bar to 6
# when the true max rank is 5.
trna$label <- paste0(
  trna$reassigned_aa, " — ",
  gsub("_(nuclear|mito)", "", gsub(" vs ", " / ", trna$pairing))
)

pE <- ggplot(trna, aes(x = reorder(label, aa_rank), y = aa_rank,
                         fill = exact_p)) +
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
# Data source: primary H(3,4) cell, definition Δβ₀ > 0 (row 1 of T9)
# All displayed stats are data-driven from the CSV (no hard-coded literals).
# ═══════════════════════════════════════════════════════════════
cat("  [F] Topology avoidance\n")
ta <- read.csv(file.path(table_dir, "T9_topology_avoidance.csv"))
df_ta <- data.frame(
  category = factor(c("Observed\n(natural)", "Possible\n(all variants)"),
                    levels = c("Observed\n(natural)", "Possible\n(all variants)")),
  rate = c(ta$rate_observed_pct[1], ta$rate_possible_pct[1]),
  group = c("Observed", "Possible")
)

fig_f_subtitle <- sprintf(
  "Observed %.1f%% vs Possible %.1f%% (RR = %.2f)\nHypergeom. p = %.1e | Permut. p <= %s (n = %s)",
  ta$rate_observed_pct[1],
  ta$rate_possible_pct[1],
  ta$risk_ratio[1],
  ta$hypergeom_p[1],
  format(ta$permutation_p[1], digits = 1),
  format(ta$n_permutations[1], big.mark = ",")
)

pF <- ggplot(df_ta, aes(x = category, y = rate, fill = group)) +
  geom_col(width = 0.5, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(rate, 1), "%")),
            vjust = -0.5, size = 5, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("Observed" = PAL_BLUE, "Possible" = PAL_GREY_LT)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "Natural Reassignments Avoid\nTopology-Breaking Changes",
    subtitle = fig_f_subtitle,
    x = NULL,
    y = expression("Fraction with " * Delta * beta[0] * " > 0 (%)")
  ) +
  theme_codon_pub() +
  theme(legend.position = "none")

save_figure(pF, file.path(output_dir, "FigF_topology_avoidance"), width = 5.5, height = 5)


# ═══════════════════════════════════════════════════════════════
# Fig G: Bit-position bias (6-bin observed vs BOTH null models)
# ═══════════════════════════════════════════════════════════════
# NOTE: This figure treats bit-position bias as EXPLORATORY / not-supported
# (Supplement §S22.1). The uniform-null apparent bit-4 excess vanishes under
# a codon-preserving null (permutation p = 1.0 in the manuscript;
# chi-square weighted p ≈ 0.027 from the JSON). Both nulls MUST be drawn so
# the reader is not misled into treating an artefactual signal as a finding.
cat("  [G] Bit-position bias (both nulls)\n")
bit_bias_json <- fromJSON(file.path(dirname(table_dir),
                                    "reassignment_analysis.json"))$bit_bias_mitochondrial
obs_vals    <- bit_bias_json$bit_counts_observed
codon_null  <- bit_bias_json$expected_counts_weighted
uniform_val <- sum(obs_vals) / 6
n_events    <- bit_bias_json$n_events

df_bb <- data.frame(
  bit = paste0("Bit ", 0:5),
  position = factor(c("Pos 1", "Pos 1", "Pos 2", "Pos 2",
                       "Pos 3 (wobble)", "Pos 3 (wobble)"),
                    levels = c("Pos 1", "Pos 2", "Pos 3 (wobble)")),
  observed = obs_vals
)
df_bb$bit <- factor(df_bb$bit, levels = paste0("Bit ", 0:5))

null_df <- data.frame(
  bit = factor(rep(paste0("Bit ", 0:5), 2), levels = paste0("Bit ", 0:5)),
  expected = c(rep(uniform_val, 6), codon_null),
  null_type = factor(rep(c("Uniform null", "Codon-preserving null"),
                          each = 6),
                     levels = c("Uniform null", "Codon-preserving null"))
)

pG <- ggplot(df_bb, aes(x = bit, y = observed)) +
  geom_col(aes(fill = position), width = 0.7,
           color = "grey30", linewidth = 0.3) +
  geom_line(data = null_df,
            aes(x = bit, y = expected,
                colour = null_type, linetype = null_type,
                group = null_type),
            linewidth = 0.8) +
  geom_point(data = subset(null_df, null_type == "Codon-preserving null"),
             aes(x = bit, y = expected, colour = null_type),
             size = 2.2, shape = 16, show.legend = FALSE) +
  geom_text(aes(label = observed), vjust = -0.4, size = 3.8,
            fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("Pos 1" = PAL_PURPLE, "Pos 2" = PAL_BLUE,
                                "Pos 3 (wobble)" = PAL_TEAL),
                    name = "Nucleotide\nposition") +
  scale_colour_manual(name = "Null model",
                      values = c("Uniform null" = PAL_RED,
                                 "Codon-preserving null" = "grey20")) +
  scale_linetype_manual(name = "Null model",
                        values = c("Uniform null" = "dashed",
                                   "Codon-preserving null" = "solid")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Bit-Position Bias in Codon Reassignment Events",
    subtitle = sprintf(
      "Apparent bit-4 excess (uniform-null χ² = %.1f, p = %.3f) vanishes under codon-preserving null (permutation p = 1.0)",
      bit_bias_json$chi2_statistic_uniform_reference,
      bit_bias_json$chi2_p_value_uniform_reference
    ),
    x = expression("GF(2)"^6 ~ "bit position"),
    y = "Number of bit-flips"
  ) +
  theme_codon_pub() +
  theme(plot.subtitle = element_text(size = 8))

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
