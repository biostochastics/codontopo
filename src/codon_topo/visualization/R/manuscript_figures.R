#!/usr/bin/env Rscript
# ============================================================
# CODON-TOPO: 5 Multi-Panel Manuscript Figures
# Clayworth & Kornilov 2026
#
# Each figure is a dense multi-panel composite packing maximum
# findings into the 5-figure manuscript budget.
#
# Fig 1: Persistent homology + Disconnection catalogue (WS1)
# Fig 2: Coloring optimality + Per-table + Rho robustness (WS1)
# Fig 3: Bit-position + Depth + Topology avoidance + tRNA (WS2/3)
# Fig 4: Synbio feasibility + Decomposition + Catalogue (WS5/6)
# Fig 5: Conditional logit diagnostics (Theory extension)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(patchwork)
  library(scales)
  library(jsonlite)
})

# Source shared theme
theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
if (!file.exists(theme_path)) {
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  theme_path <- file.path(script_dir, "theme_codon.R")
}
source(theme_path)

args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "output"
table_dir  <- if (length(args) >= 2) args[2] else "output/tables"
output_dir <- if (length(args) >= 3) args[3] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Compact theme for multi-panel (smaller base_size)
tc <- function() theme_codon_pub(base_size = 9)

cat("Generating 5 multi-panel manuscript figures...\n\n")


# ═══════════════════════════════════════════════════════════════
# FIGURE 1: Core Topology (A: Persistent Homology, B: Heatmap)
# ═══════════════════════════════════════════════════════════════
cat("  [1/5] Core topology (2-panel)\n")

# Panel A: Persistent homology
ph <- read.csv(file.path(input_dir, "persistent_homology.csv"))
disconnected <- ph %>% filter(epsilon == 1, beta_0 > 1) %>% pull(aa) %>% unique()
show_aas <- unique(c(disconnected, "Leu", "Arg", "Val", "Pro"))
ph_sub <- ph %>% filter(aa %in% show_aas)
ph_sub$highlight <- ifelse(ph_sub$aa == "Ser", "Serine", "Other")

p1a <- ggplot(ph_sub, aes(x = epsilon, y = beta_0, color = aa, group = aa)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey75", linewidth = 0.3) +
  geom_line(aes(linewidth = highlight, alpha = highlight)) +
  geom_point(aes(size = highlight)) +
  scale_linewidth_manual(values = c("Serine" = 1.2, "Other" = 0.6), guide = "none") +
  scale_alpha_manual(values = c("Serine" = 1.0, "Other" = 0.5), guide = "none") +
  scale_size_manual(values = c("Serine" = 2.5, "Other" = 1.5), guide = "none") +
  scale_color_viridis_d(option = "D", end = 0.9) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(breaks = 1:4) +
  labs(
    x = expression("Hamming threshold " * epsilon),
    y = expression(beta[0] ~ "(components)"),
    color = "AA"
  ) +
  tc() +
  theme(legend.key.size = unit(0.3, "cm"))

# Panel B: Disconnection catalogue
disc <- read.csv(file.path(input_dir, "disconnection_catalogue.csv"))
disc$text_color <- ifelse(disc$reconnect_eps >= 3.5, "white", "grey15")

p1b <- ggplot(disc, aes(x = reorder(aa, -reconnect_eps),
                          y = reorder(table_name, table_id),
                          fill = reconnect_eps)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = reconnect_eps, color = text_color),
            size = 2.5, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_viridis_c(option = "D", direction = -1,
                        name = expression(epsilon[reconnect])) +
  labs(x = "Amino acid", y = NULL) +
  tc() +
  theme(axis.text.y = element_text(size = 6.5),
        panel.grid.major.y = element_blank(),
        legend.key.height = unit(0.6, "cm"))

fig1 <- (p1a | p1b) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_figure(fig1, file.path(output_dir, "Fig1_core_topology"), width = 7, height = 5.5)


# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Coloring Optimality (A: Null, B: Per-table, C: Rho)
# ═══════════════════════════════════════════════════════════════
cat("  [2/5] Coloring optimality (3-panel)\n")

# Panel A: Null distribution — plots the actual empirical histogram of the
# 10,000 Freeland-Hurst null draws (was a Gaussian approximation before T15).
mc <- fromJSON(file.path(table_dir, "T3_coloring_optimality.json"))
obs <- mc$observed_score; mu <- mc$null_mean; sigma <- mc$null_std
stopifnot(!is.null(mc$null_samples),
          length(mc$null_samples) == mc$n_samples)
df_null <- data.frame(score = as.numeric(mc$null_samples))

p2a <- ggplot(df_null, aes(x = score)) +
  geom_histogram(bins = 60, fill = PAL_BLUE, alpha = 0.65,
                 colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = obs, color = PAL_RED, linewidth = 0.9) +
  annotate("label",
           x = obs, y = Inf,
           label = paste0("Obs: ", format(round(obs), big.mark = ",")),
           color = PAL_RED, fill = "white", linewidth = 0.2,
           size = ANNOT_SIZE_S, fontface = "bold",
           hjust = -0.08, vjust = 1.4) +
  annotate("label",
           x = mu + 0.5*sigma, y = Inf,
           label = paste0("p = ", format(mc$p_value_conservative, digits = 2),
                          " (", mc$n_beaten_observed, "/",
                          format(mc$n_samples, big.mark = ","), ")"),
           color = PAL_BLUE, fill = "white", linewidth = 0.2,
           size = ANNOT_SIZE_S, hjust = 0, vjust = 1.4) +
  labs(x = expression(italic(F) * "(code)"), y = "Null codes") +
  tc()

# Panel B: Per-table optimality
pt <- read.csv(file.path(table_dir, "T4_per_table_optimality.csv"))
pt$significant <- ifelse(pt$p_value < 0.05, "p < 0.05", "n.s.")

p2b <- ggplot(pt, aes(x = reorder(factor(table_id), quantile_pct),
                        y = quantile_pct, fill = significant)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.15) +
  geom_hline(yintercept = 5, linetype = "dashed", color = PAL_RED, linewidth = 0.4) +
  scale_fill_manual(values = c("p < 0.05" = PAL_BLUE, "n.s." = PAL_GREY_LT)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "NCBI table", y = "Quantile (%)", fill = NULL) +
  coord_flip() +
  tc() +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.2),
        axis.text.y = element_text(size = 6))

# Panel C: Rho robustness
rho <- read.csv(file.path(table_dir, "T5_rho_robustness.csv"))

p2c <- ggplot(rho, aes(x = rho, y = p_value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.01,
           fill = PAL_TEAL, alpha = 0.08) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = PAL_RED, linewidth = 0.4) +
  geom_hline(yintercept = 0.01, linetype = "dotted", color = PAL_RED, linewidth = 0.4) +
  geom_line(color = PAL_BLUE, linewidth = 0.8) +
  geom_point(color = PAL_BLUE, size = 1.8) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, max(0.065, max(rho$p_value) * 1.2))) +
  labs(x = expression(rho), y = expression(italic(p) * "-value")) +
  tc()

fig2 <- (p2a / p2c) | p2b
fig2 <- fig2 +
  plot_layout(widths = c(1, 0.8)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_figure(fig2, file.path(output_dir, "Fig2_coloring_optimality"), width = 7, height = 6)


# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Evolutionary Evidence (A: Bit-bias, B: Depth,
#            C: Topology avoidance, D: tRNA)
# ═══════════════════════════════════════════════════════════════
cat("  [3/5] Evolutionary evidence (4-panel)\n")

# Panel A: Bit-position bias (BOTH nulls — see Supplement §S22.1)
# Uniform-null apparent bit-4 excess vanishes under codon-preserving null
# (permutation p = 1.0). Both nulls MUST be drawn so the reader is not
# misled into treating an artefactual signal as a positive finding.
bit_bias_json <- fromJSON(file.path(input_dir,
                                    "reassignment_analysis.json"))$bit_bias_mitochondrial
obs_vals    <- bit_bias_json$bit_counts_observed
codon_null  <- bit_bias_json$expected_counts_weighted
uniform_val <- sum(obs_vals) / 6

bit_data <- data.frame(
  bit = paste0("Bit ", 0:5),
  position = factor(c("1st", "1st", "2nd", "2nd", "3rd (wob)", "3rd (wob)"),
                    levels = c("1st", "2nd", "3rd (wob)")),
  count = obs_vals
)
bit_data$bit <- factor(bit_data$bit, levels = paste0("Bit ", 0:5))

null_df <- data.frame(
  bit = factor(rep(paste0("Bit ", 0:5), 2), levels = paste0("Bit ", 0:5)),
  expected = c(rep(uniform_val, 6), codon_null),
  null_type = factor(rep(c("Uniform null", "Codon-preserving null"),
                          each = 6),
                     levels = c("Uniform null", "Codon-preserving null"))
)

p3a <- ggplot(bit_data, aes(x = bit, y = count)) +
  geom_col(aes(fill = position), width = 0.65,
           color = "grey30", linewidth = 0.2) +
  geom_line(data = null_df,
            aes(x = bit, y = expected,
                colour = null_type, linetype = null_type,
                group = null_type),
            linewidth = 0.55) +
  geom_point(data = subset(null_df, null_type == "Codon-preserving null"),
             aes(x = bit, y = expected, colour = null_type),
             size = 1.4, shape = 16, show.legend = FALSE) +
  geom_text(aes(label = count), vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("1st" = PAL_PURPLE, "2nd" = PAL_BLUE,
                                "3rd (wob)" = PAL_TEAL),
                    name = "Pos") +
  scale_colour_manual(name = "Null",
                      values = c("Uniform null" = PAL_RED,
                                 "Codon-preserving null" = "grey20")) +
  scale_linetype_manual(name = "Null",
                        values = c("Uniform null" = "dashed",
                                   "Codon-preserving null" = "solid")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = expression("GF(2)"^6 ~ "bit"), y = "Bit changes",
       subtitle = sprintf(
         "Uniform-null excess (χ²=%.1f, p=%.3f) vanishes under codon-preserving null (perm p=1.0)",
         bit_bias_json$chi2_statistic_uniform_reference,
         bit_bias_json$chi2_p_value_uniform_reference
       )) +
  tc() +
  theme(legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(size = 7),
        plot.subtitle = element_text(size = 6.5))

# Panel B: Depth calibration
depth <- read.csv(file.path(input_dir, "depth_calibration.csv"))

p3b <- ggplot(depth, aes(x = age_midpoint_mya, y = reconnect_eps,
                           color = aa, shape = aa)) +
  # Horizontal CI segments rendered manually so the right-most Ser (Pre-LUCA,
  # tight log-scale CI) is still visible above the marker glyph.
  #
  # Paired vertical jitter (round-2 fix): two pairs of points share coordinates
  # and previously overplotted — Ala tab26 vs Ser tab12 at (150 Mya, eps=3) and
  # Leu tab16 vs Leu tab22 at (600 Mya, eps=2). Applying identical
  # position_jitter(seed=DEFAULT_SEED) to every geom (segments + point) keeps
  # each marker aligned with its own CI whiskers while separating coincident
  # rows. width=0 preserves the log-scale x-position exactly.
  geom_segment(aes(x = age_mya_low, xend = age_mya_high,
                   y = reconnect_eps, yend = reconnect_eps),
               linewidth = 0.7, alpha = 0.85,
               lineend = "butt", show.legend = FALSE,
               position = position_jitter(width = 0, height = 0.15,
                                          seed = 135325)) +
  geom_segment(aes(x = age_mya_low, xend = age_mya_low,
                   y = reconnect_eps - 0.18, yend = reconnect_eps + 0.18),
               linewidth = 0.7, alpha = 0.85, show.legend = FALSE,
               position = position_jitter(width = 0, height = 0.15,
                                          seed = 135325)) +
  geom_segment(aes(x = age_mya_high, xend = age_mya_high,
                   y = reconnect_eps - 0.18, yend = reconnect_eps + 0.18),
               linewidth = 0.7, alpha = 0.85, show.legend = FALSE,
               position = position_jitter(width = 0, height = 0.15,
                                          seed = 135325)) +
  geom_point(size = 3, stroke = 0.6,
             position = position_jitter(width = 0, height = 0.15,
                                        seed = 135325)) +
  scale_color_manual(values = c("Ala" = PAL_PURPLE, "Leu" = PAL_BLUE,
                                 "Ser" = PAL_TEAL, "Thr" = PAL_ORANGE)) +
  scale_shape_manual(values = c("Ala" = 16, "Leu" = 17, "Ser" = 15, "Thr" = 3)) +
  scale_y_continuous(breaks = 1:6) +
  scale_x_log10(labels = comma) +
  labs(x = "Divergence (Mya, log)", y = expression(epsilon[reconnect]),
       color = "AA", shape = "AA") +
  tc() +
  theme(legend.key.size = unit(0.25, "cm"))

# Panel C: Topology avoidance — primary H(3,4) (encoding-independent),
# with Q_6 (encoding-dependent) shown alongside as sensitivity
stats_path <- file.path(dirname(table_dir), "manuscript_stats.json")
ms <- jsonlite::fromJSON(stats_path)
tk43 <- ms$topology_avoidance_k43
tq6  <- ms$topology_avoidance_q6

# Use plotmath-parseable factor labels so the Q_6 subscript actually renders
# in the facet strip rather than printing the literal string "Q[6]". The
# quoted "(primary)" / "(sensitivity)" segments stay upright; H(3,4) picks up
# plotmath's default italic H, matching the Typst rendering in the body text.
graph_levels <- c("H(3,4)~\"(primary)\"", "Q[6]~\"(sensitivity)\"")
df_ta <- data.frame(
  graph    = factor(rep(graph_levels, each = 2), levels = graph_levels),
  category = factor(rep(c("Observed", "Candidate"), 2),
                    levels = c("Observed", "Candidate")),
  rate     = c(tk43$rate_observed * 100, tk43$rate_possible * 100,
               tq6$rate_observed  * 100, tq6$rate_possible  * 100)
)

df_ta$category_label <- factor(
  ifelse(df_ta$category == "Observed",
         sprintf("Observed\n(n=%d)", tk43$observed_total),
         sprintf("Candidate\n(N=%s)", format(tk43$possible_total, big.mark = ","))),
  levels = c(
    sprintf("Observed\n(n=%d)", tk43$observed_total),
    sprintf("Candidate\n(N=%s)", format(tk43$possible_total, big.mark = ","))
  )
)

# Data-drive the subtitle from manuscript_stats.json (round-2 audit flagged
# the previous hard-coded "p<1e-4" as imprecise: the hypergeometric p is 1e-6
# scale). Split the p-value into mantissa/exponent so bquote() can render it
# as "1.3 x 10^-6" via plotmath rather than the literal string "1.3e-06".
p_parts   <- strsplit(formatC(tk43$hypergeom_p, format = "e", digits = 1), "e")[[1]]
p_man_str <- p_parts[1]
p_exp_int <- as.integer(p_parts[2])

subtitle_expr <- bquote(
  H(3,4) ~ "RR" == .(sprintf("%.2f", tk43$risk_ratio)) * "," ~
  italic(p) == .(p_man_str) %*% 10^.(p_exp_int) ~~~
  Q[6] ~ .(sprintf("%.1f", tq6$depletion_fold)) * "-fold depl."
)

p3c <- ggplot(df_ta, aes(x = category_label, y = rate, fill = category)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(rate, 1), "%")),
            vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("Observed" = PAL_BLUE, "Candidate" = PAL_GREY_LT)) +
  scale_y_continuous(limits = c(0, 95), expand = expansion(mult = c(0, 0.10))) +
  facet_wrap(~ graph, scales = "free_x", labeller = label_parsed) +
  labs(x = NULL, y = "% topology-breaking moves",
       subtitle = subtitle_expr) +
  tc() +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 6.0, color = PAL_RED,
                                      margin = margin(b = 4)),
        strip.text = element_text(size = 7.0,
                                   margin = margin(t = 2, b = 2)),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.spacing = unit(1.0, "lines"),
        axis.text.x = element_text(size = 6.5, lineheight = 0.85,
                                    margin = margin(t = 1)))

# Panel D: tRNA enrichment
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

p3d <- ggplot(trna, aes(x = reorder(label, aa_rank), y = aa_rank,
                          fill = exact_p)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.15) +
  geom_hline(yintercept = 1, color = "grey80", linewidth = 0.2) +
  scale_fill_gradient(low = PAL_BLUE, high = PAL_GOLD,
                      limits = c(0, 0.3), name = expression(italic(p))) +
  scale_y_reverse(breaks = c(1, 3, 5)) +
  labs(x = NULL, y = "Rank (1 = most enriched)") +
  coord_flip() +
  tc() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.2),
        axis.text.y = element_text(size = 6.5),
        legend.key.height = unit(0.4, "cm"))

fig3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_figure(fig3, file.path(output_dir, "Fig3_evolutionary_evidence"), width = 7, height = 6)


# ═══════════════════════════════════════════════════════════════
# FIGURE 4: Translational Applications (A: Synbio, B: Decomp,
#            C: Catalogue)
# ═══════════════════════════════════════════════════════════════
cat("  [4/5] Translational applications (3-panel)\n")

# Panel A: Synbio feasibility (grouped bar)
synbio <- read.csv(file.path(input_dir, "synbio_landscape.csv"))
synbio$filtration <- ifelse(synbio$twofold_intact == "True" & synbio$fourfold_intact == "True",
                             "Both preserved", "Filtration broken")
n_total <- nrow(synbio)
synbio_agg <- as.data.frame(table(
  score = factor(synbio$feasibility_score),
  filtration = synbio$filtration
))
names(synbio_agg)[3] <- "count"

score_levels <- levels(synbio_agg$score)
cut_idx <- which(as.numeric(score_levels) >= 0.8)[1]
vline_x <- cut_idx - 0.5

p4a <- ggplot(synbio_agg, aes(x = score, y = count, fill = filtration)) +
  annotate("rect", xmin = vline_x, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = PAL_BLUE, alpha = 0.04) +
  geom_col(position = "dodge", width = 0.65, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            position = position_dodge(width = 0.65), vjust = -0.3,
            size = 2.8, fontface = "bold") +
  geom_vline(xintercept = vline_x, linetype = "dashed", color = PAL_RED, linewidth = 0.4) +
  scale_fill_manual(values = c("Both preserved" = PAL_BLUE, "Filtration broken" = PAL_ORANGE),
                    name = "Filtration") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Feasibility score", y = "Variants") +
  tc() +
  theme(legend.key.size = unit(0.25, "cm"))

# Panel B: Score decomposition
dec <- fromJSON(file.path(table_dir, "T6_score_decomposition.json"))
df_dec <- data.frame(
  position = c("Pos 1", "Pos 2", "Pos 3\n(wobble)"),
  score = c(dec$by_nucleotide_position$pos1,
            dec$by_nucleotide_position$pos2,
            dec$by_nucleotide_position$pos3_wobble),
  fraction = c(dec$position_fractions$pos1,
               dec$position_fractions$pos2,
               dec$position_fractions$pos3_wobble)
)
df_dec$position <- factor(df_dec$position,
                          levels = df_dec$position[order(df_dec$score, decreasing = TRUE)])

p4b <- ggplot(df_dec, aes(x = position, y = score, fill = position)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(fraction * 100, 1), "%")),
            vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = c(PAL_PURPLE, PAL_BLUE, PAL_TEAL)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Grantham mismatch") +
  tc() +
  theme(legend.position = "none")

# Panel C: Prediction catalogue
cat_data <- read.csv(file.path(input_dir, "catalogue.csv"))
cat_data$evidence_strength <- factor(
  gsub("_", " ", tools::toTitleCase(gsub("_", " ", cat_data$evidence_strength))),
  levels = c("Very Strong", "Strong", "Moderate", "Weak", "Untested")
)
# Named vector: aligns levels -> pretty labels by NAME so that ggplot's
# empty-level drop (drop=TRUE, default) cannot re-assign labels positionally.
# Rename "null" -> "Falsified" (paper's own honest-null wording for WS4/KRAS).
# All six statuses present in output/catalogue.json are enumerated so no
# claim silently becomes NA and disappears from the bar heights.
status_labels <- c(
  verified     = "Verified",
  tested       = "Tested",
  pending      = "Pending",
  null         = "Falsified",
  tautological = "Tautological",
  exploratory  = "Exploratory"
)
cat_data$status <- factor(cat_data$status, levels = names(status_labels))
cat_data$workstream <- factor(cat_data$workstream,
  levels = c("WS1", "WS2", "WS3", "WS4", "WS6"))

# Only show legend entries for statuses actually present in the catalogue.
# The "pending" bucket, in particular, is now empty (KRAS moved to
# "Falsified") — an orphan legend entry would misrepresent the state.
present_statuses <- intersect(names(status_labels), levels(droplevels(cat_data$status)))

p4c <- ggplot(cat_data, aes(x = workstream, fill = status)) +
  geom_bar(width = 0.55, color = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = PAL_STATUS, labels = status_labels[present_statuses],
                    breaks = present_statuses, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Predictions", fill = "Status") +
  tc() +
  theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"))

fig4 <- p4a / (p4b | p4c) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_figure(fig4, file.path(output_dir, "Fig4_translational"), width = 7, height = 6)


# ═══════════════════════════════════════════════════════════════
# FIGURE 5: Conditional Logit Model (A: AICc, B: Ranks, C: LRT)
# ═══════════════════════════════════════════════════════════════
cat("  [5/5] Conditional logit diagnostics (3-panel)\n")

model_comp <- read.csv("output/tables/T_model_comparison.csv", stringsAsFactors = FALSE)
rank_files <- c("M1_phys", "M2_topo", "M3_phys_topo", "M4_full",
                "M2_topo_k43", "M3_phys_topo_k43")
all_ranks <- do.call(rbind, lapply(rank_files, function(m) {
  fp <- paste0("output/tables/T_ranks_", m, ".csv")
  if (file.exists(fp)) read.csv(fp, stringsAsFactors = FALSE) else NULL
}))

# Panel A: AICc — covers all 6 models (Q_6 + H(3,4) variants).
label_map <- c(
  "M1_phys"           = "M1\nPhys",
  "M2_topo"           = "M2\nTopo (Q₆)",
  "M3_phys_topo"      = "M3\nPhys+Topo (Q₆)",
  "M4_full"           = "M4\n+ tRNA",
  "M2_topo_k43"       = "M2\nTopo H(3,4)",
  "M3_phys_topo_k43"  = "M3\nPhys+Topo H(3,4)"
)
# Drop any rows whose Model isn't in the label_map (shouldn't happen, but guards
# against silent NA labels if the CSV adds new models without corresponding labels).
model_comp <- model_comp[model_comp$Model %in% names(label_map), ]
model_comp$Model_label <- label_map[model_comp$Model]
# Keep the same row order Python wrote (sorted by AICc ascending), then reverse
# for coord_flip so the best model appears at the top of the panel.
model_comp$Model_label <- factor(model_comp$Model_label,
                                 levels = rev(model_comp$Model_label))
model_comp$is_best <- model_comp$Delta_AICc == 0

p5a <- ggplot(model_comp, aes(x = Model_label, y = AICc, fill = is_best)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = paste0("Δ=", sprintf("%.1f", Delta_AICc))),
            hjust = -0.1, size = ANNOT_SIZE_S, color = "grey30") +
  scale_fill_manual(values = c("TRUE" = PAL_BLUE, "FALSE" = PAL_BLUE_PALE), guide = "none") +
  coord_flip(ylim = c(min(model_comp$AICc) * 0.97, max(model_comp$AICc) * 1.05)) +
  labs(x = NULL, y = "AICc") +
  tc() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.2),
        axis.text.y = element_text(size = 7))

# Panel B: Percentile ranks
ranks_m1 <- all_ranks[all_ranks$Model == "M1_phys", ]
ranks_m3 <- all_ranks[all_ranks$Model == "M3_phys_topo", ]
ranks_combined <- rbind(
  data.frame(Model = "M1 (Phys)", Percentile = ranks_m1$Percentile),
  data.frame(Model = "M3 (Phys+Topo)", Percentile = ranks_m3$Percentile)
)

p5b <- ggplot(ranks_combined, aes(x = Percentile, fill = Model)) +
  geom_histogram(binwidth = 10, position = "dodge",
                 color = "grey30", linewidth = 0.15, alpha = 0.85) +
  scale_fill_manual(values = c("M1 (Phys)" = PAL_BLUE_PALE, "M3 (Phys+Topo)" = PAL_BLUE)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Percentile rank", y = "Count", fill = NULL) +
  tc() +
  theme(legend.position = c(0.22, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.25, "cm"))

# Panel C: LRT — read all 5 LR tests from the CSV (3 Q_6 + 2 H(3,4) variants).
lrt_csv <- read.csv("output/tables/T_likelihood_ratio_tests.csv", stringsAsFactors = FALSE)
# Map (Restricted, Full) to a short label.
lrt_label_map <- c(
  "M1_phys|M3_phys_topo"               = "M1→M3",
  "M2_topo|M3_phys_topo"               = "M2→M3",
  "M3_phys_topo|M4_full"               = "M3→M4",
  "M1_phys|M3_phys_topo_k43"           = "M1→M3_H",
  "M2_topo_k43|M3_phys_topo_k43"       = "M2_H→M3_H"
)
lrt_csv$key <- paste(lrt_csv$Restricted, lrt_csv$Full, sep = "|")
lrt_csv$Test_label <- lrt_label_map[lrt_csv$key]
lrt_csv <- lrt_csv[!is.na(lrt_csv$Test_label), ]
lrt_csv$Test_label <- factor(lrt_csv$Test_label, levels = unname(lrt_label_map))
lrt_csv$is_sig <- lrt_csv$Significant_p05
# Pretty-print p-value
fmt_p <- function(p_str) {
  p <- suppressWarnings(as.numeric(p_str))
  if (is.na(p)) p_str
  else if (p < 1e-10) "<10⁻¹⁰"
  else if (p < 0.001) sprintf("%.1e", p)
  else sprintf("%.2f", p)
}
lrt_csv$p_pretty <- vapply(as.character(lrt_csv$p_value), fmt_p, character(1))
# Grammatical p-annotation: for strict inequalities (fmt_p returns "<...") use
# "p < ..."; otherwise use "p = ..." — avoids the ungrammatical "p=<0.001".
lrt_csv$annot <- paste0(
  ifelse(lrt_csv$is_sig, "***", "n.s."), "\n",
  ifelse(startsWith(lrt_csv$p_pretty, "<"),
         paste0("p < ", substring(lrt_csv$p_pretty, 2)),
         paste0("p = ", lrt_csv$p_pretty))
)

p5c <- ggplot(lrt_csv, aes(x = Test_label, y = LR_statistic, fill = is_sig)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = annot), vjust = -0.2, size = 2.3, fontface = "bold",
            lineheight = 0.85) +
  scale_fill_manual(values = c("TRUE" = PAL_RED, "FALSE" = PAL_GREY_LT), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "Likelihood ratio") +
  tc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

fig5 <- p5a / (p5b | p5c) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_figure(fig5, file.path(output_dir, "Fig5_condlogit"), width = 7, height = 6)


cat("\nAll 5 manuscript figures saved to:", output_dir, "\n")
