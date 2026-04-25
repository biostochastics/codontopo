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

# Panel A: Null distribution
mc <- fromJSON(file.path(table_dir, "T3_coloring_optimality.json"))
obs <- mc$observed_score; mu <- mc$null_mean; sigma <- mc$null_std
x_range <- seq(mu - 4*sigma, mu + 4*sigma, length.out = 400)
df_null <- data.frame(score = x_range, density = dnorm(x_range, mu, sigma))
df_shade <- df_null %>% filter(score <= obs)

p2a <- ggplot(df_null, aes(x = score, y = density)) +
  geom_area(fill = PAL_BLUE_PALE, alpha = 0.6) +
  geom_area(data = df_shade, fill = PAL_BLUE, alpha = 0.25) +
  geom_line(color = PAL_BLUE, linewidth = 0.7) +
  geom_vline(xintercept = obs, color = PAL_RED, linewidth = 0.9) +
  annotate("label", x = obs - 60, y = max(df_null$density) * 0.85,
           label = paste0("Obs: ", format(round(obs), big.mark = ",")),
           color = PAL_RED, fill = "white", linewidth = 0.2,
           size = ANNOT_SIZE_S, fontface = "bold", hjust = 1) +
  annotate("label", x = mu + 0.5*sigma, y = max(df_null$density) * 0.6,
           label = paste0("p = ", format(mc$p_value_conservative, digits = 2)),
           color = PAL_BLUE, fill = "white", linewidth = 0.2,
           size = ANNOT_SIZE_S, hjust = 0) +
  labs(x = expression(italic(F) * "(code)"), y = "Density") +
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

# Panel A: Bit-position bias
bit_data <- data.frame(
  bit = paste0("Bit ", 0:5),
  position = factor(c("1st", "1st", "2nd", "2nd", "3rd (wob)", "3rd (wob)"),
                    levels = c("1st", "2nd", "3rd (wob)")),
  count = c(5, 4, 0, 8, 13, 5)
)
bit_data$bit <- factor(bit_data$bit, levels = paste0("Bit ", 0:5))
expected <- sum(bit_data$count) / 6

p3a <- ggplot(bit_data, aes(x = bit, y = count, fill = position)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = expected, linetype = "dashed", color = PAL_RED, linewidth = 0.4) +
  geom_text(aes(label = count), vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("1st" = PAL_PURPLE, "2nd" = PAL_BLUE, "3rd (wob)" = PAL_TEAL),
                    name = "Pos") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = expression("GF(2)"^6 ~ "bit"), y = "Bit changes") +
  tc() +
  theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(size = 7))

# Panel B: Depth calibration
depth <- read.csv(file.path(input_dir, "depth_calibration.csv"))

p3b <- ggplot(depth, aes(x = age_midpoint_mya, y = reconnect_eps,
                           color = aa, shape = aa)) +
  # Horizontal CI segments rendered manually so the right-most Ser (Pre-LUCA,
  # tight log-scale CI) is still visible above the marker glyph.
  geom_segment(aes(x = age_mya_low, xend = age_mya_high,
                   y = reconnect_eps, yend = reconnect_eps),
               linewidth = 0.7, alpha = 0.85,
               lineend = "butt", show.legend = FALSE) +
  geom_segment(aes(x = age_mya_low, xend = age_mya_low,
                   y = reconnect_eps - 0.18, yend = reconnect_eps + 0.18),
               linewidth = 0.7, alpha = 0.85, show.legend = FALSE) +
  geom_segment(aes(x = age_mya_high, xend = age_mya_high,
                   y = reconnect_eps - 0.18, yend = reconnect_eps + 0.18),
               linewidth = 0.7, alpha = 0.85, show.legend = FALSE) +
  geom_point(size = 3, stroke = 0.6) +
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

df_ta <- data.frame(
  graph    = factor(rep(c("H(3,4) (primary)", "Q[6] (sensitivity)"), each = 2),
                    levels = c("H(3,4) (primary)", "Q[6] (sensitivity)")),
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

p3c <- ggplot(df_ta, aes(x = category_label, y = rate, fill = category)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(rate, 1), "%")),
            vjust = -0.4, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("Observed" = PAL_BLUE, "Candidate" = PAL_GREY_LT)) +
  scale_y_continuous(limits = c(0, 95), expand = expansion(mult = c(0, 0.10))) +
  facet_wrap(~ graph, scales = "free_x") +
  labs(x = NULL, y = "% topology-breaking moves",
       subtitle = sprintf("H(3,4) RR=%.2f, p<1e-4   Q_6 %.1f-fold depl.",
                          tk43$risk_ratio,
                          tq6$rate_possible / max(tq6$rate_observed, 0.001))) +
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
trna$label <- paste0(trna$reassigned_aa,
                     " (", gsub(" vs .*", "", trna$pairing), ")")
trna$exact_p_clamped <- pmin(trna$exact_p, 0.3)

p3d <- ggplot(trna, aes(x = reorder(label, aa_rank), y = aa_rank,
                          fill = exact_p_clamped)) +
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
cat_data$status <- factor(cat_data$status,
  levels = c("verified", "tested", "pending", "null"))
cat_data$workstream <- factor(cat_data$workstream,
  levels = c("WS1", "WS2", "WS3", "WS4", "WS6"))

p4c <- ggplot(cat_data, aes(x = workstream, fill = status)) +
  geom_bar(width = 0.55, color = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = PAL_STATUS,
                    labels = c("Verified", "Tested", "Pending", "Null")) +
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
  "M1_phys|M3_phys_topo"               = "M1→M3\n(+Q₆ topo)",
  "M2_topo|M3_phys_topo"               = "M2→M3\n(+phys)",
  "M3_phys_topo|M4_full"               = "M3→M4\n(+tRNA)",
  "M1_phys|M3_phys_topo_k43"           = "M1→M3 H(3,4)\n(+H(3,4) topo)",
  "M2_topo_k43|M3_phys_topo_k43"       = "M2→M3 H(3,4)\n(+phys)"
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
lrt_csv$annot <- paste0(ifelse(lrt_csv$is_sig, "***", "n.s."), "\np=", lrt_csv$p_pretty)

p5c <- ggplot(lrt_csv, aes(x = Test_label, y = LR_statistic, fill = is_sig)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = annot), vjust = -0.2, size = 2.3, fontface = "bold",
            lineheight = 0.85) +
  scale_fill_manual(values = c("TRUE" = PAL_RED, "FALSE" = PAL_GREY_LT), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "Likelihood ratio") +
  tc() +
  theme(axis.text.x = element_text(size = 6.5))

fig5 <- p5a / (p5b | p5c) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_figure(fig5, file.path(output_dir, "Fig5_condlogit"), width = 7, height = 6)


cat("\nAll 5 manuscript figures saved to:", output_dir, "\n")
