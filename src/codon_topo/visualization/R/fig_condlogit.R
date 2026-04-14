#!/usr/bin/env Rscript
# ============================================================
# Figure: Conditional logit model comparison diagnostics
# Panel A: AICc comparison (horizontal bar chart)
# Panel B: Observed move percentile ranks (histogram, M1 vs M3)
# Panel C: Likelihood ratio tests (bar chart)
#
# Clayworth & Kornilov 2026
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

# Source shared theme
theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
if (!file.exists(theme_path)) {
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  theme_path <- file.path(script_dir, "theme_codon.R")
}
source(theme_path)

# ── Load data ──
model_comp <- read.csv("output/tables/T_model_comparison.csv",
                       stringsAsFactors = FALSE)

rank_files <- c("M1_phys", "M2_topo", "M3_phys_topo", "M4_full")
all_ranks <- do.call(rbind, lapply(rank_files, function(m) {
  fp <- paste0("output/tables/T_ranks_", m, ".csv")
  if (file.exists(fp)) read.csv(fp, stringsAsFactors = FALSE)
  else NULL
}))

# ── Panel A: AICc comparison ──
# Map labels by model ID column, not row order
label_map <- c(
  "M1_phys" = "M1\nPhys",
  "M2_topo" = "M2\nTopo",
  "M3_phys_topo" = "M3\nPhys+Topo",
  "M4_full" = "M4\nFull"
)
if ("Model" %in% names(model_comp)) {
  model_comp$Model_label <- label_map[model_comp$Model]
} else {
  model_comp$Model_label <- c("M3\nPhys+Topo", "M4\nFull", "M2\nTopo", "M1\nPhys")
}
model_comp$Model_label <- factor(model_comp$Model_label,
                                  levels = rev(unname(label_map)))
model_comp$is_best <- model_comp$Delta_AICc == 0

pA <- ggplot(model_comp, aes(x = Model_label, y = AICc, fill = is_best)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = paste0("d=", sprintf("%.1f", Delta_AICc))),
            hjust = -0.12, size = ANNOT_SIZE, color = "grey30") +
  scale_fill_manual(values = c("TRUE" = PAL_BLUE, "FALSE" = PAL_BLUE_PALE),
                    guide = "none") +
  coord_flip(ylim = c(min(model_comp$AICc) * 0.97,
                       max(model_comp$AICc) * 1.02)) +
  labs(x = NULL, y = "AICc", title = "A. Model comparison (AICc)") +
  theme_codon_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3))

# ── Panel B: Percentile rank distribution (M3 vs M1) ──
ranks_m1 <- all_ranks[all_ranks$Model == "M1_phys", ]
ranks_m3 <- all_ranks[all_ranks$Model == "M3_phys_topo", ]
ranks_combined <- rbind(
  data.frame(Model = "M1 (Phys only)", Percentile = ranks_m1$Percentile),
  data.frame(Model = "M3 (Phys+Topo)", Percentile = ranks_m3$Percentile)
)

pB <- ggplot(ranks_combined, aes(x = Percentile, fill = Model)) +
  geom_histogram(binwidth = 10, position = "dodge",
                 color = "grey30", linewidth = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("M1 (Phys only)" = PAL_BLUE_PALE,
                                "M3 (Phys+Topo)" = PAL_BLUE)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("text", x = 52, y = Inf, label = "Chance", vjust = 1.5,
           hjust = 0, size = ANNOT_SIZE, color = "grey50") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Percentile rank of observed move",
       y = "Count",
       title = "B. Move prediction quality",
       fill = NULL) +
  theme_codon_pub() +
  theme(legend.position = c(0.22, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"))

# ── Panel C: LRT summary ──
lrt_data <- data.frame(
  Test = c("M1 -> M3\n(+Topo)", "M2 -> M3\n(+Phys)", "M3 -> M4\n(+tRNA)"),
  LR = c(102.2, 68.8, 0.15),
  p_val = c("<0.001", "<0.001", "0.70"),
  Sig = c("***", "***", "n.s.")
)
lrt_data$Test <- factor(lrt_data$Test, levels = lrt_data$Test)
lrt_data$is_sig <- c(TRUE, TRUE, FALSE)
lrt_data$label <- paste0(lrt_data$Sig, "\np=", lrt_data$p_val)

pC <- ggplot(lrt_data, aes(x = Test, y = LR, fill = is_sig)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = label), vjust = -0.3, size = 3.2, fontface = "bold",
            lineheight = 0.85) +
  scale_fill_manual(values = c("TRUE" = PAL_RED, "FALSE" = PAL_GREY_LT),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Likelihood ratio", title = "C. Nested model tests (LRT)") +
  theme_codon_pub()

# ── Compose ──
fig <- pA / (pB | pC) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

save_figure(fig, "output/figures/FigG_condlogit_diagnostics", width = 7, height = 7)

cat("Conditional logit diagnostics figure saved.\n")
