#!/usr/bin/env Rscript
# Disconnection catalogue heatmap (standalone)
# Input: CSV from codon_topo.visualization.data_export.export_disconnection_catalogue()
# Output: Publication-quality heatmap
#
# Usage:
#   Rscript disconnection_heatmap.R disconnection_catalogue.csv output.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
if (!file.exists(theme_path)) {
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  theme_path <- file.path(script_dir, "theme_codon.R")
}
source(theme_path)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript disconnection_heatmap.R <input.csv> <output.pdf>")

input_file  <- args[1]
output_file <- args[2]

df <- read.csv(input_file, stringsAsFactors = FALSE)

p <- ggplot(df, aes(x = aa, y = factor(table_id), fill = reconnect_eps)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_components), size = 3, fontface = "bold",
            color = ifelse(df$reconnect_eps >= 3.5, "white", "grey15")) +
  scale_fill_gradient2(
    low = PAL_ORANGE, mid = PAL_RED_LT, high = PAL_PURPLE,
    midpoint = 3, na.value = "grey90",
    name = expression("Reconnection " * epsilon)
  ) +
  labs(
    title = "Amino Acid Disconnections Across NCBI Translation Tables",
    x = "Amino acid",
    y = "NCBI table ID"
  ) +
  theme_codon_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_blank())

ggsave(output_file, p, width = 8, height = 10, dpi = DPI)
cat("Saved:", output_file, "\n")
