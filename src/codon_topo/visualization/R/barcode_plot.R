#!/usr/bin/env Rscript
# Persistent homology barcode plot (standalone)
# Input: CSV from codon_topo.visualization.data_export.export_persistent_homology()
# Output: Publication-quality barcode plot
#
# Usage:
#   Rscript barcode_plot.R persistent_homology.csv output.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(viridis)
})

theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
if (!file.exists(theme_path)) {
  script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  theme_path <- file.path(script_dir, "theme_codon.R")
}
source(theme_path)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript barcode_plot.R <input.csv> <output.pdf>")

input_file  <- args[1]
output_file <- args[2]

df <- read.csv(input_file, stringsAsFactors = FALSE)

# Filter to AAs that are ever disconnected (beta_0 > 1)
disconnected <- unique(df$aa[df$beta_0 > 1])
df_disc <- df[df$aa %in% disconnected, ]

p <- ggplot(df_disc, aes(x = epsilon, y = beta_0, color = aa, group = aa)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:6) +
  scale_color_viridis_d(option = "D", end = 0.85) +
  labs(
    title = "Persistent Homology: Disconnected Amino Acids",
    x = expression(epsilon ~ "(Hamming distance threshold)"),
    y = expression(beta[0] ~ "(connected components)"),
    color = "Amino acid"
  ) +
  theme_codon_pub()

ggsave(output_file, p, width = 8, height = 5, dpi = DPI)
cat("Saved:", output_file, "\n")
