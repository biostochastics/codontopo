# Persistent homology barcode plot
# Input: CSV from codon_topo.visualization.data_export.export_persistent_homology()
# Output: Publication-quality barcode plot (ggplot2 + ggpubr, 300 DPI)
#
# Usage:
#   Rscript barcode_plot.R persistent_homology.csv output.pdf

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript barcode_plot.R <input.csv> <output.pdf>")

input_file <- args[1]
output_file <- args[2]

df <- read.csv(input_file, stringsAsFactors = FALSE)

# Filter to AAs that are ever disconnected (beta_0 > 1)
disconnected <- unique(df$aa[df$beta_0 > 1])
df_disc <- df[df$aa %in% disconnected, ]

p <- ggplot(df_disc, aes(x = epsilon, y = beta_0, color = aa, group = aa)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:6) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Persistent Homology: Disconnected Amino Acids",
    x = expression(epsilon ~ "(Hamming distance threshold)"),
    y = expression(beta[0] ~ "(connected components)"),
    color = "Amino Acid"
  ) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "right")

ggsave(output_file, p, width = 8, height = 5, dpi = 300)
cat("Saved:", output_file, "\n")
