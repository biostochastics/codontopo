# C^3 embedding scatter plot
# Input: CSV from codon_topo.visualization.data_export.export_embedding_coords()
# Output: Publication-quality scatter plot (ggplot2 + ggpubr, 300 DPI)
#
# Usage:
#   Rscript embedding_scatter.R embedding_coords.csv output.pdf

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript embedding_scatter.R <input.csv> <output.pdf>")

input_file <- args[1]
output_file <- args[2]

df <- read.csv(input_file, stringsAsFactors = FALSE)
df <- df[df$aa != "Stop", ]

# Project C^3 to 2D using first two coordinates (Re(z1), Re(z2))
p <- ggplot(df, aes(x = z1_re, y = z2_re, color = aa, label = codon)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(size = 2, vjust = -1, show.legend = FALSE) +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = expression("Codon Embedding" ~ phi * ": GF(2)"^6 %->% C^3),
    x = expression("Re(" * z[1] * ")"),
    y = expression("Re(" * z[2] * ")"),
    color = "Amino Acid"
  ) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "right")

ggsave(output_file, p, width = 10, height = 8, dpi = 300)
cat("Saved:", output_file, "\n")
