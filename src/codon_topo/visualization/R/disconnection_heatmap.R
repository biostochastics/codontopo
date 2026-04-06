# Disconnection catalogue heatmap across NCBI tables
# Input: CSV from codon_topo.visualization.data_export.export_disconnection_catalogue()
# Output: Publication-quality heatmap (ggplot2 + ggpubr, 300 DPI)
#
# Usage:
#   Rscript disconnection_heatmap.R disconnection_catalogue.csv output.pdf

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript disconnection_heatmap.R <input.csv> <output.pdf>")

input_file <- args[1]
output_file <- args[2]

df <- read.csv(input_file, stringsAsFactors = FALSE)

p <- ggplot(df, aes(x = aa, y = factor(table_id), fill = reconnect_eps)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_components), size = 3) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  labs(
    title = "Amino Acid Disconnections Across NCBI Translation Tables",
    x = "Amino Acid",
    y = "NCBI Table ID",
    fill = expression("Reconnection" ~ epsilon)
  ) +
  theme_pubr(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_file, p, width = 8, height = 10, dpi = 300)
cat("Saved:", output_file, "\n")
