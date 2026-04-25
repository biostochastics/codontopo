#!/usr/bin/env Rscript
# Supplementary figure: Q_6 topology-avoidance across all 24 base-to-bit encodings.
#
# Produces FigS_encoding_sweep showing per-encoding observed vs candidate-landscape
# topology-breaking rates and hypergeometric p-values, demonstrating that 8 of 24
# encodings give no Q_6 depletion (motivating H(3,4) as the primary adjacency).
#
# Inputs:
#   output/topology_avoidance.json (Q6_encoding_sweep block)
#
# Outputs:
#   output/figures/FigS_encoding_sweep.{png,pdf}

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(dplyr)
})

# Locate theme
script_dir <- if (length(sys.frames()) >= 1) {
  dirname(sys.frame(1)$ofile %||% commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))[1]] %||% ".")
} else "."
theme_path <- file.path(script_dir, "theme_codon.R")
if (!file.exists(theme_path)) {
  theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
}
if (file.exists(theme_path)) source(theme_path)

args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "output"
output_dir <- if (length(args) >= 2) args[2] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating supplementary encoding-sweep figure...\n")

j <- fromJSON(file.path(input_dir, "topology_avoidance.json"))
sweep_rows <- j$Q6_encoding_sweep$rows
df <- as.data.frame(sweep_rows)
df$is_default <- df$encoding_index == 0
df$signif <- df$hypergeom_p < 0.05

# Sort by depletion fold for visual order
df$encoding_label <- factor(df$encoding_label, levels = df$encoding_label[order(df$depletion_fold)])

df$p_label <- vapply(df$hypergeom_p, function(p) {
  if (p < 1e-10) "<10^-10"
  else if (p < 0.001) sprintf("%.0e", p)
  else if (p < 0.01) sprintf("%.3f", p)
  else sprintf("%.2f", p)
}, character(1))

p1 <- ggplot(df, aes(x = encoding_label, y = depletion_fold, fill = signif)) +
  geom_col(color = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
  # p-value labels in plot reading orientation (horizontal under coord_flip);
  # placed just past the bar end with a small fixed offset.
  geom_text(aes(label = p_label),
            hjust = -0.15, size = 2.4, family = "sans") +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#3182bd"),
                    name = "p < 0.05") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Base-to-bit encoding (sorted by Q_6 depletion fold)",
       y = "Q_6 depletion fold (rate_possible / rate_observed)",
       title = "Q_6 topology-avoidance is encoding-sensitive",
       subtitle = "8/24 encodings give no significant depletion; H(3,4) is encoding-independent and is the primary test") +
  coord_flip() +
  if (exists("theme_codon_pub")) theme_codon_pub(base_size = 9) else theme_minimal(base_size = 9)

# Save
ggsave(file.path(output_dir, "FigS_encoding_sweep.png"), p1, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "FigS_encoding_sweep.pdf"), p1, width = 7, height = 5.5)
cat("  Wrote FigS_encoding_sweep.{png,pdf}\n")

# Also export a tidy CSV for the table
write.csv(df, file.path(output_dir, "../tables/TS_encoding_sweep.csv"), row.names = FALSE)
cat("  Wrote TS_encoding_sweep.csv\n")
