#!/usr/bin/env Rscript
# Supplementary figure: 2-adic Walsh spectral-depth null + 24-encoding invariance
# sweep. Companion to §S16 (Walsh-Hadamard / 2-adic spectral probe) in the
# supplement.
#
# Panel A: block-size-matched null distribution of spectral depth with the
#          standard code's observed value (=544) marked and z-score annotated.
#          If per-draw null samples are not persisted in walsh_2adic.json we
#          fall back to a normal-approximation density (mu, sigma from the
#          test2_block_size_null block); with z = -17.74 the observed marker
#          sits far in the tail either way.
# Panel B: 24-encoding z-score sweep summary showing the tight cluster of
#          z-scores in [-18.9, -16.9], demonstrating that the null rejection
#          is invariant across all 24 base-to-bit bijections.
#
# Inputs:
#   output/walsh_2adic.json
#     - test2_block_size_null.{null_mean, null_std, observed_depth, z_score, n_null}
#       (optionally null_samples[])
#     - test3_encoding_invariance.{z_score_min, z_score_max, z_score_mean,
#                                  n_encodings, n_null_per_encoding}
#       (optionally z_scores[])
#
# Outputs:
#   output/figures/FigW_walsh_spectrum.{png,pdf}
#
# R1 response: R1 explicitly asked for engagement with the 2-adic framework.
# The block-size null z = -17.74 quantifies the descriptive bridge; the
# encoding-invariance panel shows the bridge is a property of the code, not
# an artefact of one particular base-to-bit bijection.

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

# Locate theme (works from repo root or when sourced)
`%||%` <- function(a, b) if (!is.null(a)) a else b
script_frame <- tryCatch(sys.frame(1), error = function(e) NULL)
script_dir <- if (!is.null(script_frame) && !is.null(script_frame$ofile)) {
  dirname(script_frame$ofile)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  fa   <- grep("^--file=", args, value = TRUE)
  if (length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else "."
}
theme_path <- file.path(script_dir, "theme_codon.R")
if (!file.exists(theme_path)) {
  theme_path <- "src/codon_topo/visualization/R/theme_codon.R"
}
if (!file.exists(theme_path)) {
  stop("Cannot locate theme_codon.R (looked at: ", theme_path, ")")
}
source(theme_path)

args        <- commandArgs(trailingOnly = TRUE)
input_dir   <- if (length(args) >= 1) args[1] else "output"
output_dir  <- if (length(args) >= 2) args[2] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating FigW_walsh_spectrum ...\n")

json_path <- file.path(input_dir, "walsh_2adic.json")
if (!file.exists(json_path)) {
  stop("walsh_2adic.json not found at ", json_path)
}
w <- fromJSON(json_path, simplifyVector = TRUE)

# ── Panel A: block-size null ──────────────────────────────────
b       <- w$test2_block_size_null
mu      <- b$null_mean
sigma   <- b$null_std
obs     <- b$observed_depth
z_val   <- b$z_score
n_null  <- b$n_null %||% 2000
have_samples <- !is.null(b$null_samples) && length(b$null_samples) > 0

# X range: include the observed value plus generous padding around the null
x_lo <- min(obs - 5, mu - 4 * sigma)
x_hi <- max(obs + 5, mu + 4 * sigma)

if (have_samples) {
  null_df <- data.frame(spectral_depth = as.numeric(b$null_samples))
  pA_base <- ggplot(null_df, aes(x = spectral_depth)) +
    geom_histogram(bins = 40, fill = PAL_BLUE, colour = "white",
                   alpha = 0.75, linewidth = 0.2) +
    labs(y = "Count",
         subtitle = sprintf("Empirical histogram (n = %s null draws)",
                            format(n_null, big.mark = ",")))
  y_max <- max(hist(null_df$spectral_depth,
                    breaks = 40, plot = FALSE)$counts)
} else {
  x_range <- seq(x_lo, x_hi, length.out = 400)
  dens_df <- data.frame(
    spectral_depth = x_range,
    density        = dnorm(x_range, mu, sigma)
  )
  pA_base <- ggplot(dens_df, aes(x = spectral_depth, y = density)) +
    geom_area(fill = PAL_BLUE, alpha = 0.35) +
    geom_line(colour = PAL_BLUE, linewidth = 0.7) +
    labs(y = "Density",
         subtitle = sprintf("Normal approximation (n = %s null draws)",
                            format(n_null, big.mark = ",")))
  y_max <- max(dens_df$density)
}

# Observed marker + null-summary annotation. Anchor labels to a fraction of
# y_max so they render above the geoms in both panels.
label_y_obs <- y_max * 1.02
label_y_mu  <- y_max * 0.98

pA <- pA_base +
  geom_vline(xintercept = obs, colour = PAL_RED, linewidth = 1.1) +
  annotate("segment",
           x = obs, xend = mu, y = y_max * 0.55, yend = y_max * 0.55,
           colour = PAL_GREY, linewidth = 0.4,
           arrow = grid::arrow(length = unit(0.15, "cm"), ends = "both")) +
  annotate("text",
           x = (obs + mu) / 2, y = y_max * 0.62,
           label = sprintf("z = %.2f", z_val),
           colour = PAL_GREY, size = ANNOT_SIZE, fontface = "bold") +
  annotate("label",
           x = obs, y = label_y_obs,
           label = sprintf("Observed\n= %d", obs),
           colour = PAL_RED, fill = "white",
           hjust = -0.05, vjust = 1, size = ANNOT_SIZE, fontface = "bold") +
  annotate("label",
           x = mu, y = label_y_mu,
           label = sprintf("Null:\n%.1f +/- %.1f", mu, sigma),
           colour = PAL_BLUE, fill = "white",
           hjust = 0.5, vjust = 1, size = ANNOT_SIZE) +
  scale_x_continuous(limits = c(x_lo, x_hi), expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "A. 2-adic spectral depth vs block-size-matched null",
       x = expression("Walsh spectral depth (" * Sigma * " " * v[2] *
                      " of Walsh coefficients)")) +
  theme_codon_pub()

# ── Panel B: 24-encoding z-score sweep ────────────────────────
enc <- w$test3_encoding_invariance
z_min <- enc$z_score_min
z_max <- enc$z_score_max
z_med <- enc$z_score_mean   # JSON stores mean, not median; label accordingly
n_enc <- enc$n_encodings %||% 24
n_null_per <- enc$n_null_per_encoding %||% 1500
have_z_vector <- !is.null(enc$z_scores) && length(enc$z_scores) > 0

if (have_z_vector) {
  enc_df <- data.frame(
    encoding = seq_along(enc$z_scores),
    z        = as.numeric(enc$z_scores)
  )
  pB <- ggplot(enc_df, aes(x = z, y = 0)) +
    geom_jitter(width = 0, height = 0.15, size = 2.6,
                shape = 21, fill = PAL_TEAL, colour = "grey20",
                alpha = 0.85, stroke = 0.3) +
    geom_vline(xintercept = z_med, colour = PAL_GREY, linetype = "dashed",
               linewidth = 0.4) +
    annotate("text", x = z_med, y = 0.42,
             label = sprintf("mean = %.2f", z_med),
             colour = PAL_GREY, size = ANNOT_SIZE_S, hjust = 0.5) +
    scale_y_continuous(limits = c(-0.5, 0.5), breaks = NULL) +
    labs(y = NULL)
} else {
  # Fallback: 3-point min/mean/max summary + range bar
  enc_df <- data.frame(
    label = factor(c("min", "mean", "max"), levels = c("min", "mean", "max")),
    z     = c(z_min, z_med, z_max),
    shape = c("min", "mean", "max")
  )
  # All three labels placed above their markers with a callout leader; markers
  # sit on the range bar (y=0). Values above, "min"/"mean"/"max" tags below.
  enc_df$y_pt   <- 0
  enc_df$y_lab  <- c(0.55, 0.30, 0.55)   # min & max lifted higher to avoid overlap

  pB <- ggplot(enc_df, aes(x = z, y = y_pt)) +
    annotate("segment",
             x = z_min, xend = z_max, y = 0, yend = 0,
             colour = PAL_TEAL, linewidth = 6, alpha = 0.25) +
    annotate("text",
             x = (z_min + z_max) / 2, y = -0.28,
             label = sprintf("range across %d encodings", n_enc),
             colour = "grey30", size = ANNOT_SIZE_S - 0.3,
             fontface = "italic") +
    geom_segment(aes(x = z, xend = z, y = 0, yend = y_lab - 0.05),
                 colour = "grey60", linewidth = 0.3) +
    geom_point(aes(shape = label, fill = label), size = 4.5,
               colour = "grey20", stroke = 0.5) +
    geom_label(aes(y = y_lab, label = sprintf("%s\n%.2f", label, z),
                   fill = label, colour = label),
               size = ANNOT_SIZE_S, fontface = "bold",
               label.size = 0.3, alpha = 0.9) +
    scale_shape_manual(values = c(min = 21, mean = 22, max = 21),
                       guide = "none") +
    scale_fill_manual(values = c(min = PAL_BLUE_PALE,
                                 mean = "#FDE0DE",
                                 max  = PAL_BLUE_PALE),
                      guide = "none") +
    scale_colour_manual(values = c(min = PAL_BLUE,
                                   mean = PAL_RED,
                                   max  = PAL_BLUE),
                        guide = "none") +
    scale_y_continuous(limits = c(-0.6, 0.9), breaks = NULL) +
    labs(y = NULL)
}

pB <- pB +
  # Reference vertical at z = 0 for orientation
  geom_vline(xintercept = 0, colour = "grey60",
             linetype = "dotted", linewidth = 0.35) +
  # p < 0.001 reference at z = -3.29 (two-sided) for context; the cluster sits
  # far to the left of this, so annotate it clearly.
  geom_vline(xintercept = -3.29, colour = PAL_RED,
             linetype = "dotdash", linewidth = 0.35, alpha = 0.6) +
  annotate("text", x = -3.29, y = -0.45,
           label = "p = 0.001\n(z = -3.29)",
           colour = PAL_RED, size = ANNOT_SIZE_S, hjust = -0.05, vjust = 1) +
  scale_x_continuous(
    limits = c(min(-5, z_min - 3), 1.5),
    breaks = c(0, -3.29, -5, -10, -15, -20)
  ) +
  labs(
    title = sprintf(
      "B. Encoding invariance: %d base-to-bit bijections", n_enc
    ),
    subtitle = sprintf(
      "z-score against block-size null (n = %s per encoding); range [%.2f, %.2f], mean %.2f",
      format(n_null_per, big.mark = ","), z_min, z_max, z_med
    ),
    x = "z-score (spectral depth vs block-size-matched null)"
  ) +
  theme_codon_pub() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())

# ── Combine and save ──────────────────────────────────────────
p_final <- pA / pB + plot_layout(heights = c(1.5, 1))

out_stem <- file.path(output_dir, "FigW_walsh_spectrum")
save_figure(p_final, out_stem, width = 7.2, height = 7.0)

cat("Done.\n")
