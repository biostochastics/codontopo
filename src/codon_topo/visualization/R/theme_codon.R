#!/usr/bin/env Rscript
# ============================================================
# CODON-TOPO: Shared publication theme and color palette
# Clayworth & Kornilov 2026
#
# Source this file in any figure script:
#   source("src/codon_topo/visualization/R/theme_codon.R")
#
# Provides:
#   theme_codon_pub()    — ggpubr-based theme for all figures
#   PAL_*                — named color constants
#   scale_fill_codon()   — discrete fill scale
#   scale_color_codon()  — discrete color scale
#   label_sig()          — format p-values with significance stars
#   THEME_DIR            — directory containing this file
# ============================================================

# Resolve directory of this script (works with source() and Rscript)
THEME_DIR <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback: try commandArgs
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      "src/codon_topo/visualization/R"
    }
  }
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(viridis)
})

# ── Publication constants ─────────────────────────────────────
DPI   <- 300
FIG_W <- 7      # default width (inches)
FIG_H <- 5      # default height (inches)

# ── Color palette (colorblind-safe, RdBu-derived + accents) ──
PAL_BLUE      <- "#2166AC"
PAL_BLUE_LT   <- "#92C5DE"
PAL_BLUE_PALE <- "#D1E5F0"
PAL_RED       <- "#B2182B"
PAL_RED_LT    <- "#D6604D"
PAL_ORANGE    <- "#E08214"
PAL_GOLD      <- "#DFC27D"
PAL_TEAL      <- "#35978F"
PAL_GREEN     <- "#1B7837"
PAL_PURPLE    <- "#762A83"
PAL_GREY      <- "#636363"
PAL_GREY_LT   <- "#BDBDBD"
PAL_BG        <- "white"

# Discrete palette (ordered for max perceptual separation)
PAL_DISCRETE <- c(
  PAL_BLUE, PAL_RED, PAL_TEAL, PAL_ORANGE,
  PAL_PURPLE, PAL_GREEN, PAL_GOLD, PAL_GREY
)

# Status palette (for verified/tested/pending/null)
PAL_STATUS <- c(
  "verified" = PAL_BLUE,
  "tested"   = PAL_TEAL,
  "pending"  = PAL_GOLD,
  "null"     = PAL_RED
)

# Evidence strength palette
PAL_EVIDENCE <- c(
  "Very strong" = PAL_BLUE,
  "Strong"      = PAL_TEAL,
  "Moderate"    = PAL_ORANGE,
  "Weak"        = PAL_GREY_LT
)

# Significance palette
PAL_SIG <- c("Significant" = PAL_BLUE, "n.s." = PAL_GREY_LT)

# ── Theme ─────────────────────────────────────────────────────
theme_codon_pub <- function(base_size = 11) {
  theme_pubr(base_size = base_size) %+replace%
    theme(
      # Titles
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                       hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey35",
                                       hjust = 0, margin = margin(b = 8)),
      plot.title.position = "plot",
      # Axes
      axis.title       = element_text(size = base_size, face = "plain"),
      axis.text        = element_text(size = base_size - 1, color = "grey25"),
      axis.line        = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks       = element_line(color = "grey40", linewidth = 0.3),
      # Legend
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      legend.key.size  = unit(0.45, "cm"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.position  = "right",
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      # Margins
      plot.margin      = margin(12, 12, 8, 8)
    )
}

# ── Scale helpers ─────────────────────────────────────────────
scale_fill_codon <- function(...) {
  scale_fill_manual(values = PAL_DISCRETE, ...)
}
scale_color_codon <- function(...) {
  scale_color_manual(values = PAL_DISCRETE, ...)
}

# ── Annotation helpers ────────────────────────────────────────
label_sig <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "n.s.")))
}

# Consistent annotation text size
ANNOT_SIZE   <- 3.2
ANNOT_SIZE_S <- 2.8

# Save helper: always PNG + PDF
save_figure <- function(p, path_stem, width = FIG_W, height = FIG_H) {
  ggsave(paste0(path_stem, ".png"), p, width = width, height = height,
         dpi = DPI, bg = "white")
  ggsave(paste0(path_stem, ".pdf"), p, width = width, height = height,
         bg = "white")
  cat("  Saved:", paste0(path_stem, ".{png,pdf}"), "\n")
}
