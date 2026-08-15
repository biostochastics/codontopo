#!/usr/bin/env Rscript
# Supplementary figure: Slavov (Tsour et al. 2026) SAAP codon-distance
# grouped bars. Companion to §S17 of the supplement.
#
# R3 explicitly asked for a *visual* contrast of the 65.0% vs 39.5% Hamming-1
# excess (§S17 currently ships only a 3-row table). The single-panel grouped
# bar chart lets a reviewer see the enrichment at each of distances 1, 2, 3
# in one glance.
#
# Inputs:
#   output/slavov_saap_codon_distances.json
#     - n_parseable_events (5611)
#     - n_high_confidence_saap_events (5873)
#     - event_distance_distribution.{1,2,3}.{count, fraction}
#     - baseline_distance_distribution.{1,2,3}.{count, fraction}
#     - event_enrichment_binomial_p_one_sided (single-NT test; p < 1e-100)
#
# Outputs:
#   output/figures/FigSl_slavov_saap.{png,pdf}
#
# Naming: keeps the two-letter alphabetic supplement-figure convention
# (FigA, FigB, ..., FigW_walsh_spectrum). "Sl" == Slavov.

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

# Locate theme (works from repo root or when sourced from anywhere)
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

args       <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "output"
output_dir <- if (length(args) >= 2) args[2] else "output/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating FigSl_slavov_saap ...\n")

json_path <- file.path(input_dir, "slavov_saap_codon_distances.json")
if (!file.exists(json_path)) {
  stop("slavov_saap_codon_distances.json not found at ", json_path)
}
s <- fromJSON(json_path, simplifyVector = TRUE)

# ── Extract per-distance data (data-driven from JSON) ─────────
if (is.null(s$event_distance_distribution) ||
    is.null(s$baseline_distance_distribution)) {
  stop("JSON is missing event_distance_distribution / ",
       "baseline_distance_distribution. Re-run scripts/slavov_saap_analysis.py ",
       "to regenerate.")
}
distances <- sort(as.integer(names(s$event_distance_distribution)))
obs_frac  <- sapply(as.character(distances),
                    function(d) s$event_distance_distribution[[d]]$fraction)
obs_count <- sapply(as.character(distances),
                    function(d) s$event_distance_distribution[[d]]$count)
base_frac <- sapply(as.character(distances),
                    function(d) s$baseline_distance_distribution[[d]]$fraction)
base_count <- sapply(as.character(distances),
                     function(d) s$baseline_distance_distribution[[d]]$count)

n_events  <- s$n_parseable_events
n_hc      <- s$n_high_confidence_saap_events %||% n_events
n_pairs   <- s$n_all_pairs_baseline %||% sum(base_count)
p_binom   <- s$event_enrichment_binomial_p_one_sided

# Sanity checks (event counts must sum to n_parseable_events)
if (abs(sum(obs_count) - n_events) > 0) {
  warning("Event counts across distances (", sum(obs_count),
          ") do not sum to n_parseable_events (", n_events, ")")
}
if (sum(base_count) != n_pairs) {
  warning("Baseline pair counts (", sum(base_count),
          ") do not sum to n_all_pairs_baseline (", n_pairs, ")")
}

# ── Assemble long-format df for grouped bars ──────────────────
lvl_obs  <- sprintf("Tsour/Slavov SAAP events (n = %s)",
                    format(n_events, big.mark = ","))
lvl_base <- sprintf("All %d AA-to-AA codon pairs (baseline)", n_pairs)

df <- data.frame(
  distance = rep(distances, 2),
  source   = factor(rep(c(lvl_obs, lvl_base), each = length(distances)),
                    levels = c(lvl_obs, lvl_base)),
  pct      = c(obs_frac, base_frac) * 100,
  count    = c(obs_count, base_count)
)

# Per-bar annotation: event counts for observed, pair counts for baseline
df$label <- ifelse(
  df$source == lvl_obs,
  sprintf("n = %s\n(%.1f%%)", format(df$count, big.mark = ","), df$pct),
  sprintf("%d pairs\n(%.1f%%)", df$count, df$pct)
)

# Format p-value for subtitle: cap at < 10^-100 to match §S17 wording
p_label <- if (is.finite(p_binom) && p_binom > 0) {
  sprintf("binomial p = %.2g", p_binom)
} else {
  "binomial p < 1e-100"
}

# ── Plot ──────────────────────────────────────────────────────
dodge_w <- 0.75
y_max   <- max(df$pct) * 1.18   # headroom for annotations

p <- ggplot(df, aes(x = factor(distance), y = pct, fill = source)) +
  geom_col(position = position_dodge(width = dodge_w),
           width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label, y = pct + 1.5),
            position = position_dodge(width = dodge_w),
            size = ANNOT_SIZE, lineheight = 0.9, vjust = 0) +
  scale_fill_manual(values = setNames(c(PAL_BLUE, PAL_GREY_LT),
                                      c(lvl_obs, lvl_base)),
                    name = NULL) +
  scale_y_continuous(limits = c(0, y_max),
                     breaks = seq(0, 100, by = 10),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Slavov SAAP events cluster at codon-Hamming distance 1",
    subtitle = sprintf(
      "Tsour et al. 2026 (>1,000 human samples): %s high-confidence events; %s\nsingle-NT enrichment vs %.1f%% baseline (%s)",
      format(n_hc, big.mark = ","),
      sprintf("%s with unambiguous AA pairs", format(n_events, big.mark = ",")),
      base_frac[1] * 100,
      p_label
    ),
    x = "Minimum source-target codon nucleotide-Hamming distance",
    y = "Percent of events / codon pairs",
    caption = paste0(
      "Blue: observed event distribution across the ", format(n_events, big.mark = ","),
      " SAAP events. Grey: null distribution of min-NT distance across all ",
      n_pairs, " AA-to-AA pairs\n",
      "in the standard code. R3 response: visual companion to Table S17 ",
      "(figure and table read the same JSON, cannot drift)."
    )
  ) +
  theme_codon_pub() +
  theme(
    legend.position     = "top",
    legend.direction    = "horizontal",
    legend.text         = element_text(size = 9),
    plot.caption        = element_text(size = 8, colour = "grey40",
                                       hjust = 0, margin = margin(t = 10)),
    panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.3)
  )

save_figure(p, file.path(output_dir, "FigSl_slavov_saap"),
            width = 7.5, height = 5.2)
cat("Done.\n")
