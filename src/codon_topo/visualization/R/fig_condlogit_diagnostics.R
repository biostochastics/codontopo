#!/usr/bin/env Rscript
# ============================================================
# CODON-TOPO: Fig G — Conditional-logit diagnostics (supplement §S21.5)
# Clayworth & Kornilov 2026
#
# Companion to main-text Fig 5. Fig 5 reports AICc comparison +
# observed-move percentile-rank distribution + LRT — so this figure
# must NOT duplicate that content. Panels here:
#   A. Delta_phys vs Delta_topo joint distribution (with Spearman rho)
#   B. M3/M4 coefficient estimates with 95% Wald CIs (forest plot)
#   C. Posterior-predictive replication of observed topology-breaking
#      rate under M3 (histogram vs observed)
#
# Data source: output/condlogit_diagnostics.json (built by
#   scripts/compute_condlogit_diagnostics.py), plus
#   output/manuscript_stats.json for cross-checks / labels.
#
# Output: output/figures/FigG_condlogit_diagnostics.png
# ============================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
})

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      "src/codon_topo/visualization/R"
    }
  }
)
source(file.path(script_dir, "theme_codon.R"))

# ── Locate repo root and load data ────────────────────────────
repo_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))
diag_path <- file.path(repo_root, "output", "condlogit_diagnostics.json")
stats_path <- file.path(repo_root, "output", "manuscript_stats.json")
out_path <- file.path(repo_root, "output", "figures", "FigG_condlogit_diagnostics.png")

if (!file.exists(diag_path)) {
  stop("Missing ", diag_path, ". Run scripts/compute_condlogit_diagnostics.py first.")
}
diag <- fromJSON(diag_path, simplifyDataFrame = FALSE)
stats <- fromJSON(stats_path, simplifyDataFrame = FALSE)
cl <- stats$condlogit

# ── Panel A: Delta_phys vs Delta_topo joint distribution ──────
scatter <- diag$phys_topo_scatter_sample
rho <- diag$phys_topo_rho
if (is.null(rho)) rho <- cl$phys_topo_rho
n_pool <- scatter$n_pool
n_sample <- scatter$n_sample

df_scatter <- data.frame(
  delta_phys = unlist(scatter$delta_phys),
  delta_topo = unlist(scatter$delta_topo)
)

# Jitter delta_topo (integer-valued) for visibility; delta_phys is continuous.
pA <- ggplot(df_scatter, aes(x = delta_phys, y = delta_topo)) +
  geom_jitter(width = 0, height = 0.15, alpha = 0.10, size = 0.5,
              colour = PAL_BLUE) +
  geom_smooth(method = "loess", se = FALSE, colour = PAL_RED,
              linewidth = 0.7, formula = y ~ x, span = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Spearman ~ rho == %.2f", rho),
           parse = TRUE, hjust = 1.05, vjust = 1.4, size = ANNOT_SIZE,
           colour = "grey20") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("italic(n) == '%s of %s'",
                           format(n_sample, big.mark = ","),
                           format(n_pool, big.mark = ",")),
           parse = TRUE, hjust = 1.05, vjust = 3.0, size = ANNOT_SIZE_S,
           colour = "grey40") +
  labs(
    title = "A. Feature joint distribution",
    subtitle = expression(Delta[phys] ~ "vs" ~ Delta[topo] ~ "across candidate moves"),
    x = expression(Delta[phys] ~ "(local Grantham cost change)"),
    y = expression(Delta[topo] ~ "(" * Q[6] ~ "component-count change)")
  ) +
  theme_codon_pub()

# ── Panel B: M3 + M4 coefficient forest (95% Wald CIs) ────────
model_ses <- diag$model_ses
# Only include models whose weight_labels != NULL (M3 and M4 have 2 and 3 params)
build_forest_df <- function(name) {
  block <- model_ses[[name]]
  if (is.null(block) || length(block$weights_normalized) == 0) return(NULL)
  data.frame(
    model = name,
    term = unlist(block$weight_labels),
    est_norm = unlist(block$weights_normalized),
    se_norm = unlist(block$weights_se_normalized),
    est_raw = unlist(block$weights_raw),
    se_raw = unlist(block$weights_se_raw),
    stringsAsFactors = FALSE
  )
}

forest_df <- do.call(rbind, lapply(c("M3_phys_topo", "M4_full"), build_forest_df))

# Present on the normalized (z-scored feature) scale — coefficients are
# on a common scale so magnitudes are directly comparable.
forest_df$lo_norm <- forest_df$est_norm - 1.96 * forest_df$se_norm
forest_df$hi_norm <- forest_df$est_norm + 1.96 * forest_df$se_norm

# Nicer term labels for display
term_labels <- c(
  delta_phys = "beta[phys]",
  delta_topo = "beta[topo]",
  delta_trna = "beta[tRNA]"
)
forest_df$term_lbl <- factor(
  unname(term_labels[forest_df$term]),
  levels = c("beta[phys]", "beta[topo]", "beta[tRNA]")
)
model_labels <- c(M3_phys_topo = "M3 (phys + topo)",
                  M4_full      = "M4 (phys + topo + tRNA)")
forest_df$model_lbl <- factor(unname(model_labels[forest_df$model]),
                              levels = c("M3 (phys + topo)",
                                         "M4 (phys + topo + tRNA)"))

pB <- ggplot(forest_df,
             aes(x = est_norm, y = term_lbl, colour = model_lbl)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo_norm, xmax = hi_norm),
                width = 0.18, linewidth = 0.6,
                orientation = "y",
                position = position_dodge(width = 0.5)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_colour_manual(values = c("M3 (phys + topo)" = PAL_BLUE,
                                 "M4 (phys + topo + tRNA)" = PAL_TEAL)) +
  labs(
    title = "B. Coefficient estimates",
    subtitle = "Standardized-feature scale (95% Wald CI)",
    x = "Coefficient (normalized)",
    y = NULL,
    colour = NULL
  ) +
  theme_codon_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9))

# ── Panel C: posterior-predictive histogram ───────────────────
pp <- diag$posterior_predictive
sim_rates <- unlist(pp$simulated_rates)
obs_rate <- pp$observed_topo_breaking_rate
sim_mean <- pp$simulated_mean_rate
sim_sd <- pp$simulated_std_rate
pp_p <- pp$posterior_predictive_p
n_sim <- pp$n_simulations

df_pp <- data.frame(rate = sim_rates)

# Compute y-limit that leaves headroom for the annotations above the histogram
sim_hist <- hist(sim_rates, breaks = seq(0, max(sim_rates) + 0.005, length.out = 41),
                 plot = FALSE)
y_max <- max(sim_hist$counts) * 1.25

pC <- ggplot(df_pp, aes(x = rate)) +
  geom_histogram(bins = 40, fill = PAL_BLUE_LT, colour = PAL_BLUE,
                 linewidth = 0.25) +
  geom_vline(xintercept = obs_rate, linetype = "solid", colour = PAL_RED,
             linewidth = 0.8) +
  geom_vline(xintercept = sim_mean, linetype = "dashed", colour = "grey30",
             linewidth = 0.5) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0))) +
  # Right-aligned annotation lines in the top-right corner.
  annotate("text", x = Inf, y = y_max * 0.97,
           label = sprintf("— observed = %.3f", obs_rate),
           hjust = 1.05, vjust = 1, size = ANNOT_SIZE, colour = PAL_RED) +
  annotate("text", x = Inf, y = y_max * 0.87,
           label = sprintf("- - sim mean = %.3f (SD %.3f)", sim_mean, sim_sd),
           hjust = 1.05, vjust = 1, size = ANNOT_SIZE_S, colour = "grey25") +
  annotate("text", x = Inf, y = y_max * 0.79,
           label = sprintf("pp p = %.2f", pp_p),
           hjust = 1.05, vjust = 1, size = ANNOT_SIZE_S, colour = "grey25",
           fontface = "italic") +
  labs(
    title = "C. Posterior-predictive check",
    subtitle = sprintf("M3 simulated topology-breaking rate (%s draws)",
                        format(n_sim, big.mark = ",")),
    x = "Simulated topology-breaking rate",
    y = "Frequency"
  ) +
  theme_codon_pub()

# ── Assemble ──────────────────────────────────────────────────
final <- (pA | pB | pC) +
  plot_layout(widths = c(1.1, 1.0, 1.1)) +
  plot_annotation(
    caption = "Companion to Fig 5 (AICc + rank + LRT). Panels here are non-overlapping diagnostics."
  ) &
  theme(plot.caption = element_text(size = 8, colour = "grey40",
                                     hjust = 0))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
ggsave(out_path, final, width = 13, height = 4.3, dpi = DPI, bg = "white")
cat("Saved:", out_path, "\n")
