## ============================================================
## 07_correlation_plots.R
## Replicate-vs-replicate correlation scatterplots (log10 TPM + 1),
## generated for every pairwise combination of biological
## replicates within each of the three conditions (NC, CL, INT).
##
## CONFIG: this assumes a single TPM table containing columns
## A1-A3 (NC), B1-B3 (CL), and C1-C3 (INT). If your NC-CL and
## NC-INT TPM tables are separate files, read/merge them below
## before building `tpm_all`.
## ============================================================

source(here::here("scripts", "00_setup.R"))

## Load and, if needed, merge TPM tables so all conditions are present
tpm_NC_CL  <- read.table(TPM_NC_CL_FILE,  header = TRUE, row.names = 1)
tpm_NC_INT <- read.table(TPM_NC_INT_FILE, header = TRUE, row.names = 1)

## Combine on shared gene rownames; adjust join logic if your files
## don't share identical row ordering/IDs.
common_genes <- intersect(rownames(tpm_NC_CL), rownames(tpm_NC_INT))
tpm_all <- cbind(
  tpm_NC_CL[common_genes, c("A1", "A2", "A3", "B1", "B2", "B3")],
  tpm_NC_INT[common_genes, c("C1", "C2", "C3")]
)

pseudo_count <- 1
log10_tpm <- log10(tpm_all + pseudo_count)

## ── Reusable replicate-correlation plot ───────────────────────
plot_replicate_correlation <- function(col_x, col_y, condition_label) {
  
  df <- data.frame(x = log10_tpm[[col_x]], y = log10_tpm[[col_y]])
  pearson_r <- cor(df$x, df$y, method = "pearson")
  cat(sprintf("[%s] %s vs %s: Pearson r = %.4f\n", condition_label, col_x, col_y, pearson_r))
  
  axis_min <- floor(min(c(df$x, df$y)))
  axis_max <- ceiling(max(c(df$x, df$y)))
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(alpha = 0.7, size = 2.5, color = "#2C7BB6") +
    geom_smooth(method = "lm", se = TRUE, color = "#D7191C", fill = "#D7191C", alpha = 0.15, linewidth = 0.8) +
    annotate("text", x = axis_min + 0.1 * (axis_max - axis_min), y = axis_max - 0.05 * (axis_max - axis_min),
             label = sprintf("r = %.4f", pearson_r), hjust = 0, size = 4.5, fontface = "bold", color = "black") +
    coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
    labs(
      title    = sprintf("%s: %s vs %s", condition_label, col_x, col_y),
      subtitle = sprintf("log10(TPM + %g) | Pearson r = %.4f", pseudo_count, pearson_r),
      x = sprintf("%s log10(TPM + 1)", col_x),
      y = sprintf("%s log10(TPM + 1)", col_y)
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "grey40", size = 11),
          panel.grid.minor = element_blank())
  
  ## NOTE: print(p) intentionally omitted here. Calling print() on many
  ## ggplot objects back-to-back in a loop (9 iterations here: 3
  ## replicate pairs x 3 conditions) can corrupt RStudio's graphics
  ## device viewport stack, producing an error like:
  ##   "no applicable method for 'depth' applied to an object of class NULL"
  ## ggsave() renders independently of the interactive plot pane, so
  ## every PDF still gets saved correctly without print(). Open the
  ## saved files in results/figures/ to view them instead of relying on
  ## the interactive pop-up.
  out_file <- file.path(FIG_DIR, sprintf("correlation_%s_%s_vs_%s.pdf", condition_label, col_x, col_y))
  ggsave(out_file, plot = p, width = 6, height = 6, dpi = 300)
  cat(sprintf("Plot saved to: %s\n", out_file))
  
  invisible(pearson_r)
}

## ── Generate every pairwise replicate combination per condition ──
replicate_cols <- list(
  NC  = c("A1", "A2", "A3"),
  CL  = c("B1", "B2", "B3"),
  INT = c("C1", "C2", "C3")
)

for (condition_label in names(replicate_cols)) {
  cols <- replicate_cols[[condition_label]]
  pairs <- combn(cols, 2, simplify = FALSE)
  for (pair in pairs) {
    plot_replicate_correlation(pair[1], pair[2], condition_label)
  }
}