## ============================================================
## 06_quadrant_plot.R
## Quadrant plot comparing log2FC in NC-CL vs NC-INT, with GO
## highlight genes overlaid. Reads intermediate objects saved by
## 02_GO_enrichment_MF.R.
##
## Includes a Fisher's exact test for the CL up / INT down quadrant
## overlap (bottom-right of the plot) against the shared
## protein-coding universe, testing whether that overlap is larger than
## expected by chance. 
## ============================================================

source(here::here("scripts", "00_setup.R"))

pc_annot        <- readRDS(file.path(TABLES_DIR, "protein_coding_annotated.rds"))
res_NC_CL_pc    <- pc_annot$NC_CL_pc
res_NC_INT_pc   <- pc_annot$NC_INT_pc
highlight_genes <- readRDS(file.path(TABLES_DIR, "go_highlight_genes.rds"))

plot_df <- inner_join(
  res_NC_CL_pc %>%
    dplyr::select(SYMBOL, log2FoldChange, padj) %>%
    rename(LFC_CL = log2FoldChange, padj_CL = padj) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE),
  res_NC_INT_pc %>%
    dplyr::select(SYMBOL, log2FoldChange, padj) %>%
    rename(LFC_INT = log2FoldChange, padj_INT = padj) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE),
  by = "SYMBOL"
) %>%
  mutate(
    group = case_when(
      padj_CL < PADJ_THRESH & LFC_CL >=  LFC_THRESH & padj_INT < PADJ_THRESH & LFC_INT >=  LFC_THRESH ~ "Up in both",
      padj_CL < PADJ_THRESH & LFC_CL <= -LFC_THRESH & padj_INT < PADJ_THRESH & LFC_INT <= -LFC_THRESH ~ "Down in both",
      padj_CL < PADJ_THRESH & LFC_CL >=  LFC_THRESH & padj_INT < PADJ_THRESH & LFC_INT <= -LFC_THRESH ~ "CL up / INT down",
      padj_CL < PADJ_THRESH & LFC_CL <= -LFC_THRESH & padj_INT < PADJ_THRESH & LFC_INT >=  LFC_THRESH ~ "CL down / INT up",
      TRUE ~ "NS"
    ),
    group = factor(group, levels = c("CL up / INT down", "Up in both", "Down in both", "CL down / INT up", "NS"))
  )

group_counts <- plot_df %>% filter(group != "NS") %>% count(group)

## ── Statistical test: is the CL up / INT down quadrant's overlap ──
## bigger than chance? (bottom-right quadrant of the plot)
NC_CL_up    <- subset(res_NC_CL_pc,  padj < PADJ_THRESH & log2FoldChange >=  LFC_THRESH)$SYMBOL
NC_INT_down <- subset(res_NC_INT_pc, padj < PADJ_THRESH & log2FoldChange <= -LFC_THRESH)$SYMBOL

quadrant_universe <- intersect(res_NC_CL_pc$SYMBOL, res_NC_INT_pc$SYMBOL)

quadrant_fisher <- rbind(
  `CL up / INT down` = compare_sets(NC_CL_up, NC_INT_down, quadrant_universe)
)
quadrant_fisher$padj <- p.adjust(quadrant_fisher$p_value, method = "BH")

cat("\nFisher's exact test -- CL up / INT down quadrant overlap vs chance:\n")
print(quadrant_fisher)
write.csv(quadrant_fisher, file.path(TABLES_DIR, "quadrant_fisher_tests.csv"), row.names = TRUE)

## Pull the quadrant's p-value for the plot annotation below 
get_p <- function(grp) {
  p <- quadrant_fisher[grp, "padj"]
  if (is.null(p) || length(p) == 0 || is.na(p)) "NA" else formatC(p, format = "e", digits = 1)
}

label_genes <- plot_df %>%
  filter(group == "CL up / INT down") %>%
  mutate(mean_abs_lfc = (abs(LFC_CL) + abs(LFC_INT)) / 2) %>%
  slice_max(mean_abs_lfc, n = 20)

highlight_df <- plot_df %>% filter(SYMBOL %in% highlight_genes)

get_n <- function(grp) {
  n <- group_counts %>% filter(group == grp) %>% pull(n)
  if (length(n) == 0) 0 else n
}

group_colours <- c("CL up / INT down" = "#E63946", "Up in both" = "#F4A261",
                   "Down in both" = "#457B9D", "CL down / INT up" = "#2A9D8F", "NS" = "grey80")
group_sizes   <- c("CL up / INT down" = 2.5, "Up in both" = 2.0,
                   "Down in both" = 2.0, "CL down / INT up" = 2.0, "NS" = 0.8)
group_alpha   <- c("CL up / INT down" = 1.0, "Up in both" = 0.9,
                   "Down in both" = 0.9, "CL down / INT up" = 0.9, "NS" = 0.3)

p_quad <- ggplot(plot_df, aes(x = LFC_CL, y = LFC_INT, colour = group, size = group, alpha = group)) +
  geom_hline(yintercept = 0,  colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0,  colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept =  LFC_THRESH, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = -LFC_THRESH, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept =  LFC_THRESH, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = -LFC_THRESH, colour = "grey40", linewidth = 0.3, linetype = "dashed") +
  geom_point(data = filter(plot_df, group == "NS"), shape = 16) +
  geom_point(data = filter(plot_df, group != "NS"), shape = 16) +
  geom_point(data = highlight_df, aes(x = LFC_CL, y = LFC_INT), shape = 21, size = 4,
             colour = "black", fill = "purple", alpha = 1, inherit.aes = FALSE) +
  geom_text_repel(data = label_genes, aes(label = SYMBOL), size = 4.2, colour = "#E63946",
                  fontface = "italic", box.padding = 0.6, point.padding = 0.4,
                  segment.colour = "grey50", segment.size = 0.3, max.overlaps = Inf,
                  force = 2, min.segment.length = 0) +
  geom_text_repel(data = highlight_df, aes(x = LFC_CL, y = LFC_INT, label = SYMBOL), size = 4.2,
                  colour = "purple", fontface = "italic", box.padding = 0.6, point.padding = 0.4,
                  segment.colour = "grey50", segment.size = 0.3, max.overlaps = Inf,
                  force = 2, min.segment.length = 0, inherit.aes = FALSE) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Up in both\nn = %d", get_n("Up in both")),
           hjust = 1.05, vjust = 1.5, size = 4.5, colour = "#F4A261", fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf, label = sprintf("Down in both\nn = %d", get_n("Down in both")),
           hjust = -0.05, vjust = -0.5, size = 4.5, colour = "#457B9D", fontface = "bold") +
  annotate("text", x = Inf, y = -Inf, label = sprintf("CL up / INT down\nn = %d, p.adj = %s", get_n("CL up / INT down"), get_p("CL up / INT down")),
           hjust = 1.05, vjust = -0.5, size = 4.5, colour = "#E63946", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = sprintf("CL down / INT up\nn = %d", get_n("CL down / INT up")),
           hjust = -0.05, vjust = 1.5, size = 4.5, colour = "#2A9D8F", fontface = "bold") +
  scale_colour_manual(values = group_colours) +
  scale_size_manual(values = group_sizes) +
  scale_alpha_manual(values = group_alpha) +
  labs(x = expression("NC-CL  log"[2]*"(fold change)"),
       y = expression("NC-INT  log"[2]*"(fold change)"),
       title = "Opposing regulation between NC-CL and NC-INT", colour = NULL) +
  guides(size = "none", alpha = "none") +
  theme_classic(base_size = 18) +
  theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 19),
        axis.title   = element_text(size = 18),
        axis.text    = element_text(size = 15),
        legend.position = "top",
        legend.text  = element_text(size = 14),
        axis.line    = element_line(colour = "grey40"))

print(p_quad)
ggsave(file.path(FIG_DIR, "quadrant_plot_NC_CL_vs_NC_INT.pdf"), p_quad, width = 10, height = 11)
ggsave(file.path(FIG_DIR, "quadrant_plot_NC_CL_vs_NC_INT.png"), p_quad, width = 10, height = 11, dpi = 300)