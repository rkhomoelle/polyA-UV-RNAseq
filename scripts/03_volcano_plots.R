## ============================================================
## 03_volcano_plots.R
## Volcano plot + DEG up/down/ns splitting for each DESeq2
## comparison. Reads res_NC_CL / res_NC_INT from the .rds files
## produced by 01_deseq2.R (no DESeq2 re-run here).
##
## Volcano plots label genes belonging to significant GO MF
## terms from the up- and down-regulated GO MF analyses (from
## 02_GO_enrichment_MF.R). This
## script must run AFTER 02_GO_enrichment_MF.R.
## ============================================================

source(here::here("scripts", "00_setup.R"))

res_NC_CL  <- readRDS(RDS_NC_CL)
res_NC_INT <- readRDS(RDS_NC_INT)

## Gene sets belonging to significant GO MF terms, produced by
## 02_GO_enrichment_MF.R (run that script first if this errors).
go_significant_genes <- readRDS(file.path(TABLES_DIR, "go_significant_genes.rds"))
go_entrez_sets       <- readRDS(file.path(TABLES_DIR, "go_entrez_sets.rds"))

## Genes from any GO MF term matching this pattern are always included
## as volcano labels for NC-CL, in addition to the general significant-
## term set above -- this re-queries enrichGO directly (bypassing the
## simplify() step in 02) so a specific term you want won't get dropped
## by that semantic-similarity collapsing. Adjust the pattern if you
## want a different/additional term labeled.
NC_CL_TERM_PATTERN <- "transmembrane transporter"

transporter_genes_NC_CL <- union(
  get_genes_for_go_term(go_entrez_sets$NC_CL_down_entrez, go_entrez_sets$universe_entrez, NC_CL_TERM_PATTERN),
  get_genes_for_go_term(go_entrez_sets$NC_CL_up_entrez,   go_entrez_sets$universe_entrez, NC_CL_TERM_PATTERN)
)
cat(sprintf("Genes matching '%s' (NC-CL): %d\n", NC_CL_TERM_PATTERN, length(transporter_genes_NC_CL)))

## ── Reusable volcano + DEG-split routine ──────────────────────
## Labels only genes belonging to significant GO MF terms (go_genes),
## excluding any manually-highlighted genes (those get their own boxed
## label below instead).
make_volcano_plot <- function(res, comparison_label, highlight_genes = c("GAPDH", "ACTB"),
                               go_genes = character(0)) {

  res_annot <- annotate_results(res)

  ## Genes belonging to significant GO MF terms, restricted to genes
  ## actually present in this comparison's results and excluding
  ## the manually-highlighted genes (those get their own boxed labels below)
  all_labels <- setdiff(intersect(go_genes, res_annot$SYMBOL), highlight_genes)

  cat(sprintf("[%s] Volcano labels -- GO-term genes: %d\n", comparison_label, length(all_labels)))

  keyvals <- ifelse(
    res_annot$log2FoldChange < -LFC_THRESH & res_annot$padj < PADJ_THRESH, "#377EB8",
    ifelse(res_annot$log2FoldChange > LFC_THRESH & res_annot$padj < PADJ_THRESH, "#E41A1C", "grey")
  )
  keyvals[is.na(keyvals)] <- "grey"
  names(keyvals)[keyvals == "#E41A1C"] <- "Up"
  names(keyvals)[keyvals == "#377EB8"] <- "Down"
  names(keyvals)[keyvals == "grey"]    <- "NS"

  p_volcano <- EnhancedVolcano(res_annot,
    lab            = res_annot$SYMBOL,
    selectLab      = all_labels,
    x              = "log2FoldChange",
    y              = "padj",
    title          = comparison_label,
    axisLabSize    = 12,
    titleLabSize   = 16,
    subtitle       = NULL,
    pCutoff        = PADJ_THRESH,
    FCcutoff       = LFC_THRESH,
    pointSize      = 2.0,
    labSize        = 3.0,
    labCol         = "black",
    colCustom      = keyvals,
    colAlpha       = 0.5,
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    colConnectors  = "grey40",
    legendLabels   = NULL,
    legendPosition = "none"
  )

  highlight_df <- res_annot[res_annot$SYMBOL %in% highlight_genes, ]
  p_volcano_labeled <- p_volcano +
    geom_label_repel(
      data = highlight_df,
      aes(x = log2FoldChange, y = -log10(padj), label = SYMBOL),
      box.padding = 0.5, point.padding = 0.3,
      segment.color = "black", segment.size = 0.5,
      fill = "white", color = "black", fontface = "bold",
      size = 3, max.overlaps = Inf, min.segment.length = 0
    )
  print(p_volcano_labeled)
  ggsave(file.path(FIG_DIR, sprintf("volcano_%s.pdf", comparison_label)), p_volcano_labeled, width = 5, height = 5)

  ## DEG up/down/ns split
  degs <- split_degs(res_annot, PADJ_THRESH, LFC_THRESH)
  write.csv(degs$up,   file.path(TABLES_DIR, sprintf("DEGs_%s_up.csv", comparison_label)),   row.names = FALSE)
  write.csv(degs$down, file.path(TABLES_DIR, sprintf("DEGs_%s_dn.csv", comparison_label)),   row.names = FALSE)
  write.csv(degs$ns,   file.path(TABLES_DIR, sprintf("DEGs_%s_ns.csv", comparison_label)),   row.names = FALSE)

  list(res_annot = res_annot, degs = degs)
}

## ── Run for both comparisons ──────────────────────────────────
## NC-CL: label genes from significant GO terms in either direction,
## plus genes specifically from NC_CL_TERM_PATTERN above
## NC-INT: only "down" GO MF was run (no "up" GO run), so only that set is available
out_NC_CL <- make_volcano_plot(
  res_NC_CL, "NC_vs_CL",
  go_genes = union(
    union(go_significant_genes$NC_CL_up, go_significant_genes$NC_CL_down),
    transporter_genes_NC_CL
  )
)

out_NC_INT <- make_volcano_plot(
  res_NC_INT, "NC_vs_INT",
  go_genes = go_significant_genes$NC_INT_down
)
