## ============================================================
## utils.R
## Shared helper functions used by multiple scripts.
## Sourced automatically by 00_setup.R -- do not source directly.
## ============================================================

## ── Run DESeq2 from a raw counts file ────────────────────────
## conditions[1] is used as the reference level.
run_deseq2 <- function(count_file, conditions) {
  counts <- read.table(count_file, header = TRUE, row.names = 1)
  
  # Collapse Ensembl version suffixes (e.g. ENSG00000123.4 -> ENSG00000123)
  gene_ids_clean   <- sub("\\..*", "", rownames(counts))
  counts_collapsed <- rowsum(counts, group = gene_ids_clean)
  
  sample_info <- data.frame(
    row.names = colnames(counts_collapsed),
    condition = factor(conditions)
  )
  stopifnot(all(colnames(counts_collapsed) == rownames(sample_info)))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_collapsed,
    colData   = sample_info,
    design    = ~ condition
  )
  dds$condition <- relevel(dds$condition, ref = conditions[1])
  dds <- DESeq(dds)
  results(dds)
}

## ── Annotate a DESeq2 results object with SYMBOL + GENETYPE ──
annotate_results <- function(res) {
  gene_symbols <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = rownames(res),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENETYPE")
  )
  
  res_df <- as.data.frame(res)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- res_df %>% filter(!is.na(log2FoldChange))
  
  merge(
    res_df,
    gene_symbols,
    by.x = "ensembl_gene_id",
    by.y = "ENSEMBL",
    all.x = TRUE
  )
}

## ── Print up/down DEG counts ──────────────────────────────────
summarise_degs <- function(res_annot, dataset_name,
                           padj_thresh = 0.05, lfc_thresh = 1) {
  up   <- sum(res_annot$padj < padj_thresh & res_annot$log2FoldChange >=  lfc_thresh, na.rm = TRUE)
  down <- sum(res_annot$padj < padj_thresh & res_annot$log2FoldChange <= -lfc_thresh, na.rm = TRUE)
  cat(sprintf(
    "%s -- Significant DEGs (padj < %g, |LFC| >= %g):\n  Up:    %d\n  Down:  %d\n  Total: %d\n\n",
    dataset_name, padj_thresh, lfc_thresh, up, down, up + down
  ))
}

## ── Split an annotated results data frame into up/down/ns ────
split_degs <- function(res_annot, padj_thresh = 0.05, lfc_thresh = 1) {
  up <- res_annot %>%
    filter(!is.na(padj), padj < padj_thresh, log2FoldChange > lfc_thresh) %>%
    mutate(group = "upregulated")
  
  dn <- res_annot %>%
    filter(!is.na(padj), padj < padj_thresh, log2FoldChange < -lfc_thresh) %>%
    mutate(group = "downregulated")
  
  ns <- res_annot %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    filter(padj >= padj_thresh | abs(log2FoldChange) <= lfc_thresh) %>%
    mutate(group = "not_significant")
  
  cat(sprintf("Upregulated:     %d genes\n", nrow(up)))
  cat(sprintf("Downregulated:   %d genes\n", nrow(dn)))
  cat(sprintf("Not significant: %d genes\n\n", nrow(ns)))
  
  list(up = up, down = dn, ns = ns)
}

## ── Fisher's Exact Test between two gene sets against a universe ──
compare_sets <- function(A, B, universe) {
  A <- intersect(A, universe)
  B <- intersect(B, universe)
  
  a <- length(intersect(A, B))
  b <- length(setdiff(A, B))
  c <- length(setdiff(B, A))
  d <- length(setdiff(universe, union(A, B)))
  
  cont_matrix <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  test <- fisher.test(cont_matrix)
  
  data.frame(
    overlap    = a,
    odds_ratio = as.numeric(test$estimate),
    p_value    = test$p.value
  )
}

## ── Convert SYMBOL -> ENTREZID ────────────────────────────────
symbols_to_entrez <- function(symbols) {
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = symbols,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
    pull(ENTREZID)
}

## ── Run enrichGO Molecular Function + simplify ────────────────
run_go_mf <- function(gene_entrez, universe_entrez) {
  go <- enrichGO(
    gene          = gene_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

## ── Save GO results: CSV + dot plot + cnet plot ───────────────
save_go_results <- function(go_obj, label, fig_dir, table_dir) {
  df <- as.data.frame(go_obj)
  cat(sprintf("\n%s -- significant MF GO terms: %d\n", label, nrow(df)))
  
  write.csv(df, file.path(table_dir, sprintf("GO_MF_%s_results.csv", label)), row.names = FALSE)
  
  if (nrow(df) == 0) {
    message(sprintf("No significant GO terms for %s -- skipping plot", label))
    return(invisible(NULL))
  }
  
  p_dot <- dotplot(go_obj, showCategory = 20) +
    labs(title = sprintf("GO Molecular Function\n%s", gsub("_", " ", label))) +
    theme_classic(base_size = 16) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title.x = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 13)
    )
  ggsave(file.path(fig_dir, sprintf("GO_MF_%s_dotplot.pdf", label)), plot = p_dot, width = 10, height = 8)
  ggsave(file.path(fig_dir, sprintf("GO_MF_%s_dotplot.png", label)), plot = p_dot, width = 10, height = 8, dpi = 300)
  
  invisible(go_obj)
}

## ── Dot plot panel that degrades gracefully if empty ──────────
make_dotplot <- function(go_obj, title) {
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant terms", size = 5) +
      labs(title = title) +
      theme_void() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
  } else {
    dotplot(go_obj, showCategory = 15) +
      labs(title = title) +
      theme_classic(base_size = 15) +
      theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 15),
        axis.text.x  = element_text(size = 13),
        axis.text.y  = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)
      )
  }
}

## ── Extract the set of genes belonging to any significant term ──
## in a GO enrichment result (assumes readable = TRUE, so geneID
## is already SYMBOL rather than ENTREZID).
extract_go_genes <- function(go_obj) {
  if (is.null(go_obj)) return(character(0))
  df <- as.data.frame(go_obj)
  if (nrow(df) == 0) return(character(0))
  df$geneID %>%
    paste(collapse = "/") %>%
    strsplit("/") %>%
    unlist() %>%
    unique()
}

## ── Extract genes belonging to GO terms matching a description ──
## pattern (e.g. "transmembrane transporter"), rather than ALL
## significant terms. Case-insensitive grep against Description.
extract_go_genes_for_terms <- function(go_obj, pattern) {
  if (is.null(go_obj)) return(character(0))
  df <- as.data.frame(go_obj)
  if (nrow(df) == 0) return(character(0))
  df %>%
    filter(grepl(pattern, Description, ignore.case = TRUE)) %>%
    pull(geneID) %>%
    paste(collapse = "/") %>%
    strsplit("/") %>%
    unlist() %>%
    unique()
}

## ── Get genes for a specific/named GO term, bypassing simplify() ──
## run_go_mf()'s simplify() step can collapse a term (e.g. "transmembrane
## transporter activity") into a more specific child term, dropping it
## from the final results even if it was significant. This re-runs
## enrichGO fresh (same params, no simplify) and filters directly by
## Description, so a term you specifically want won't be lost to that
## collapsing.
get_genes_for_go_term <- function(gene_entrez, universe_entrez, pattern) {
  if (length(gene_entrez) == 0) return(character(0))
  go <- enrichGO(
    gene          = gene_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  extract_go_genes_for_terms(go, pattern)
}

## ── Validate that UTR5 + CDS + UTR3 roughly reconstructs the total ──
## transcript length. Flags rows where the parts deviate from
## transcript_length by more than `tolerance` (fractional, default 10%),
## which can indicate a BioMart double-counted UTR segment (e.g. a
## duplicated exon row) rather than a genuinely long UTR. Saves flagged
## rows to a CSV and returns them invisibly.
validate_transcript_lengths <- function(df, comparison_label, table_dir, tolerance = 0.10) {
  checked <- df %>%
    mutate(
      sum_parts     = utr5_length + cds_length + utr3_length,
      deviation     = sum_parts - transcript_length,
      pct_deviation = deviation / transcript_length
    )
  
  flagged <- checked %>%
    filter(abs(pct_deviation) > tolerance) %>%
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, description,
                  transcript_length, cds_length, utr5_length, utr3_length,
                  sum_parts, deviation, pct_deviation) %>%
    arrange(desc(abs(deviation)))
  
  cat(sprintf(
    "[%s] Length validation -- %d of %d transcripts have UTR5+CDS+UTR3 deviating >%.0f%% from transcript_length\n",
    comparison_label, nrow(flagged), nrow(checked), tolerance * 100
  ))
  
  if (nrow(flagged) > 0) {
    out_file <- file.path(table_dir, sprintf("length_validation_flagged_%s.csv", comparison_label))
    write.csv(flagged, out_file, row.names = FALSE)
    cat(sprintf("[%s] Flagged transcripts saved to %s -- top 5 by |deviation|:\n", comparison_label, out_file))
    print(head(flagged, 5))
  }
  
  invisible(flagged)
}

## ── Ensure the shared GTF is available locally, download once ─
ensure_gtf <- function(gtf_file, gtf_url) {
  if (!file.exists(gtf_file)) {
    message("GTF not found locally -- downloading (one-time)...")
    download.file(url = gtf_url, destfile = gtf_file)
  }
  gtf_file
}