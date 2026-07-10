## ============================================================
## 05_GO_enrichment_MF.R
## GO Molecular Function enrichment on protein-coding DEGs.
## ============================================================

source(here::here("scripts", "00_setup.R"))

res_NC_CL  <- readRDS(RDS_NC_CL)
res_NC_INT <- readRDS(RDS_NC_INT)

res_NC_CL_annot  <- annotate_results(res_NC_CL)
res_NC_INT_annot <- annotate_results(res_NC_INT)

summarise_degs(res_NC_CL_annot,  "NC vs CL (all genes)")
summarise_degs(res_NC_INT_annot, "NC vs INT (all genes)")

## ── Filter to protein-coding ───────────────────────────────────
res_NC_CL_pc  <- res_NC_CL_annot  %>% filter(GENETYPE == "protein-coding", !is.na(SYMBOL))
res_NC_INT_pc <- res_NC_INT_annot %>% filter(GENETYPE == "protein-coding", !is.na(SYMBOL))

summarise_degs(res_NC_CL_pc,  "NC vs CL (protein-coding)")
summarise_degs(res_NC_INT_pc, "NC vs INT (protein-coding)")

## ── DEG subsets + shared universe ──────────────────────────────
NC_CL_up    <- subset(res_NC_CL_pc,  padj < PADJ_THRESH & log2FoldChange >=  LFC_THRESH)$SYMBOL
NC_CL_down  <- subset(res_NC_CL_pc,  padj < PADJ_THRESH & log2FoldChange <= -LFC_THRESH)$SYMBOL
NC_INT_down <- subset(res_NC_INT_pc, padj < PADJ_THRESH & log2FoldChange <= -LFC_THRESH)$SYMBOL

## Universe: protein-coding genes detected in both datasets
universe <- intersect(res_NC_CL_pc$SYMBOL, res_NC_INT_pc$SYMBOL)

NC_CL_up_clean    <- intersect(NC_CL_up,    universe)
NC_CL_down_clean  <- intersect(NC_CL_down,  universe)
NC_INT_down_clean <- intersect(NC_INT_down, universe)

## ── Fisher's Exact Test: NC-CL up/down vs NC-INT down ─────────
results_fisher <- rbind(
  CL_up_vs_INT_down   = compare_sets(NC_CL_up_clean,   NC_INT_down_clean, universe),
  CL_down_vs_INT_down = compare_sets(NC_CL_down_clean,  NC_INT_down_clean, universe)
)
results_fisher$padj <- p.adjust(results_fisher$p_value, method = "BH")
print(results_fisher)
write.csv(results_fisher, file.path(TABLES_DIR, "fisher_tests_NC_CL_vs_NC_INT.csv"), row.names = TRUE)

## ── Overlap genes: NC-CL up & NC-INT down ─────────────────────
overlap_genes <- intersect(NC_CL_up_clean, NC_INT_down_clean)
cat(sprintf("Overlap genes (NC-CL up / NC-INT down): %d\n", length(overlap_genes)))
write.table(overlap_genes, file.path(TABLES_DIR, "overlap_NC_CL_up_NC_INT_down.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## ── GO Molecular Function enrichment ──────────────────────────
universe_entrez     <- symbols_to_entrez(universe)
NC_CL_up_entrez     <- symbols_to_entrez(NC_CL_up_clean)
NC_CL_down_entrez   <- symbols_to_entrez(NC_CL_down_clean)
NC_INT_down_entrez  <- symbols_to_entrez(NC_INT_down_clean)

cat(sprintf("NC-CL up:    %d / %d genes with ENTREZID\n", length(NC_CL_up_entrez),    length(NC_CL_up_clean)))
cat(sprintf("NC-CL down:  %d / %d genes with ENTREZID\n", length(NC_CL_down_entrez),  length(NC_CL_down_clean)))
cat(sprintf("NC-INT down: %d / %d genes with ENTREZID\n", length(NC_INT_down_entrez), length(NC_INT_down_clean)))

go_NC_CL_up    <- run_go_mf(NC_CL_up_entrez,    universe_entrez)
go_NC_CL_down  <- run_go_mf(NC_CL_down_entrez,  universe_entrez)
go_NC_INT_down <- run_go_mf(NC_INT_down_entrez, universe_entrez)

save_go_results(go_NC_CL_up,    "NC_CL_up",    FIG_DIR, TABLES_DIR)
save_go_results(go_NC_CL_down,  "NC_CL_down",  FIG_DIR, TABLES_DIR)
save_go_results(go_NC_INT_down, "NC_INT_down", FIG_DIR, TABLES_DIR)

## Enriched vs depleted, side by side (NC-CL only, since NC-INT has no "up" GO run)
p_NC_CL_compare <- plot_grid(
  make_dotplot(go_NC_CL_up,   "NC-CL -- Enriched"),
  make_dotplot(go_NC_CL_down, "NC-CL -- Depleted"),
  ncol = 2, align = "h"
)
ggsave(file.path(FIG_DIR, "GO_MF_NC_CL_enriched_vs_depleted.pdf"), p_NC_CL_compare, width = 16, height = 8)
ggsave(file.path(FIG_DIR, "GO_MF_NC_CL_enriched_vs_depleted.png"), p_NC_CL_compare, width = 16, height = 8, dpi = 300)
print(p_NC_CL_compare)

## Shared GO terms between NC-CL up and NC-INT down
df_NC_CL_up    <- as.data.frame(go_NC_CL_up)
df_NC_INT_down <- as.data.frame(go_NC_INT_down)

if (nrow(df_NC_CL_up) > 0 && nrow(df_NC_INT_down) > 0) {
  shared_terms <- intersect(df_NC_CL_up$Description, df_NC_INT_down$Description)
  cat(sprintf("\nShared GO MF terms (NC-CL up & NC-INT down): %d\n", length(shared_terms)))
  print(shared_terms)
  write.csv(data.frame(shared_GO_term = shared_terms),
            file.path(TABLES_DIR, "GO_MF_shared_NC_CL_up_NC_INT_down.csv"), row.names = FALSE)
}

## ── Extract highlight genes for the quadrant plot (06) ────────
## Update the grep pattern to match whichever GO terms you want to
## highlight in the quadrant plot.
highlight_genes <- df_NC_INT_down %>%
  filter(grepl("ribosome|NADH dehydrogenase|electron transfer", Description)) %>%
  pull(geneID) %>%
  paste(collapse = "/") %>%
  strsplit("/") %>%
  unlist() %>%
  unique() %>%
  intersect(overlap_genes)

cat(sprintf("\nOverlap genes in GO terms of interest: %d\n", length(highlight_genes)))
print(highlight_genes)

saveRDS(highlight_genes, file.path(TABLES_DIR, "go_highlight_genes.rds"))
saveRDS(list(NC_CL_pc = res_NC_CL_pc, NC_INT_pc = res_NC_INT_pc),
        file.path(TABLES_DIR, "protein_coding_annotated.rds"))

## ── Save gene sets belonging to significant GO MF terms ───────
## Used by 03_volcano_MA_plots.R to label these genes on the
## volcano plots. NC-INT has no "up" entry since no GO MF run
## was done on NC-INT up genes.
go_significant_genes <- list(
  NC_CL_up    = extract_go_genes(go_NC_CL_up),
  NC_CL_down  = extract_go_genes(go_NC_CL_down),
  NC_INT_down = extract_go_genes(go_NC_INT_down)
)
saveRDS(go_significant_genes, file.path(TABLES_DIR, "go_significant_genes.rds"))

## Also save the full GO objects themselves (not just extracted gene
## lists) so downstream scripts can pull genes from a SPECIFIC term by
## name/pattern.
saveRDS(
  list(go_NC_CL_up = go_NC_CL_up, go_NC_CL_down = go_NC_CL_down, go_NC_INT_down = go_NC_INT_down),
  file.path(TABLES_DIR, "go_objects.rds")
)

## Save the ENTREZID gene sets + universe too, so downstream scripts can
## re-query a SPECIFIC GO term by name (see get_genes_for_go_term() in
## utils.R) without needing to re-derive them from the DESeq2 results.
saveRDS(
  list(
    NC_CL_up_entrez    = NC_CL_up_entrez,
    NC_CL_down_entrez  = NC_CL_down_entrez,
    NC_INT_down_entrez = NC_INT_down_entrez,
    universe_entrez    = universe_entrez
  ),
  file.path(TABLES_DIR, "go_entrez_sets.rds")
)
