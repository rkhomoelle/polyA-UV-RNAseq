## ============================================================
## 05_length_analysis.R
## Pulls transcript/CDS/UTR lengths from BioMart and merges with
## DESeq2 results, producing per-group (up/down/ns) CSVs for Prism.
##
## Run for:
##   NC vs CL:  up, down, and not-significant genes
##   NC vs INT: down and not-significant genes only
##     (no "up" group for NC vs INT, per the original analysis)
##
## Includes a length-validation check: flags any transcript where
## UTR5 + CDS + UTR3 deviates substantially from transcript_length,
## which can indicate BioMart double-counted a UTR segment rather than
## the transcript genuinely having a very long UTR.
## ============================================================

source(here::here("scripts", "00_setup.R"))

res_NC_CL  <- readRDS(RDS_NC_CL)
res_NC_INT <- readRDS(RDS_NC_INT)

## ── One BioMart pull covering genes from both comparisons ────
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

all_gene_ids <- union(rownames(res_NC_CL), rownames(res_NC_INT))

transcript_info <- getBM(
  attributes = c(
    "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
    "transcript_length", "cds_length",
    "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end"
  ),
  filters = "ensembl_gene_id",
  values  = all_gene_ids,
  mart    = mart
)

names(transcript_info)[names(transcript_info) == "external_gene_name"] <- "description"
names(transcript_info)[names(transcript_info) == "5_utr_start"] <- "utr5_start"
names(transcript_info)[names(transcript_info) == "5_utr_end"]   <- "utr5_end"
names(transcript_info)[names(transcript_info) == "3_utr_start"] <- "utr3_start"
names(transcript_info)[names(transcript_info) == "3_utr_end"]   <- "utr3_end"

## Sum UTR segments per transcript (BioMart returns one row per segment)
transcript_info_clean <- transcript_info %>%
  mutate(
    utr5_seg_len = ifelse(!is.na(utr5_start) & !is.na(utr5_end), abs(utr5_end - utr5_start) + 1, 0),
    utr3_seg_len = ifelse(!is.na(utr3_start) & !is.na(utr3_end), abs(utr3_end - utr3_start) + 1, 0)
  ) %>%
  group_by(ensembl_gene_id, description, ensembl_transcript_id, transcript_length, cds_length) %>%
  summarise(
    utr5_length = sum(utr5_seg_len, na.rm = TRUE),
    utr3_length = sum(utr3_seg_len, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("Rows after collapsing UTRs: %d unique transcripts\n", nrow(transcript_info_clean)))

## ── Merge one DESeq2 comparison with transcript lengths, split by group ──
run_length_analysis <- function(res, comparison_label, groups_to_export = c("up", "down", "ns")) {
  
  res_df <- as.data.frame(res) %>% rownames_to_column("ensembl_gene_id")
  
  merged <- res_df %>%
    inner_join(transcript_info_clean, by = "ensembl_gene_id") %>%
    dplyr::select(ensembl_gene_id, description, log2FoldChange, pvalue, padj,
                  ensembl_transcript_id, transcript_length, cds_length, utr5_length, utr3_length)
  
  cat(sprintf("[%s] DESeq2 genes: %d | after merge: %d rows | unique genes: %d | unique transcripts: %d\n",
              comparison_label, nrow(res_df), nrow(merged),
              n_distinct(merged$ensembl_gene_id), n_distinct(merged$ensembl_transcript_id)))
  
  ## Longest transcript per gene, requiring a CDS annotation
  merged_longest <- merged %>%
    filter(!is.na(cds_length), cds_length > 0) %>%
    group_by(ensembl_gene_id) %>%
    slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  cat(sprintf("[%s] Genes after requiring CDS annotation: %d (of %d)\n",
              comparison_label, n_distinct(merged_longest$ensembl_gene_id), n_distinct(merged$ensembl_gene_id)))
  
  ## Sanity check: UTR5 + CDS + UTR3 should roughly reconstruct
  ## transcript_length. Large deviations can indicate a BioMart
  ## double-counted UTR segment rather than a genuinely long UTR.
  validate_transcript_lengths(merged_longest, comparison_label, TABLES_DIR)
  
  degs <- split_degs(merged_longest, PADJ_THRESH, LFC_THRESH)
  names(degs) <- c("up", "down", "ns")   # split_degs() already returns this order
  
  select_cols <- c("ensembl_gene_id", "description", "group", "log2FoldChange", "padj",
                   "ensembl_transcript_id", "transcript_length", "utr5_length", "cds_length", "utr3_length")
  
  for (g in groups_to_export) {
    df <- degs[[g]]
    write.csv(df, file.path(TABLES_DIR, sprintf("%s_DEGs_%s_len.csv", comparison_label, g)), row.names = FALSE)
    write.csv(df %>% dplyr::select(all_of(select_cols)),
              file.path(TABLES_DIR, sprintf("prism_%s_%s.csv", comparison_label, g)), row.names = FALSE)
  }
  
  prism_all <- bind_rows(degs[groups_to_export]) %>% dplyr::select(all_of(select_cols))
  write.csv(prism_all, file.path(TABLES_DIR, sprintf("prism_%s_all.csv", comparison_label)), row.names = FALSE)
  
  cat(sprintf("[%s] Median lengths per group:\n", comparison_label))
  print(
    prism_all %>%
      group_by(group) %>%
      summarise(n = n(),
                transcript_bp = median(transcript_length, na.rm = TRUE),
                utr5_bp = median(utr5_length, na.rm = TRUE),
                cds_bp = median(cds_length, na.rm = TRUE),
                utr3_bp = median(utr3_length, na.rm = TRUE),
                .groups = "drop")
  )
  
  invisible(prism_all)
}

## NC vs CL: up, down, and not-significant
run_length_analysis(res_NC_CL, "NC_vs_CL", groups_to_export = c("up", "down", "ns"))

## NC vs INT: down and not-significant only (no "up" group used in this comparison)
run_length_analysis(res_NC_INT, "NC_vs_INT", groups_to_export = c("down", "ns"))
