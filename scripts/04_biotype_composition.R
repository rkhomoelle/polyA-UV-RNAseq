## ============================================================
## 04_biotype_composition.R
## RNA biotype composition from TPM values (all genes), per
## condition/replicate. Run for both NC-CL and NC-INT.
## ============================================================

source(here::here("scripts", "00_setup.R"))

ensure_gtf(GTF_FILE, GTF_URL)
gtf <- import(GTF_FILE)
gtf_genes <- as.data.frame(gtf) %>%
  filter(type == "gene") %>%
  dplyr::select(gene_id, gene_biotype) %>%
  distinct() %>%
  mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

biotype_recode_map <- c(
  "protein_coding"                     = "Protein Coding",
  "lncRNA"                             = "lncRNA",
  "processed_pseudogene"               = "Proc. Pseudogene",
  "transcribed_processed_pseudogene"   = "Trans. Proc. Pseudogene",
  "transcribed_unprocessed_pseudogene" = "Trans. Unproc. Pseudogene",
  "snoRNA"                             = "snoRNA",
  "miRNA"                              = "miRNA",
  "unannotated"                        = "Unannotated",
  "other"                              = "Other",
  "misc_RNA"                           = "miscRNA",
  "Mt_rRNA"                            = "Mt rRNA",
  "Mt_tRNA"                            = "Mt tRNA"
)

## Ordered color ramp, palest/most-neutral first through boldest/most
## saturated last. 
biotype_color_ramp <- c(
  "#F2E4C9",  # palest cream/tan
  "#C97B63",  # warm terracotta
  "#E0A34E",  # amber/gold
  "#5B7FA6",  # medium blue
  "#A23B32",  # deep red
  "#D9A5A0",  # dusty rose
  "#3D5A80",  # navy
  "#8C6B7A",  # mauve
  "#6B4C5B",  # deep plum
  "#8FAAC7",  # steel blue
  "#B58B5A",  # tan-brown
  "#C9C2B8"   # neutral warm grey
)

## Assign biotype_color_ramp to the levels present in a biotype_summary
## data frame, ranked by total TPM (largest first). Returns a named
## vector suitable for scale_fill_manual(values = ...), plus the
## abundance-ranked level order (for consistent legend/stack ordering).
order_biotype_colors <- function(biotype_summary) {
  ranked_levels <- biotype_summary %>%
    group_by(gene_biotype) %>%
    summarise(total = sum(total_TPM), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(gene_biotype) %>%
    as.character()
  
  list(
    levels = ranked_levels,
    colors = setNames(biotype_color_ramp[seq_along(ranked_levels)], ranked_levels)
  )
}

## ── Whole-transcriptome biotype composition from TPM ──────────
## sample_cols: named vector, e.g. c(A1="NC",A2="NC",A3="NC",B1="CL",B2="CL",B3="CL")
plot_biotype_composition_tpm <- function(tpm_file, sample_cols, comparison_label) {
  
  tpm_df <- read.delim(tpm_file, header = TRUE, sep = "\t", check.names = FALSE)
  tpm_df[ncol(tpm_df)] <- NULL   # drop trailing empty column, matches original file format
  
  tpm_df <- tpm_df %>%
    rename(gene_id = Id) %>%
    mutate(gene_id_base = str_remove(gene_id, "\\.\\d+$"))
  
  tpm_annot <- tpm_df %>%
    left_join(gtf_genes, by = c("gene_id_base" = "gene_id")) %>%
    mutate(gene_biotype = replace_na(gene_biotype, "unannotated"))
  
  tpm_long <- tpm_annot %>%
    pivot_longer(
      cols      = all_of(names(sample_cols)),
      names_to  = "sample",
      values_to = "TPM"
    ) %>%
    mutate(
      condition = sample_cols[sample],
      rep_label = paste0("Rep ", str_extract(sample, "\\d+"))
    )
  
  biotype_summary <- tpm_long %>%
    group_by(rep_label, condition, gene_biotype) %>%
    summarise(total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    group_by(rep_label, condition) %>%
    mutate(proportion = total_TPM / sum(total_TPM)) %>%
    ungroup() %>%
    mutate(gene_biotype = fct_lump_n(gene_biotype, n = 6, w = total_TPM, other_level = "other")) %>%
    mutate(gene_biotype = recode(as.character(gene_biotype), !!!biotype_recode_map)) %>%
    group_by(rep_label, condition, gene_biotype) %>%
    summarise(total_TPM = sum(total_TPM), proportion = sum(proportion), .groups = "drop")
  
  ## Rank biotypes by actual computed abundance and assign the color
  ## ramp accordingly (largest = palest, smallest = boldest) -- see
  ## order_biotype_colors() above.
  biotype_order <- order_biotype_colors(biotype_summary)
  biotype_summary <- biotype_summary %>%
    mutate(gene_biotype = factor(gene_biotype, levels = biotype_order$levels))
  
  p <- ggplot(biotype_summary, aes(x = rep_label, y = proportion, fill = gene_biotype)) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(~condition, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = biotype_order$colors) +
    labs(title = sprintf("RNA biotype composition -- %s", comparison_label),
         x = NULL, y = "Proportion of total TPM", fill = "Biotype") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey92"))
  
  print(p)
  ggsave(file.path(FIG_DIR, sprintf("biotype_composition_%s.pdf", comparison_label)), p, width = 5, height = 4)
  
  ## Exploratory stats (small n, treat as exploratory)
  conditions_present <- unique(biotype_summary$condition)
  if (length(conditions_present) == 2) {
    biotype_stats <- biotype_summary %>%
      group_by(gene_biotype) %>%
      summarise(
        mean_group1 = mean(total_TPM[condition == conditions_present[1]]),
        mean_group2 = mean(total_TPM[condition == conditions_present[2]]),
        fold_change = mean_group2 / mean_group1,
        p_value = tryCatch(
          wilcox.test(total_TPM[condition == conditions_present[1]],
                      total_TPM[condition == conditions_present[2]])$p.value,
          error = function(e) NA_real_
        ),
        .groups = "drop"
      ) %>%
      arrange(p_value)
    write.csv(biotype_stats, file.path(TABLES_DIR, sprintf("biotype_stats_%s.csv", comparison_label)), row.names = FALSE)
  }
  
  invisible(biotype_summary)
}

## CONFIG: sample-to-condition maps for each comparison
sample_cols_NC_CL  <- c(A1 = LEVEL_NC, A2 = LEVEL_NC, A3 = LEVEL_NC,
                        B1 = LEVEL_CL, B2 = LEVEL_CL, B3 = LEVEL_CL)
sample_cols_NC_INT <- c(A1 = LEVEL_NC,  A2 = LEVEL_NC,  A3 = LEVEL_NC,
                        C1 = LEVEL_INT, C2 = LEVEL_INT, C3 = LEVEL_INT)

plot_biotype_composition_tpm(TPM_NC_CL_FILE,  sample_cols_NC_CL,  "NC_vs_CL")
plot_biotype_composition_tpm(TPM_NC_INT_FILE, sample_cols_NC_INT, "NC_vs_INT")