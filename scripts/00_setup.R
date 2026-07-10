## ============================================================
## 00_setup.R
## Shared configuration: paths, libraries, condition labels.
## Every other script in scripts/ sources this file first.
## ============================================================

## ── Libraries ────────────────────────────────────────────────
## (loaded here so every script has a consistent, complete set;
##  individual scripts do not re-library() these)
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(rtracklayer)
  library(biomaRt)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(cowplot)
  library(ggVennDiagram)
  library(EnhancedVolcano)
  library(tibble)
})

source(here::here("scripts", "utils.R"))

## ── Project paths ────────────────────────────────────────────
## CONFIG: update these if your folder layout differs.
## Everything is relative to the project root, so the project
## should be opened as an RStudio Project (.Rproj) or run with
## working directory set to the repo root.
DATA_DIR    <- here::here("data")
TABLES_DIR  <- here::here("results", "tables")
FIG_DIR     <- here::here("results", "figures")

dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)

## ── Raw count / TPM file paths ───────────────────────────────
## CONFIG: point these at your actual files inside data/.
COUNTS_NC_CL_FILE  <- file.path(DATA_DIR, "A-B_compiled_fastq.GRCh38.counts_gene.txt")
COUNTS_NC_INT_FILE <- file.path(DATA_DIR, "A-C_compiled_fastq.GRCh38.counts_gene.txt")

TPM_NC_CL_FILE  <- file.path(DATA_DIR, "A-B_compiled_fastq.GRCh38.TPM_gene.txt")
TPM_NC_INT_FILE <- file.path(DATA_DIR, "A-C_compiled_fastq.GRCh38.TPM_gene.txt")

GTF_FILE <- file.path(DATA_DIR, "Homo_sapiens.GRCh38.115.gtf.gz")
GTF_URL  <- "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"

## ── RDS locations for consolidated DESeq2 output ─────────────
## 01_deseq2.R writes these; every downstream script reads them
## instead of re-running DESeq2.
RDS_NC_CL  <- file.path(TABLES_DIR, "res_NC_vs_CL.rds")
RDS_NC_INT <- file.path(TABLES_DIR, "res_NC_vs_INT.rds")

## ── Condition labels (standardized across all scripts) ───────
LEVEL_NC  <- "NC"   # noncrosslinked (control), samples A1-A3
LEVEL_CL  <- "CL"   # crosslinked, samples B1-B3
LEVEL_INT <- "INT"  # APGC interphase, samples C1-C3

## ── Shared thresholds ─────────────────────────────────────────
PADJ_THRESH <- 0.05
LFC_THRESH  <- 1
