# data/

## What's on GEO

Processed data was deposited to GEO as **nine per-replicate files**, one
per sample, named:

```
A1_fastq.GRCh38.counts_gene.txt   A2_fastq.GRCh38.counts_gene.txt   A3_fastq.GRCh38.counts_gene.txt
B1_fastq.GRCh38.counts_gene.txt   B2_fastq.GRCh38.counts_gene.txt   B3_fastq.GRCh38.counts_gene.txt
C1_fastq.GRCh38.counts_gene.txt   C2_fastq.GRCh38.counts_gene.txt   C3_fastq.GRCh38.counts_gene.txt
```

(A = NC, B = CL, C = INT.) Each file is plain text, tab-delimited, one
row per gene, with columns **`Id`, `Symbol`, `Count`, `FPKM`, `TPM`** —
counts, TPM, and FPKM are all in the same file, despite the
`counts_gene.txt` name in the filename. `Symbol` isn't used by the
compile step or any script below — the pipeline keys off `Id`
(Ensembl gene ID) throughout, with gene symbols added later via
`AnnotationDbi`/`org.Hs.eg.db` in `02_GO_enrichment_MF.R` and
`03_volcano_plots.R`.

## What the scripts expect

The scripts in `scripts/` don't read these nine files directly — they
expect four **wide** tables (one row per gene, one column per sample),
matching the paths in `scripts/00_setup.R`:

| File | Columns |
|---|---|
| `A-B_compiled_fastq.GRCh38.counts_gene.txt` | gene ID, `A1`, `A2`, `A3`, `B1`, `B2`, `B3` (raw counts, NC vs CL) |
| `A-C_compiled_fastq.GRCh38.counts_gene.txt` | gene ID, `A1`, `A2`, `A3`, `C1`, `C2`, `C3` (raw counts, NC vs INT) |
| `A-B_compiled_fastq.GRCh38.TPM_gene.txt` | gene ID, `A1`, `A2`, `A3`, `B1`, `B2`, `B3` (TPM, NC vs CL) |
| `A-C_compiled_fastq.GRCh38.TPM_gene.txt` | gene ID, `A1`, `A2`, `A3`, `C1`, `C2`, `C3` (TPM, NC vs INT) |

The A replicates (NC) are reused across both the B (CL) and C (INT)
comparisons — three A columns, not six.

## Compiling the GEO files into the four wide tables

This snippet extracts the count and TPM columns from each per-replicate
file, renames them to the sample name, and joins them into the four
wide tables above. Run it once after downloading the nine GEO files
into `data/`, before running `01_deseq2.R` onward.

```r
library(dplyr)
library(purrr)

data_dir <- "data"

GENE_ID_COL <- "Id"      # matches the columns confirmed in these GEO files
COUNT_COL   <- "Count"
TPM_COL     <- "TPM"

read_sample_column <- function(sample, value_col) {
  f <- file.path(data_dir, sprintf("%s_fastq.GRCh38.counts_gene.txt", sample))
  df <- read.delim(f, header = TRUE, check.names = FALSE)
  out <- df[, c(GENE_ID_COL, value_col)]
  names(out) <- c("gene_id", sample)
  out
}

compile_wide_table <- function(samples, value_col, out_file) {
  samples %>%
    map(read_sample_column, value_col = value_col) %>%
    reduce(full_join, by = "gene_id") %>%
    write.table(file.path(data_dir, out_file), sep = "\t", row.names = FALSE, quote = FALSE)
}

compile_wide_table(c("A1","A2","A3","B1","B2","B3"), COUNT_COL, "A-B_compiled_fastq.GRCh38.counts_gene.txt")
compile_wide_table(c("A1","A2","A3","C1","C2","C3"), COUNT_COL, "A-C_compiled_fastq.GRCh38.counts_gene.txt")
compile_wide_table(c("A1","A2","A3","B1","B2","B3"), TPM_COL,   "A-B_compiled_fastq.GRCh38.TPM_gene.txt")
compile_wide_table(c("A1","A2","A3","C1","C2","C3"), TPM_COL,   "A-C_compiled_fastq.GRCh38.TPM_gene.txt")
```

After this runs, `data/` should contain the four compiled files listed
above (FPKM isn't used by any script, so it's fine to leave out of the
compile step).

## Other files

- `Homo_sapiens.GRCh38.115.gtf.gz` — downloaded automatically on first run if missing

## GEO accession

Processed per-replicate count/TPM/FPKM files are available from GEO
accession **GSE338149**
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338149)