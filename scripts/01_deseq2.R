## ============================================================
## 01_deseq2.R
## Runs DESeq2 ONCE per comparison and saves the results object
## to results/tables/ as an .rds file. Every other script reads
## these .rds files instead of re-running DESeq2.
## Comparisons:
##   NC vs CL  (A vs B): noncrosslinked control vs crosslinked
##   NC vs INT (A vs C): noncrosslinked control vs APGC interphase
## In both cases NC (A) is the reference level.
## ============================================================

source(here::here("scripts", "00_setup.R"))

## ── NC vs CL (A-B) ────────────────────────────────────────────
res_NC_CL <- run_deseq2(
  COUNTS_NC_CL_FILE,
  conditions = c(LEVEL_NC, LEVEL_NC, LEVEL_NC, LEVEL_CL, LEVEL_CL, LEVEL_CL)
)
cat("=== NC vs CL ===\n")
summary(res_NC_CL)
saveRDS(res_NC_CL, RDS_NC_CL)

## ── NC vs INT (A-C) ───────────────────────────────────────────
res_NC_INT <- run_deseq2(
  COUNTS_NC_INT_FILE,
  conditions = c(LEVEL_NC, LEVEL_NC, LEVEL_NC, LEVEL_INT, LEVEL_INT, LEVEL_INT)
)
cat("=== NC vs INT ===\n")
summary(res_NC_INT)
saveRDS(res_NC_INT, RDS_NC_INT)

cat(sprintf("\nSaved:\n  %s\n  %s\n", RDS_NC_CL, RDS_NC_INT))
