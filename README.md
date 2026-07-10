# polyA-UV-RNAseq

Citation information will be added upon publication. This manuscript is currently under review.

In this study, we analyzed three main datasets. (1) Poly(A) selected RNA from non-crosslinked HeLa cells, (2) Poly(A) selected RNA from UV-crosslinked HeLa cells, and (3) Poly(A) selected RNA from the AGPC interphase of UV-crosslinked HeLa cells. 

This repository contains the analysis scripts used to generate the
RNA-seq figures and tables in [citation/DOI here]. It covers
differential expression, RNA biotype composition, transcript length
analysis, GO Molecular Function enrichment, and replicate correlation
for two comparisons:

- **NC vs CL**: noncrosslinked control vs crosslinked
- **NC vs INT**: noncrosslinked control vs APGC interphase

NC is the reference/control condition in both comparisons. Each
condition has 3 biological replicates.

## Citation

<!-- Add: full citation once published, e.g. -->
<!-- Author et al. (Year). Title. Journal. DOI: xx.xxxx/xxxxx -->


## Data availability

Processed per-replicate count/TPM/FPKM files are available from GEO
accession **GSE338149**
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338149).

See `data/README.md` for the exact file format and a script to compile
the nine GEO files into the tables this pipeline expects.

## Repository structure

```
.
├── README.md
├── data/                    # input files (not tracked in git — see data/README.md)
├── scripts/
│   ├── 00_setup.R           # paths, libraries, condition labels — sourced by every other script
│   ├── utils.R              # shared helper functions
│   ├── 01_deseq2.R          # runs DESeq2 once per comparison, saves .rds
│   ├── 02_GO_enrichment_MF.R
│   ├── 03_volcano_plots.R     # labels genes from significant GO MF terms
│   ├── 04_biotype_composition.R
│   ├── 05_length_analysis.R
│   ├── 06_quadrant_plot.R
│   └── 07_correlation_plots.R
```

## Requirements

- R version [4.4.0]
- Packages: DESeq2, tidyverse, dplyr, ggplot2, ggrepel, scales,
  RColorBrewer, rtracklayer, biomaRt, org.Hs.eg.db, AnnotationDbi,
  clusterProfiler, cowplot, ggVennDiagram, EnhancedVolcano, tibble, here

## Setup

1. Clone this repository.
2. Download the processed data files from GEO (accession above) and
   compile them into the four tables the scripts expect — see
   `data/README.md` for the exact steps and a ready-to-run compile
   script.
3. Open this folder as an RStudio Project (or set your working
   directory to the repo root) — all paths are relative to the project
   root via the `here` package.

## Running the pipeline

Scripts are numbered in the order they should be run — each one after
`01_deseq2.R` reads its `.rds` output rather than re-running DESeq2:

```r
source("scripts/01_deseq2.R")
source("scripts/02_GO_enrichment_MF.R")   # must run before 03 and 06
source("scripts/03_volcano_plots.R")      # labels GO-term genes from 02
source("scripts/04_biotype_composition.R")
source("scripts/05_length_analysis.R")
source("scripts/06_quadrant_plot.R")
source("scripts/07_correlation_plots.R")  # independent of the others
```

Figures are written to `results/figures/`, tables to `results/tables/`.

## License and AI Disclosure

### Code Origin
The R scripts in this repository were generated using Claude Sonnet 5
(Anthropic). The code has been verified by the authors to accurately
reproduce the figures in the main text.

### License
To maximize reproducibility and remove all legal barriers, this
repository is made available under the [MIT License](LICENSE).
