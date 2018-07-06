---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# gscars

The goal of gscars is to ...

## Installation

You can install gscars from github with:


```r
# install.packages("devtools")
devtools::install_github("mcieslik-mctp/gscars")
```

## Example

This is a basic example which shows you how to solve a common problem:


```r
library(gscars)
## initial barcode counting based on the EMA algorithm:
countBarcodes("inst/extdata/test_1.fq", "barcode_counts.bin", "inst/extdata/4M-with-alts-february-2016.txt")
#> [1] "barcode_counts.bin"
## preprocess fastq files for alignment using BWA
preprocessFastq("barcode_counts.bin", "inst/extdata/test_1.fq", "inst/extdata/test_2.fq", "1.fq", "2.fq")
#> [1] "1.fq" "2.fq"
```
