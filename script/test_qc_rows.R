argv <- commandArgs(trailingOnly = TRUE)

row.score.dir <- argv[1]
out.file <- argv[2]

library(data.table)
library(dplyr)
library(stringr)

.files <- list.files(row.score.dir, full.names = TRUE, pattern = ".rows.gz")

## Find features to remove if
## (1) SD < .01 or
## (2) mean < 1e-4
## (3) nnz < 10
.remove <-
    lapply(.files, fread) %>%
    lapply(function(x) x[`sd` < 1e-2 | `mean` < 1e-4 | `nnz` < 10, .(`name`)]) %>%
    do.call(what = rbind) %>%
    unique() %>%
    unlist()

.genes <- lapply(.files, fread) %>%
    do.call(what = rbind) %>%
    filter(!(`name` %in% .remove)) %>%
    select(`name`) %>%
    unique

fwrite(.genes, out.file, col.names = FALSE)
