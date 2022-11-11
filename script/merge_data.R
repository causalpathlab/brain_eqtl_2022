argv <- commandArgs(trailingOnly = TRUE)

data.dir <- argv[1] # "result/step1/qc/"
out.hdr <- argv[2]  # "result/step1/merged"

library(dplyr)
library(stringr)

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

.batches <- list.files(data.dir, pattern = ".mtx.gz$") %>%
    str_remove(".mtx.gz$") %>%
    str_remove("^SM_")

.headers <-
    list.files(data.dir, pattern = ".mtx.gz$", full.names = TRUE) %>%
    str_remove(".mtx.gz$")

.data <-
    mmutilR::rcpp_mmutil_merge_file_sets(
                 r_headers=.headers,
                 r_batches=.batches,
                 r_mtx=str_c(.headers, ".mtx.gz"),
                 r_row=str_c(.headers, ".rows.gz"),
                 r_col=str_c(.headers, ".cols.gz"),
                 output=out.hdr)
