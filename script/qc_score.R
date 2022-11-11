argv <- commandArgs(trailingOnly = TRUE)
hdr <- argv[1]
row.file <- argv[2]
col.file <- argv[3]

library(mmutilR)
library(data.table)
.sc <- fileset.list(hdr)
.scores <- rcpp_mmutil_compute_scores(.sc$mtx, .sc$row, .sc$col)

fwrite(setDT(.scores$row), file=row.file)
fwrite(setDT(.scores$col), file=col.file)
