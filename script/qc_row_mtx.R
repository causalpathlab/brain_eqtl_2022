argv <- commandArgs(trailingOnly = TRUE)

valid.row.file <- argv[1] # "result/step1/valid_features.txt.gz"
data.hdr <- argv[2]       # "result/step1/temp/SM_IVA6E"
qc.hdr <- argv[3]         # output

dir.create(dirname(qc.hdr), recursive = TRUE, showWarnings = FALSE)

.rows <- unlist(data.table::fread(valid.row.file, header = FALSE))

library(mmutilR)

.raw <- fileset.list(data.hdr)

.qc.data <- rcpp_mmutil_copy_selected_rows(
    mtx_file = .raw$mtx,
    row_file = .raw$row,
    col_file = .raw$col,
    r_selected = .rows,
    output = qc.hdr)

