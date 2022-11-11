##sc.data.file <- "result/step1/merged.mtx.gz"
##annot.file <- "result/step2/celltypes.txt.gz"
##ct.name <- "Mic"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

sc.data.file <- argv[1] # e.g., "result/step1/merged.mtx.gz"
annot.file <- argv[2]   # e.g., "result/step2/celltypes.txt.gz"
ct.name <- argv[3]      # e.g., "Mic"
out.file <- argv[4]     # e.g., ""

library(data.table)
library(tidyverse)
library(mmutilR)

out.hdr <- gsub(".mtx.gz$", "", out.file)
dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
out.data <- fileset.list(out.hdr)

annot.dt <- fread(annot.file, sep = "\t", header = TRUE)

sc.data <- fileset.list(gsub(".mtx.gz$", "", sc.data.file))
cells <- fread(sc.data$col, header=FALSE, col.names = "cell")
cells[, c("barcode","projid","batch") := tstrsplit(`cell`, split="_")]
cells[, projid := as.integer(projid)]

cells.select <- left_join(cells, annot.dt, by = c("barcode", "projid", "batch"))
r_selected <- unlist(cells.select[celltype == ct.name, .(cell)], use.names=FALSE)

out.data <-
    rcpp_mmutil_copy_selected_columns(sc.data$mtx,
                                      sc.data$row,
                                      sc.data$col,
                                      r_selected,
                                      out.hdr)
