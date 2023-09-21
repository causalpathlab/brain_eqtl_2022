## data.file <- "result/step2/sorted/Mic.mtx.gz"
## out.file <- "result/step3/pb/Mic.rds"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

data.file <- argv[1] # e.g., data.file <- "result/step2/sorted/Mic.mtx.gz"
out.file <- argv[2]  # e.g., out.file <- "result/step3/Mic.RDS"

if(file.exists(out.file)) q()

dir.create(dirname(out.file),
           recursive = TRUE,
           showWarnings = FALSE)

ct <- gsub(".mtx.gz$","",basename(data.file))

library(mmutilR)
library(data.table)
library(dplyr)
setDTthreads(1)

get.cols <- function(.data){
    .cols <- fread(.data$col, header=FALSE, col.names="V1")
    .cols[, c("barcode", "projid", "batch") := tstrsplit(`V1`, split="_")]
    .cols[, projid := as.integer(`projid`)]
    return(.cols)
}

.data <- fileset.list(gsub(".mtx.gz$", "", data.file))
.cols <- get.cols(.data)
c2indv <- as.data.frame(.cols[, c(1, 3)])
pb.data <- make.cocoa(.data, ct, cell2indv = c2indv)

saveRDS(pb.data, out.file)
