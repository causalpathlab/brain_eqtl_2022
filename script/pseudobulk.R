## data.file <- "result/step2/sorted/Mic.mtx.gz"
## pheno.file <- "data/metadata_PFC_all_individuals_092520.tsv.gz"
## out.file <- "result/step3/Mic.RDS"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

data.file <- argv[1] # e.g., data.file <- "result/step2/sorted/Mic.mtx.gz"
pheno.file <- argv[2] # e.g., pheno.file <- "data/metadata_PFC_all_individuals_092520.tsv.gz"
out.file <- argv[3]  # e.g., out.file <- "result/step3/Mic.RDS"

if(file.exists(out.file)) q()
dir.create(dirname(out.file),
           recursive = TRUE,
           showWarnings = FALSE)

ct <- gsub(".mtx.gz$","",basename(data.file))

library(mmutilR)
library(data.table)
library(dplyr)
setDTthreads(1)

pheno.dt <- fread(pheno.file, header=TRUE)

get.cols <- function(.data){
    .cols <- fread(.data$col, header=FALSE, col.names="V1")
    .cols[, c("barcode", "projid", "batch") := tstrsplit(`V1`, split="_")]
    .cols[, projid := as.integer(`projid`)]
    return(.cols)
}

.data <- fileset.list(gsub(".mtx.gz$", "", data.file))
.cols <- get.cols(.data)

scores <- rcpp_mmutil_compute_scores(.data$mtx, .data$row, .data$col)
cell.scores <- setDT(scores$col)
gene.scores <- setDT(scores$row)

######################################
## match with phenotype information ##
######################################

.pheno.cols <- get.cols(.data) %>% left_join(pheno.dt)

indv.pheno <- .pheno.cols[, head(.SD, 1), by = .(projid)]
indv.pheno[, AD := 0]
indv.pheno[braaksc >= 4, AD := 1] # BRAAK Stage 4,5,6 -> AD in PFC
indv.pheno[, APOE := 0]
indv.pheno[apoe_genotype %in% c("24","34","44"), APOE := 1] # apoe4 carrier

#################################
## 1. vanilla pseudo-bulk data ##
#################################

c2indv <- as.data.frame(.cols[, c(1, 3)])
pb.data <- make.cocoa(.data, ct, cell2indv = c2indv)

out.indv <- data.table(sample = colnames(pb.data$sum))
out.indv[, c("projid","celltype") := tstrsplit(`sample`,split="_")]
out.indv[, projid := as.integer(projid)]
out.indv <- left_join(out.indv, indv.pheno) %>% as.data.table()

#############################
## 2. AD-aware pseudo-bulk ##
#############################

ad.data <- make.cocoa(.data, ct, cell2indv = c2indv,
                      indv2exp = indv.pheno[, .(projid, AD)],
                      knn=50, .rank=50, impute.by.knn = TRUE,
                      .em.iter=10, num.threads = 1)

##########################################
## 3. pairwise individual matching data ##
##########################################

pine.data <- make.pine(.data, ct,
                       cell2indv = c2indv,
                       knn.cell = 50, knn.indv = 1,
                       .rank = 50, .em.iter=10, num.threads = 1,
                       impute.by.knn = TRUE)

out.list <- list(PB = pb.data,
                 AD = ad.data,
                 PINE = pine.data,
                 pheno = out.indv)

saveRDS(out.list, out.file)
