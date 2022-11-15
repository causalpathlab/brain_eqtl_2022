argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## expr.file <- "result/step3/pb/Exc-L3-5-RORB-PLCH1.rds"
## feature.file <- "result/step1/features_annotated_GRCh37.txt.gz"
## NPC <- 50
## out.dir <- "temp/Mic"

if(length(argv) < 4) q()

expr.file <- argv[1]    # e.g., "result/step3/pb/Mic.rds"
feature.file <- argv[2] # e.g., "result/step1/features_annotated_GRCh37.txt.gz"
NPC <- as.integer(argv[3])
out.dir <- argv[4]

library(data.table)
library(dplyr)
library(rsvd)

.mkdir <- function(...) {
    dir.create(..., recursive = TRUE, showWarnings = FALSE)
}

`%&%` <- function(a,b) paste0(a,b)

.safe.lm <- function(Y, C){
    Y.resid <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Y.fitted <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for(j in 1:ncol(Y)){
        .lm <- lm(Y[, j] ~ C, na.action = "na.exclude")
        Y.resid[,j] <- residuals(.lm)
        Y.fitted[,j] <- fitted(.lm)
    }
    list(fitted = Y.fitted, residuals = Y.resid)
}

.quantile.norm <- function(.mat) {
    stopifnot(is.matrix(.mat))
    ngenes <- nrow(.mat)
    qq <- qnorm((1:ngenes)/(ngenes + 1))
    .rows <- rownames(.mat)
    ret <- .mat
    for(k in 1:ncol(ret)){
        ret[order(ret[,k]), k] <- qq
    }
    return(ret)
}


run.qc <- function(expr.mat,
                   feature.info,
                   num.pc = NPC,
                   do.quantile.norm = TRUE){

    if(do.quantile.norm){
        expr.mat <- .quantile.norm(expr.mat)
    }

    features <- data.table(row=rownames(expr.mat))
    features[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`row`, split="_")]
    features <- left_join(features, feature.info)

    if(num.pc > 0){
        expr.qc <- matrix(NA, nrow(expr.mat), ncol(expr.mat))
        rownames(expr.qc) <- rownames(expr.mat)
        colnames(expr.qc) <- colnames(expr.mat)

        ## adjust confounding effects
        for(chr in 1:22) {
            chr.loc <- which(features$chromosome_name == chr)
            x1 <- expr.mat[chr.loc, , drop = FALSE]
            x0 <- expr.mat[-chr.loc, , drop = FALSE]
            x0[is.na(x0)] <- 0
            chr.svd <- rsvd::rsvd(x0, k = pmin(num.pc, ncol(x0) - 1))
            x1.resid <- .safe.lm(t(x1), chr.svd$v)$residuals
            expr.qc[chr.loc, ] <- t(x1.resid)
        }
        return(expr.qc)
    } else {
        return(expr.mat)
    }
}

.sort.cols <- function(.mat, pheno) {
    .dt <- data.table(col = colnames(.mat))
    .dt[, c("projid", "ct") := tstrsplit(`col`, split="_")]
    .loc <- match(pheno$projid, as.integer(.dt$projid))
    .mat[, .loc, drop = FALSE]
}

.sort.col.pairs <- function(.mat, .knn){
    .pairs <- data.table(pp = colnames(.mat))
    .pairs[, c("obs.name", "matched.name", "ct") := tstrsplit(`pp`,split="_")]
    .pairs[, .loc := 1:.N]
    .sorted <- left_join(.knn, .pairs)
    .mat[, .sorted$.loc, drop = FALSE]
}

.sort.rows <- function(.mat, feature.info){

    rows <- data.table(row = rownames(.mat))
    rows[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`row`, split="_")]
    rows <- rows %>%
        mutate(r = 1:n()) %>%
        left_join(feature.info) %>%
        filter(`chromosome_name` %in% 1:22) %>%
        as.data.table()

    ret <- rows %>%
        select(`chromosome_name`, `tss`, `tes`,
               `ensembl_gene_id`, `hgnc_symbol`) %>%
        cbind(.mat[rows$r, , drop = FALSE]) %>%
        arrange(`chromosome_name`, `tss`) %>%
        rename(`#chromosome_name` = `chromosome_name`) %>% 
        as.data.table()

    return(ret)
}

write.bed.gz <- function(.dt, out.file){
    if(!file.exists(out.file)){
        out.raw.file <- gsub(".gz$", "", out.file)
        fwrite(.dt, out.raw.file, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
        ## Rsamtools::bgzip(out.raw.file, dest=out.file)
        ## Rsamtools::indexTabix(out.file, format="vcf")
        system(paste0("bgzip ", out.raw.file))
        system(paste0("tabix -p bed ", out.file))
        unlink(out.raw.file)
        message("Wrote: ", out.file)
    }
}

################################################################

.mkdir(out.dir)

ct <- gsub(".rds$","",basename(expr.file))
out.hdr <- out.dir %&% "/" %&% ct %&% "_" %&% "PC" %&% NPC
ad.hdr <- out.dir %&% "/" %&% ct %&% "_" %&% "AD_all"
pine.hdr <- out.dir %&% "/" %&% ct %&% "_" %&% "PINE_all"

out.files <- list(mu = out.hdr %&% "_all.bed.gz",
                  cond.AD = out.hdr %&% "_AD.bed.gz",
                  cond.noAD = out.hdr %&% "_noAD.bed.gz",
                  cond.APOE = out.hdr %&% "_APOE.bed.gz",
                  cond.noAPOE = out.hdr %&% "_noAPOE.bed.gz",
                  cond.male = out.hdr %&% "_male.bed.gz",
                  cond.female = out.hdr %&% "_female.bed.gz",
                  PINE = pine.hdr %&% ".bed.gz",
                  AD = ad.hdr %&% ".bed.gz")

if(all(file.exists(unlist(out.files)))) {
    message("all files exist!")
    q()
}

feature.info <- fread(feature.file)

expr <- readRDS(expr.file)

message("Read expression data")

mu.mat <-
    run.qc(expr$PB$mu, feature.info, NPC, TRUE) %>% 
    .sort.cols(pheno = expr$pheno)

mu.dt <- mu.mat %>% 
    .sort.rows(feature.info)

write.bed.gz(mu.dt, out.files$mu)

message("Average expression matrix with PC correction, ", NPC)

ad.dt <-
    run.qc(expr$AD$resid.mu, feature.info, 0, TRUE) %>%
    .sort.cols(pheno = expr$pheno) %>%
    .sort.rows(feature.info)

write.bed.gz(ad.dt, out.files$AD)

message("AD-specific expression matrix with no PC correction")

################################################################

mu.mat[, expr$pheno$AD == 1, drop = FALSE] %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$cond.AD)

mu.mat[, expr$pheno$AD == 0, drop = FALSE] %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$cond.noAD)

mu.mat[, expr$pheno$APOE == 1, drop = FALSE] %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$cond.APOE)

mu.mat[, expr$pheno$APOE == 0, drop = FALSE] %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$cond.noAPOE)

mu.mat[, expr$pheno$msex == 1, drop = FALSE] %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$cond.male)

mu.mat[, expr$pheno$msex == 0, drop = FALSE] %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$cond.female)

message("Built condition-specific expression matrices after PC correction")

################################################################

knn.info <- setDT(expr$PINE$knn)

run.qc(expr$PINE$delta, feature.info, 0, TRUE) %>% 
    .sort.col.pairs(knn.info) %>%
    .sort.rows(feature.info) %>%
    write.bed.gz(out.files$PINE)

message("Done")
