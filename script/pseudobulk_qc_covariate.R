argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## expr.file <- "result/step3/pb/Mic.rds"
## feature.file <- "result/step1/features_annotated_GRCh37.txt.gz"
## pheno.file <- "data/metadata_PFC_all_individuals_092520.tsv.gz"
## NPC <- 75
## out.dir <- "temp/Mic"

if(length(argv) < 4) q()

expr.file <- argv[1]    # e.g., "result/step3/pb/Mic.rds"
feature.file <- argv[2] # e.g., "result/step1/features_annotated_GRCh37.txt.gz"
NPC <- as.integer(argv[3])
out.dir <- argv[4]

library(data.table)
library(dplyr)
library(rsvd)
library(preprocessCore) # quantile norm

.mkdir <- function(...) {
    dir.create(..., recursive = TRUE, showWarnings = FALSE)
}

`%&%` <- function(a,b) paste0(a,b)

.safe.lm <- function(Y, C){
    Y.resid <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Y.fitted <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for(j in 1:ncol(Y)){
        y.j <- Y[, j]
        if(sum(is.finite(y.j)) < 2) next
        .lm <- lm(y.j ~ C, na.action = "na.exclude")
        Y.resid[,j] <- residuals(.lm)
        Y.fitted[,j] <- fitted(.lm)
    }
    list(fitted = Y.fitted, residuals = Y.resid)
}

.quantile.norm.std <- function(.mat) {
    stopifnot(is.matrix(.mat))
    ret <- .mat
    for(k in 1:ncol(ret)){
        x.k <- .mat[, k]
        .pos.k <- which(is.finite(x.k))
        x.k.valid <- x.k[.pos.k]
        ngenes <- length(x.k.valid)
        if(ngenes < 1) next
        qq <- qnorm((1:ngenes)/(ngenes + 1))
        x.k.valid[order(x.k.valid)] <- qq
        ret[.pos.k, k] <- x.k.valid
    }
    return(ret)
}

take.marginal.stat <- function(xx, yy, se.min=1e-8) {
    .xx <- apply(xx, 2, scale)
    .yy <- apply(yy, 2, scale)
    rm.na.zero <- function(xx) {
        return(replace(xx, is.na(xx), 0))
    }
    n.obs <- crossprod(!is.na(.xx), !is.na(.yy))
    beta.mat <- crossprod(rm.na.zero(.xx), rm.na.zero(.yy))/n.obs
    ## calibrate residuals
    resid.se.mat <- matrix(NA, ncol(.xx), ncol(.yy))
    for (k in 1:ncol(.yy)) {
        beta.k <- beta.mat[, k]
        yy.k <- .yy[, k]
        err.k <- sweep(sweep(.xx, 2, beta.k, `*`), 1, yy.k, `-`)
        se.k <- apply(err.k, 2, sd, na.rm = TRUE)
        resid.se.mat[, k] <- se.k + se.min
    }
    y.cols <- 1:ncol(yy)
    colnames(beta.mat) <- y.cols
    colnames(n.obs) <- y.cols
    colnames(resid.se.mat) <- y.cols
    rownames(resid.se.mat) <- rownames(beta.mat)
    ## combine the results
    .melt.mat <- function(.mat, ...) {
        rownames(.mat) <- 1:nrow(.mat)
        colnames(.mat) <- 1:ncol(.mat)
        reshape2::melt(.mat, ...) %>%
            dplyr::rename(x.col = Var1, y.col = Var2) %>%
            as.data.table()
    }
    zscore.pvalue <- function(z) {
        2*pnorm(abs(z), lower.tail = FALSE)
    }
    beta.dt <- .melt.mat(beta.mat, value.name = "beta")
    resid.se.dt <- .melt.mat(resid.se.mat, value.name = "resid.se")
    nobs.dt <- .melt.mat(n.obs, value.name = "n")
    out <- beta.dt %>%
        left_join(nobs.dt, by = c("x.col", "y.col")) %>%
        left_join(resid.se.dt, by = c("x.col", "y.col")) %>%
        dplyr::mutate(se = resid.se/sqrt(n)) %>%
        dplyr::mutate(p.val = zscore.pvalue(beta/se)) %>%
        dplyr::mutate(beta = round(beta, 4)) %>%
        dplyr::select(-resid.se) %>%
        dplyr::mutate(se = round(se, 4)) %>%
        dplyr::mutate(x.col = as.integer(x.col)) %>%
        dplyr::mutate(y.col = as.integer(y.col))

    return(out)
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

## Remove genes with too few expression values
## - we would consider zero expression is missing
filter.mat <- function(.mat, .sum, missing.cutoff = 1, missing.rate = 0.8) {
    .row.missing.rate <- apply(t(.sum) < missing.cutoff, 2, mean)
    .qc <- .mat
    .qc[.sum < missing.cutoff] <- NA
    .qc[.row.missing.rate < missing.rate, , drop = FALSE]
}

################################################################
## Mon, 20 Feb 2023 22:54:01
##
## Would it be possible for you to still run the method with
## 1) no PCs,
## 2) with all PCs
## 3) with selected PCs?
## --> Remove PCs if correlated with AD-related vars.
##
## Potential causal DAG:
##
##       Environmental
##           |
##           v
## APOE --> gene --> brain disorders
##            ^
##            |
##       other variants

.mkdir(out.dir)

ct <- gsub(".rds$","",basename(expr.file))
out.hdr <- out.dir %&% "/" %&% ct

out.files <- list(
    vanilla = out.hdr %&% "_allIndv_noPC.bed.gz",
    all.pc = out.hdr %&% "_allIndv_allPC.bed.gz",
    selected.pc = out.hdr %&% "_allIndv_selectedPC.bed.gz",
    cond.noAD = out.hdr %&% "_noAD_selectedPC.bed.gz",
    cond.AD = out.hdr %&% "_AD_selectedPC.bed.gz",
    cond.male = out.hdr %&% "_male_selectedPC.bed.gz",
    cond.female = out.hdr %&% "_female_selectedPC.bed.gz",
    pcs = out.hdr %&% "_PC.csv.gz",
    associations = out.hdr %&% "_associations.csv.gz"
)

feature.info <- fread(feature.file)
expr <- readRDS(expr.file)

message("Read expression data")

mu.qc <- filter.mat(expr$PB$ln.mu, expr$PB$sum, missing.cutoff = .5)

## Remove column-wise bias
.bias <- apply(mu.qc, 2, mean, na.rm=TRUE)
mu.qc.adj <- sweep(mu.qc, 2, .bias, `-`) %>%
    .sort.cols(pheno = expr$pheno)

mu.qc.std <- t(apply(t(mu.qc.adj), 2, scale))
rownames(mu.qc.std) <- rownames(mu.qc.adj)
colnames(mu.qc.std) <- colnames(mu.qc.adj)

message("Basic Q/C to remove column-wise bias")

## 3. compute principal components filling in zero
mu.qc.nz <- mu.qc.std
mu.qc.nz[is.na(mu.qc.nz)] <- min(mu.qc.std, na.rm=TRUE)
X <- apply(t(mu.qc.nz), 2, scale) ## sample x gene
.svd <- rsvd(X, k = NPC)

## read meta data
pheno <- expr$pheno[, .(study,
                        `AD`,
                        `APOE`,
                        `cognep_random_slope`,
                        `age_death`,
                        `educ`,
                        `msex`,
                        `braaksc`,
                        `amyloid`,
                        `nft`,
                        `gpath`,
                        `pmi`)] %>%
    mutate(study = if_else(`study` == "ROS", 0, 1))

phi <- as.matrix(pheno)

pheno.assoc <- take.marginal.stat(phi, .svd$u) %>%
    left_join(data.table(y.col = 1:NPC, pc = 1:NPC)) %>%
    left_join(data.table(x.col = 1:ncol(phi), pheno = colnames(phi))) %>%
    select(-x.col, -y.col)

pheno.assoc[order(p.val), head(.SD, 1), by = .(pheno)]

disease.vars <- c("AD",
                  "APOE",
                  "cognep_random_slope",
                  "braaksc",
                  "amyloid",
                  "nft",
                  "gpath")

disease.pcs <-
    unlist(pheno.assoc[pheno %in% disease.vars & p.val < 0.05, .(pc)]) %>%
    unique()

################################################################
## 1) no PCs,
Y <- mu.qc.std

## 2) with all PCs
C <- .svd$u
Y.allPC <- t(.safe.lm(t(Y), C)$residuals)
colnames(Y.allPC) <- colnames(Y)
rownames(Y.allPC) <- rownames(Y)

message("Adjusted all the PCs")

## 3) with selected PCs?
## --> Remove PCs if correlated with AD-related vars.
C.safe <- C
if(length(disease.pcs) > 0){
    C.safe <- C[, -unique(disease.pcs), drop = FALSE]
}
Y.selectedPC <- t(.safe.lm(t(Y), C.safe)$residuals)
colnames(Y.selectedPC) <- colnames(Y)
rownames(Y.selectedPC) <- rownames(Y)

message("Adjusted some of the selected PCs")

.sort.rows(Y, feature.info) %>%
    write.bed.gz(out.files$vanilla)

.sort.rows(Y.allPC, feature.info) %>%
    write.bed.gz(out.files$all.pc)

.sort.rows(Y.selectedPC, feature.info) %>%
    write.bed.gz(out.files$selected.pc)

.sort.rows(Y.selectedPC[, expr$pheno$AD == 1, drop = FALSE], feature.info) %>%
    write.bed.gz(out.files$cond.AD)

.sort.rows(Y.selectedPC[, expr$pheno$AD == 0, drop = FALSE], feature.info) %>%
    write.bed.gz(out.files$cond.noAD)

.sort.rows(Y.selectedPC[, expr$pheno$msex == 1, drop = FALSE], feature.info) %>%
    write.bed.gz(out.files$cond.male)

.sort.rows(Y.selectedPC[, expr$pheno$msex == 0, drop = FALSE], feature.info) %>%
    write.bed.gz(out.files$cond.female)

fwrite(cbind(colnames(Y), as.data.table(.svd$u)), out.files$pcs)
fwrite(pheno.assoc, out.files$associations)

message("Done")
