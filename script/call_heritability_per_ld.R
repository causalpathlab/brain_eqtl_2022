argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## nPC <- 10
## geno.hdr <- "result/step4/rosmap"
## expr.file <- "result/step3/log_mean.bed.gz"
## svd.file <- "result/step3/svd.rds"
## out.file <- "output.txt.gz"

if(length(argv) < 7) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
expr.file <- argv[4]
svd.file <- argv[5]
nPC <- as.integer(argv[6])
out.file <- argv[7]

temp.dir <- paste0(out.file, "_temp")

################################################################

library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)

`%&%` <- function(a,b) paste0(a,b)

subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

    .error <- function(e) {
        print(e)
        cat("Failed to read plink!\n", file=stderr())
        return(NULL)
    }

    dir.create(temp.dir, recursive=TRUE, showWarnings=FALSE)

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num <- as.integer(gsub(pattern="chr", replacement="", chr))

        plink.cmd <- sprintf("./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --hwe 1e-4 --chr %d --from-bp %d --to-bp %d --out %s",
                             plink.hdr, chr.num, plink.lb, plink.ub, paste0(temp.dir, "/plink"))
        system(plink.cmd)

        unlink(paste0(temp.dir, "/plink.bk"))
        .bed <- bigsnpr::snp_readBed(paste0(temp.dir, "/plink.bed"))
        plink <- bigsnpr::snp_attach(.bed)
        ret <- list(bed=plink$genotypes[,],
                    fam=plink$fam,
                    map=plink$map)
        unlink(paste0(temp.dir, "/plink.bk"))
        return(ret)
    }

    plink <- tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error=.error)
    return(plink)
}

.safe.lm <- function(Y, C){
    Y.resid <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Y.fitted <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for(j in 1:ncol(Y)){
        .lm <- lm(Y[, j] ~ C, na.action = "na.exclude")
        Y.resid[,j] <- residuals(.lm)
        Y.fitted[,j] <- fitted(.lm)
    }
    rownames(Y.resid) <- rownames(Y)
    rownames(Y.fitted) <- rownames(Y)
    colnames(Y.resid) <- colnames(Y)
    colnames(Y.fitted) <- colnames(Y)
    list(fitted = Y.fitted, residuals = Y.resid)
}

match.with.plink <- function(Y, plink){

    .y.info <- data.table(projid = as.integer(rownames(Y)))
    .y.info[, y.row := 1:.N]

    .match <-
        plink$fam %>%
        mutate(projid = as.integer(sample.ID)) %>%
        mutate(x.row = 1:n()) %>%
        left_join(.y.info, by = "projid") %>%
        na.omit()

    X <- plink$bed

    ret <- list(x = X[.match$x.row, , drop = FALSE],
                y = Y[.match$y.row, , drop = FALSE])

    return(ret)
}

.quantile.norm <- function(.mat) {
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
    rownames(ret) <- rownames(.mat)
    colnames(ret) <- colnames(.mat)
    return(ret)
}

################################################################

svd.theta <- function(.svd, k, zz){
    lambda <- 1/nrow(.svd$u)
    vd2inv <- sweep(.svd$v[, 1:k, drop = F], 2, 1/(.svd$d[1:k]^2 + lambda), `*`)
    vz.k <- t(.svd$v[, 1:k, drop = F]) %*% zz
    .theta.k <- vd2inv %*% vz.k
}

fast.z <- function (x, y)
{
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0), replace(y, is.na(y), 0))/sqrt(n.obs)
    ret[is.na(ret)] <- 0
    return(ret)
}

check.sd <- function(x) {
    .sd <- sd(x, na.rm = T)
    if(is.na(.sd) || .sd <= 0) return(FALSE)
    return(TRUE)
}

safe.corr <- function(y1, y2, ...){
    if(check.sd(y1) && check.sd(y2)){
        rr <- cor(y1, y2, ...)
        return(rr)
    } else {
        return(NA)
    }
}

safe.pve <- function(y.target, y.pred){

    if(check.sd(y.target) && check.sd(y.pred)){
        .lm <- .safe.lm(matrix(y.target), matrix(y.pred))
        .pve <- var(.lm$fitted, na.rm=T)/var(y.target, na.rm=T)
        return(.pve)
    } else {
        return(NA)
    }
}

cv.svd.comp <- function(xx, yy,
                        nfold = 5,
                        ncomp = c(5, 10, 20, 30, 40, 50),
                        ...){

    nn <- nrow(xx)
    pp <- ncol(xx)

    cv.idx <- sample(1:nfold, nn, replace=TRUE)
    cv.result <- data.table()

    for(r in 1:nfold){
        .test <- which(cv.idx == r)
        .train <- which(cv.idx != r)
        xx.test <- xx[.test, , drop = F]
        yy.test <- yy[.test, , drop = F]
        xx.train <- xx[.train, , drop = F]
        yy.train <- yy[.train, , drop = F]

        .svd <- rsvd::rsvd(xx.train, k = max(ncomp))

        for(k in ncomp){
            .zz <- fast.z(xx.train, yy.train)

            .theta.k <- svd.theta(.svd, k, .zz)

            yy.hat <- xx.test %*% .theta.k

            corr <- sapply(1:ncol(yy.train), function(j){
                safe.corr(yy.test[, j], yy.hat[, j],
                          use = "pairwise.complete.obs",
                          ...)
            })

            pve <- sapply(1:ncol(yy.train), function(j){
                safe.pve(yy.test[, j], yy.hat[, j])
            })

            .dt <- data.table(corr = corr, pve = pve,
                              r = r, svd.k = k, y.col = 1:ncol(yy))
            cv.result <- rbind(cv.result, .dt)
        }
    }
    return(cv.result)
}

cv.summary <- function(cv.result){
    ret <- cv.result[, .(rr = mean(`corr`, na.rm = T),
                         rr.se = sd(`corr`, na.rm = T) / sqrt(.N),
                         pve = mean(`pve`, na.rm = T),
                         pve.se = sd(`pve`, na.rm = T) / sqrt(.N)),
                     by = .(svd.k, y.col)]
    ret[order(`rr`, decreasing = TRUE),
        head(.SD, 1),
        by = .(y.col)]
}

################################################################

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", `start`, "-", `stop`)]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

.svd.covar <- readRDS(svd.file)
.cmd <- "tabix " %&% expr.file %&% " " %&% .query %&% " " %&% " -h"
expr.dt <- fread(cmd=.cmd, sep="\t", header=T)
nPC <- min(nPC, ncol(.svd.covar$u))
COVAR <- .svd.covar$u[, 1:nPC, drop = FALSE]

message("read all the data")

genes <- unique(expr.dt$hgnc_symbol)

output <- data.table()

for(gene in genes){

    .temp <- expr.dt[hgnc_symbol == gene]

    y.ct <- .temp$celltype
    Y <- as.matrix(t(.temp[, -(1:6)]))
    colnames(Y) <- y.ct
    Y <- .quantile.norm(Y)

    observed <- apply(!is.na(Y), 2, mean)
    if(sum(observed >= .10) < 1) next

    Y <- Y[, observed >= .10, drop = F]

    covar <- apply(COVAR[rownames(Y), , drop = FALSE], 2, scale)

    .lm <- .safe.lm(Y, covar)

    .data <- match.with.plink(.lm$residuals, plink)
    .data$x <- apply(.data$x, 2, scale)

    xx <- .data$x
    xx[is.na(xx)] <- 0
    yy <- .data$y

    cv.result <- data.table()
    for(r in 1:35){
        set.seed(r)
        .temp <- cv.svd.comp(xx, yy, nfold = 3, method = "spearman")
        cv.result <- rbind(cv.result, .temp)
    }

    y.dt <- data.table(celltype=colnames(.data$y),
                       y.col=1:ncol(.data$y))

    .summary <-
        cv.summary(cv.result) %>%
        left_join(y.dt, by = "y.col") %>%
        dplyr::mutate(ld = ld.index, gene) %>%
        dplyr::select(-`y.col`)

    message("Computed: ", gene)

    output <- rbind(output, .summary)
}

if(nrow(output) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
} else {

    out.dt <-
        output %>%
        dplyr::select(`ld`, `gene`, `celltype`,
                      `svd.k`, `rr`, `rr.se`, `pve`, `pve.se`) %>%
        dplyr::mutate(nPC = nPC) %>%
        as.data.table()

    fwrite(out.dt, file=out.file, sep="\t")
}
unlink(temp.dir, recursive=TRUE)
