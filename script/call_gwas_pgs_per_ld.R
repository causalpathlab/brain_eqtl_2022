argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.index <- 1609
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## gwas.file <- "data/gwas/AD.vcf.gz"
## out.file <- "output.pgs.gz"

if(length(argv) < 5) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
gwas.file <- argv[4]
out.file <- argv[5]

cv.file <- gsub(".pgs.gz", ".cv.gz", out.file)
temp.dir <- paste0(out.file, "_temp")
max.K <- 100
nrepeat <- 20
nfold <- 5

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
    vd2inv <- sweep(.svd$v[, 1:k, drop = F], 2, 1/.svd$d[1:k] + lambda, `*`)
    vz.k <- t(.svd$v[, 1:k, drop = F]) %*% zz
    .theta.k <- vd2inv %*% vz.k
}

svd.pgs <- function(.svd, k, zz){
    lambda <- 1/nrow(.svd$u)
    udinv.k <- sweep(.svd$u[, 1:k, drop = F], 2, 1/.svd$d[1:k] + lambda, `*`)
    vz.k <- t(.svd$v[, 1:k, drop = F]) %*% zz
    udinv.k %*% vz.k
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

safe.scale <- function(.mat){
    ret <- apply(.mat, 2, scale)
    ret[is.na(ret)] <- 0
    ret
}

take.svd.sub <- function(.svd, .sub){
    ret <- .svd
    ret$u <- ret$u[.sub, , drop = F]
    return(ret)
}

################################################################

read.gwas <- function(gwas.file, .query){

    .ret.1 <- fread(cmd="tabix " %&% gwas.file %&% " chr" %&% .query %&% " -h", sep = "\t")
    .ret <- fread(cmd="tabix " %&% gwas.file %&% " " %&% .query %&% " -h", sep = "\t")

    print(gwas.file)

    .ret <- rbind(.ret, .ret.1) %>%
        as.data.table()

    .trait <- gsub(".vcf.gz$", "", basename(gwas.file))
    .ret[, trait := .trait]
    return(.ret[, .(position, variant_id, p_value, effect_allele, other_allele,
                    beta, standard_error)])
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

gwas.dt <- read.gwas(gwas.file, .query)

if(nrow(gwas.dt) < 1) {
    fwrite(data.table(), file=out.file, sep="\t")
    q()
}

## match with plink
merged <- plink$map %>%
    mutate(plink.pos = 1:n()) %>%
    as.data.table() %>%
    merge(gwas.dt, by.x = "physical.pos", by.y = "position") %>%
    na.omit() %>%
    as.data.table()

if(nrow(merged) < 1) {
    fwrite(data.table(), file=out.file, sep="\t")
    q()
}

merged[effect_allele == allele1, beta.flip := beta]
merged[effect_allele == allele2, beta.flip := -beta]
merged[is.na(beta.flip), beta.flip := 0]

zz <- matrix(merged$beta.flip / merged$standard_error, ncol = 1)
zz[zz > 8] <- 8
zz[zz < -8] <- -8

X <- safe.scale(plink$bed[, merged$plink.pos, drop = F])

## Full polygenic risk prediction
max.K <- min(min(dim(X)), max.K)

## Tune SVD PCs
nn <- nrow(X)
cv.dt <- data.table()

for(r in 1:nrepeat){

    set.seed(r)

    .idx <- sample(nfold, nn, replace=TRUE)

    for(k in 1:nfold){

        x.train <- X[.idx != k, , drop = F]
        x.test <- X[.idx == k, , drop = F]

        svd.train <- rsvd::rsvd(x.train/sqrt(nrow(x.train)), k = max.K)
        svd.test <- rsvd::rsvd(x.test/sqrt(nrow(x.test)), k = max.K)
        message("computed SVD")

        n.train <- nrow(x.train)
        n.test <- nrow(x.test)

        p.train <- sqrt(n.train/(n.train + n.test))
        p.test <- sqrt(n.test/(n.train + n.test))

        for(k in 1:max.K){
            .theta <- svd.theta(svd.train, k, zz * p.train)
            y.hat <- .quantile.norm(x.test %*% .theta)
            y.test <- .quantile.norm(svd.pgs(svd.test, k, zz * p.test))
            rr <- as.numeric(safe.corr(y.hat, y.test))
            cv.dt <- rbind(cv.dt, data.table(k, rr))
        }
    }
}

cv.summary <- cv.dt[, .(rr.mean = mean(rr), rr.sd = sd(rr)), by = .(k)]

k.opt <- unlist(cv.summary[which.max(`rr.mean`), .(k)])[1]

.svd <- rsvd::rsvd(X/sqrt(nn), k = max.K)
y.opt <- as.numeric(svd.pgs(.svd, k.opt, zz))

out.dt <-
    data.table(ld = ld.index,
               iid = plink$fam$sample.ID,
               y = y.opt)

fwrite(out.dt, out.file, sep = "\t", row.names = F, col.names = T)

cv.out.dt <- cbind(ld = ld.index, cv.summary)
fwrite(cv.out.dt, cv.file, sep = "\t", row.names = F, col.names = T)
