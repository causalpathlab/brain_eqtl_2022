options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.index <- 1609
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## qtl.dir <- "result/step4/qtl/PC37/"
## gwas.pgs.dir <- "result/step4/gwas/AD/"
## out.file <- "output.txt.gz"

if(length(argv) < 6) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
qtl.dir <- argv[4]
gwas.pgs.dir <- argv[5]
out.file <- argv[6]

lfsr.cutoff <- .05

################################################################

temp.dir <- paste0(out.file, "_temp")

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

safe.scale <- function(.mat){
    ret <- apply(.mat, 2, scale)
    ret[is.na(ret)] <- 0
    ret
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
        dplyr::mutate(y.col = as.integer(y.col)) %>%
        as.data.table()

    return(out)
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

compute.twas <- function(.svd, .mu, zz, pve.cutoff = .9){
    pve <- cumsum(.svd$d^2) / sum(.svd$d^2)
    k <- min(min(which(pve > pve.cutoff)), length(.svd$d))

    lambda <- 1/nrow(.svd$u)
    vd <- sweep(.svd$v[, 1:k, drop = F], 2, .svd$d[1:k], `*`)

    mu.vd <- t(.mu) %*% vd

    .num <- as.numeric(t(.mu) %*% zz)
    .denom <- apply(mu.vd, 1, function(x) sqrt(sum(x^2) + lambda))
    data.table(col = colnames(.mu),
               twas.num = .num,
               twas.denom = .denom,
               twas.z = .num/.denom)
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

qtl.file <- qtl.dir %&% "/" %&% ld.index %&% ".txt.gz"
qtl.dt <- fread(qtl.file)

qtl.mu.dt <- qtl.dt %>%
    filter(lfsr < lfsr.cutoff) %>%
    mutate(mu = `mean` * `alpha`)

if(nrow(qtl.mu.dt) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
    message("Empty QTL statistics after Q/C")
    q()
}

message("Read QTL results")

qtl.mu <-
    dcast(qtl.mu.dt,
          physical.pos ~ celltype + gene,
          fill = 0,
          value.var = "mu")

.match <- match(qtl.mu$physical.pos, plink$map$physical.pos)
xx <- safe.scale(plink$bed)[, .match, drop = F]
qtl.pgs <- xx %*% as.matrix(qtl.mu[, -1])
colnames(qtl.pgs) <- colnames(qtl.mu)[-1]
rownames(qtl.pgs) <- plink$fam$sample.ID

message("Estimated polygenic scores derived from QTL")

max.K <- 150
x.svd <- rsvd::rsvd(xx, k=min(max.K, ncol(xx)))

gwas.cv.file <- gwas.pgs.dir %&% "/" %&% ld.index %&% ".cv.gz"
gwas.cv.dt <- fread(gwas.cv.file)

max.K <- which.max(gwas.cv.dt$rr.mean)
geno.pc <- x.svd$u[, 1:min(ncol(x.svd$u), max.K), drop = F]

message("Remove ", ncol(geno.pc), " population PC(s)")

qtl.pgs <- .safe.lm(qtl.pgs, geno.pc)$residuals

message("Adjusted population structures in TWAS")

gwas.pgs.file <- gwas.pgs.dir %&% "/" %&% ld.index %&% ".pgs.gz"
gwas.pgs.dt <- fread(gwas.pgs.file)

gwas.pgs <-
    gwas.pgs.dt[match(plink$fam$sample.ID, `iid`), .(y)] %>%
    as.matrix()

message("Read GWAS PGS results")

col.dt <-
    data.table(col = colnames(qtl.pgs)) %>%
    mutate(x.col = 1:n()) %>%
    as.data.table()

col.dt[, c("celltype", "gene") := tstrsplit(`col`, split="_")]
col.dt[, col:= NULL]

twas.stat <-
    take.marginal.stat(qtl.pgs, gwas.pgs) %>%
    left_join(col.dt) %>%
    select(-x.col, -y.col) %>%
    mutate(ld = ld.index) %>%
    select(`ld`, `gene`, `celltype`, `beta`, `se`, `n`, `p.val`) %>%
    as.data.table()

fwrite(twas.stat, out.file, sep = "\t", col.names = T, row.names = F)

message("Done")
