options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.index <- 1609
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## qtl.dir <- "result/step4/qtl/PC37/"
## gwas.stat.file <- "data/gwas/AD.vcf.gz"
## gwas.pgs.dir <- "result/step4/gwas/AD/"
## out.file <- "output.txt.gz"

if(length(argv) < 7) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
qtl.dir <- argv[4]
gwas.stat.file <- argv[5]
gwas.pgs.dir <- argv[6]
out.file <- argv[7]

lfsr.cutoff <- .05
top.geno.PC <- 10

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

fast.z <- function (x, y)
{
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0), replace(y, is.na(y), 0))/sqrt(n.obs)
    ret[is.na(ret)] <- 0
    return(ret)
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

svd.pgs <- function(.svd, k, zz){
    lambda <- 1/nrow(.svd$u)
    udinv.k <- sweep(.svd$u[, 1:k, drop = F], 2, 1/.svd$d[1:k] + lambda, `*`)
    vz.k <- t(.svd$v[, 1:k, drop = F]) %*% zz
    udinv.k %*% vz.k
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
               stwas.num = .num,
               stwas.denom = .denom,
               stwas.z = .num/.denom)
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

if(nrow(qtl.dt) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
    message("Empty QTL statistics after Q/C")
    q()
}

qtl.mu.dt <- qtl.dt %>%
    filter(lfsr < lfsr.cutoff) %>%
    mutate(mu = `mean` * `alpha`) %>%
    mutate(beta.z = `beta` / `se`)

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

gwas.pgs.file <- gwas.pgs.dir %&% "/" %&% ld.index %&% ".pgs.gz"
gwas.pgs.dt <- fread(gwas.pgs.file)

gwas.pgs.raw <-
    gwas.pgs.dt[match(plink$fam$sample.ID, `iid`), .(y)] %>%
    as.matrix()

message("Read GWAS PGS results")

## knn.qtl <- FNN::get.knn(qtl.cf.pgs, 1)
## gwas.pgs.cf <- gwas.pgs.raw[unlist(knn.qtl$nn.index), , drop = F]
## message("Remove QTL PC(s) from GWAS PGS")

.svd <- rsvd::rsvd(safe.scale(plink$bed), k = top.geno.PC)
knn.geno <- FNN::get.knn(safe.scale(.svd$u), 1)
gwas.pgs.cf <- gwas.pgs.raw[unlist(knn.geno$nn.index), , drop = F]
gwas.pgs <- .safe.lm(Y = gwas.pgs.raw, C = gwas.pgs.cf)$residuals

message("Remove top genotype PC(s) from GWAS PGS")

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

message("PGS-based TWAS calculation")

gwas.dt <-
    read.gwas(gwas.stat.file, .query) %>%
    mutate(physical.pos = position) %>%
    as.data.table()

message("Read GWAS summary statistics")

#####################################
## conventional summary-based TWAS ##
#####################################

create.twas.dt <- function(qtl.mu, gwas.dt, plink){

    ## match with plink
    merged <-
        left_join(plink$map, gwas.dt) %>%
        na.omit() %>% 
        as.data.table()

    if(nrow(merged) < 1) {
        return(data.table())
    }

    merged[effect_allele == allele1, beta.flip := beta]
    merged[effect_allele == allele2, beta.flip := -beta]
    merged[is.na(beta.flip), beta.flip := 0]
    merged[, z := `beta.flip`/`standard_error`]

    zz <- merged[order(abs(`z`), decreasing = T), head(.SD, 1), by = .(physical.pos)]
    zz <- as.data.table(qtl.mu)[, .(physical.pos)] %>%
        left_join(zz)
    zz <- as.matrix(zz[, .(z)])
    zz[is.na(zz)] <- 0

    .match <- match(qtl.mu$physical.pos, plink$map$physical.pos)
    xx <- safe.scale(plink$bed)[, .match, drop = F]
    .mu <- as.matrix(qtl.mu[, -1, drop = F])
    .svd <- rsvd::rsvd(xx/sqrt(nrow(xx)))

    compute.twas(.svd, .mu, zz, pve.cutoff = .9)
}

stwas.stat <- create.twas.dt(qtl.mu, gwas.dt, plink)
stwas.stat[, c("celltype", "gene") := tstrsplit(`col`, split="_")]
stwas.stat[, `col` := NULL]

out.dt <- left_join(twas.stat, stwas.stat)

fwrite(out.dt, out.file, sep = "\t", col.names = T, row.names = F)

message("Done")
