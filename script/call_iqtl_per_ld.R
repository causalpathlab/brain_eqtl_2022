argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 207
## geno.hdr <- "result/step4/rosmap"
## expr.file <- "result/step3/log_mean.bed.gz"
## herit.dir <- "result/step4/heritability/"
## svd.file <- "result/step3/svd.rds"
## pheno.file <- "data/metadata_PFC_all_individuals_092520.tsv.gz"
## max.K <- 37
## out.file <- "output.txt.gz"

if(length(argv) < 9) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
expr.file <- argv[4]
herit.dir <- argv[5]
svd.file <- argv[6]
pheno.file <- argv[7]
max.K <- argv[8]
out.file <- argv[9]

temp.dir <- paste0(out.file, "_temp")

################################################################

library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)
library(susieR)

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

fast.z <- function (x, y)
{
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0), replace(y, is.na(y), 0))/sqrt(n.obs)
    ret[is.na(ret)] <- 0
    return(ret)
}

svd.theta <- function(.svd, k, zz){
    lambda <- 1/nrow(.svd$u)
    vd2inv <- sweep(.svd$v[, 1:k, drop = F], 2, 1/(.svd$d[1:k]^2 + lambda), `*`)
    vz.k <- t(.svd$v[, 1:k, drop = F]) %*% zz
    .theta.k <- vd2inv %*% vz.k
}

svd.pgs <- function(.svd, k, zz){
    lambda <- 1/nrow(.svd$u)
    udinv.k <- sweep(.svd$u[, 1:k, drop = F], 2, 1/(.svd$d[1:k] + sqrt(lambda)), `*`)
    vz.k <- t(.svd$v[, 1:k, drop = F]) %*% zz
    udinv.k %*% vz.k
}

################################################################

safe.scale <- function(.mat){
    ret <- apply(.mat, 2, scale)
    ret[is.na(ret)] <- 0
    ret
}

################################################################

adj.by.svd <- function(.data){

    imp.r <- function(r){
        ww <- .data$w[,r]

        x0 <- .data$x[ww == 0, , drop = F]
        y0 <- .data$y[ww == 0, , drop = F]
        x1 <- .data$x[ww == 1, , drop = F]
        y1 <- .data$y[ww == 1, , drop = F]

        z1 <- fast.z(x1, y1)
        z0 <- fast.z(x0, y0)

        .svd1 <- rsvd::rsvd(x1, k=max(.data$svd.k))
        .svd0 <- rsvd::rsvd(x0, k=max(.data$svd.k))


        y0.for.w1 <- sapply(1:length(.data$svd.k),
                            function(j) svd.pgs(.svd1,
                                                .data$svd.k[j],
                                                z0[,j]))

        y1.for.w0 <- sapply(1:length(.data$svd.k),
                            function(j) svd.pgs(.svd0,
                                                .data$svd.k[j],
                                                z1[,j]))

        y1.adj <- sapply(1:ncol(y1), function(j){
            .safe.lm(y1[,j,drop=F], y0.for.w1[,j,drop=F])$residuals
        })

        y0.adj <- sapply(1:ncol(y1), function(j){
            .safe.lm(y0[,j,drop=F], y1.for.w0[,j,drop=F])$residuals
        })

        list(y0 = y0.adj, y1 = y1.adj, x0 = x0, x1 = x1)
    }

    ret <- lapply(1:ncol(.data$w), imp.r)
    names(ret) <- colnames(.data$w)
    return(ret)
}

################################################################

cis.dist <- 5e5

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", pmax(`start` - cis.dist, 0), "-", as.integer(`stop` + cis.dist))]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

.svd.covar <- readRDS(svd.file)
.cmd <- "tabix " %&% expr.file %&% " " %&% .query %&% " " %&% " -h"
expr.dt <- fread(cmd=.cmd, sep="\t", header=T)

herit.file <- herit.dir %&% "/PC" %&% max.K %&% "/" %&% ld.index %&% ".txt.gz"
herit.dt <- fread(herit.file)

pheno.dt <- fread(pheno.file, header = T)

## pheno.dt[, var_braaksc := scale(`braaksc`)]
## pheno.dt[, var_amyloid := scale(sqrt(`amyloid`))]
## pheno.dt[, var_nft := scale(sqrt(`nft`))]
## pheno.dt[, var_msex := scale(`msex`)]
## pheno.dt[, var_age := scale(`age_death`)]
## pheno.dt[, var_braaksc := scale(`braaksc`)]

message("read all the data")

genes <- intersect(unique(expr.dt$hgnc_symbol), herit.dt$gene)

output <- data.table()

for(g in genes){

    .temp <- expr.dt[hgnc_symbol == g]

    tss <- min(.temp$tss)
    tes <- max(.temp$tes)

    y.ct <- .temp$celltype
    Y <- as.matrix(t(.temp[, -(1:6)]))
    colnames(Y) <- y.ct
    Y <- .quantile.norm(Y)

    observed <- apply(!is.na(Y), 2, mean)
    if(sum(observed >= .10) < 1) next

    Y <- Y[, observed >= .10, drop = F]

    #################################
    ## Adjust technical covariates ##
    #################################
    covar <- apply(.svd.covar$u[rownames(Y), 1:max.K, drop = F], 2, scale)
    .lm <- .safe.lm(Y, covar)
    .data <- match.with.plink(.lm$residuals, plink)
    .data$x <- safe.scale(.data$x)

    x.dt <- as.data.table(plink$map)
    x.dt <- x.dt[, .(`chromosome`,
                     `physical.pos`,
                     `allele1`,
                     `allele2`)] %>%
        cbind(x.col=1:ncol(.data$x))

    y.dt <- data.table(celltype=colnames(.data$y), y.col=1:ncol(.data$y))

    ################
    ## conditions ##
    ################
    W <- pheno.dt[match(rownames(.data$y), `projid`),
                  .(projid, msex, niareagansc)]

    W[, AD := 0]
    W[`niareagansc` <= 2, AD := 1]
    .data$w <- as.matrix(W[, .(msex, AD)])

    ######################
    ## optimal SVD rank ##
    ######################

    .data$svd.k <- herit.dt[gene == g][match(colnames(.data$y),`celltype`)]$svd.k

    cond.data <- adj.by.svd(.data)

    .out <- data.table()

    for(ii in 1:length(cond.data)){
        .w <- names(cond.data)[ii]
        .dt1 <- take.marginal.stat(cond.data[[ii]]$x1, cond.data[[ii]]$y1)
        .dt0 <- take.marginal.stat(cond.data[[ii]]$x0, cond.data[[ii]]$y0)

        if(nrow(.dt1) > 0){
            .dt1[, cond := .w]
            .dt1[, W := 1]
            .out <- rbind(.out, .dt1)
        }

        if(nrow(.dt0) > 0){
            .dt0[, cond := .w]
            .dt0[, W := 0]
            .out <- rbind(.out, .dt0)
        }
    }

    .out <- .out %>%
        left_join(x.dt) %>% 
        left_join(y.dt) %>% 
        dplyr::select(-y.col, -x.col) %>%
        dplyr::mutate(gene = g) %>%
        dplyr::filter(physical.pos >= (tss - cis.dist)) %>% 
        dplyr::filter(physical.pos <= (tes + cis.dist)) %>% 
        na.omit() %>%
        as.data.table()

    message("Computed: ", g)

    output <- rbind(output, .out)
}

if(nrow(output) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
} else {

    out.dt <-
        output %>%
        dplyr::select(`chromosome`, `physical.pos`, `cond`, `W`,
                      `gene`, `celltype`,
                      `beta`, `se`, `n`, `p.val`) %>%
        arrange(chromosome, physical.pos) %>%
        dplyr::rename(`#chromosome` = `chromosome`) %>%
        as.data.table()

    fwrite(out.dt, file=out.file, sep="\t")
}
unlink(temp.dir, recursive=TRUE)
