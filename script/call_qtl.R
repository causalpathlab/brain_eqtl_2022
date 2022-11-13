argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 3
## geno.hdr <- "result/step4/rosmap"
## expr.file <- "result/step3/qc/Mic/Mic_AD_all.bed.gz"
## out.file <- "output.txt.gz"

CIS.DIST <- 5e5 # Max distance between SNPs and a gene
PIP.CUTOFF <- 0

if(length(argv) < 5) q()
ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
expr.file <- argv[4]
out.file <- argv[5]

temp.dir <- paste0(out.file, "_temp")

################################################################
library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)
library(fqtl)
library(susieR)

subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {
    require(fqtl)
    require(dplyr)
    .error <- function(e) {
        print(e)
        cat("Failed to read plink!\n", file=stderr())
        return(NULL)
    }

    dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num <- as.integer(gsub(pattern = "chr", replacement = "", chr))
        plink.cmd <- sprintf("./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --hwe 1e-6 --chr %d --from-bp %d --to-bp %d --out %s",
                             plink.hdr, chr.num, plink.lb, plink.ub, paste0(temp.dir, "/plink"))
        system(plink.cmd)

        plink <- read.plink(paste0(temp.dir, "/plink"))
        colnames(plink$BIM) <- c("chr", "rs", "missing", "snp.loc", "plink.a1", "plink.a2")
        colnames(plink$FAM) <- c("fam", "iid", "father", "mother", "sex.code", ".pheno")

        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 <- "T"
        }

        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 <- "T"
        }
        return(plink)
    }
    plink <- tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error = .error)
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
    list(fitted = Y.fitted, residuals = Y.resid)
}

take.marginal.stat <- function(xx, yy, se.min=1e-8) {
    .xx <- scale(xx)
    .yy <- scale(yy)
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

match.with.plink <- function(Y, plink){

    .y.info <- data.table(y.name = rownames(Y))
    .y.info[, y.row := 1:.N]

    ny <- strsplit(.y.info$`y.name`, split="_")[[1]]
    if(length(ny) > 2){
        .y.info[, c("projid", "matched", "ct") := tstrsplit(`y.name`, split="_")]
        .y.info[, projid := as.integer(projid)]
        .y.info[, matched := as.integer(matched)]
        .fam <- plink$FAM[, .(iid)] %>% mutate(x.row = 1:n())
        .x.1 <- .fam %>% rename(projid = iid)
        .x.2 <- .fam %>% rename(matched = iid)
        .match <-
            left_join(.y.info, .x.1, by = "projid") %>%
            left_join(.x.2, by = "matched", suffix=c("",".matched")) %>%
            na.omit()

        X <- plink$BED

        x1 <- X[.match$x.row, , drop = FALSE]
        x0 <- X[.match$x.row.matched, , drop = FALSE]
        ret <- list(x = x1 - x0,
                    y = Y[.match$y.row, , drop = FALSE])

    } else {
        .y.info[, c("projid", "ct") := tstrsplit(`y.name`, split="_")]
        .y.info[, projid := as.integer(projid)]
        ## match individuals
        .match <-
            plink$FAM %>%
            mutate(projid = as.integer(iid)) %>%
            mutate(x.row = 1:n()) %>%
            left_join(.y.info, by = "projid") %>%
            na.omit()

        X <- plink$BED

        ret <- list(x = X[.match$x.row, , drop = FALSE],
                    y = Y[.match$y.row, , drop = FALSE])
    }

    ############################
    ## More Q/C to remove PC1 ##
    ############################

    xx <- apply(ret$x, 2, scale)
    xx[is.na(xx)] <- 0
    pc1 <- rsvd::rsvd(xx, k=1)
    ret$y <- .safe.lm(ret$y, pc1$u)$residuals

    return(ret)
}

run.susie <- function(X, Y){

    susie.dt <- data.table()

    for(k in 1:ncol(Y)){
        yy.k <- Y[,k,drop=FALSE]
        .susie.k <- susie(X, yy.k,
                          L = 5,
                          estimate_residual_variance = FALSE,
                          coverage = .9,
                          na.rm = TRUE,
                          refine = TRUE)
        susie.dt.k <-
            data.table(theta = susie_get_posterior_mean(.susie.k),
                       theta.sd = susie_get_posterior_sd(.susie.k),
                       pip = susie_get_pip(.susie.k),
                       lfsr = min(susie_get_lfsr(.susie.k)))
        susie.dt.k[, x.col := 1:.N]
        susie.dt.k[, y.col := k]
        rm(.susie.k); gc()
        susie.dt <- rbind(susie.dt, susie.dt.k)
    }
    return(susie.dt)
}

################################################################

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", pmax(`start` - CIS.DIST, 0), "-", `stop` + - CIS.DIST)]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)


expr.dt <- fread(cmd = paste0("tabix ", expr.file, " ", .query, " -h"))

if(nrow(expr.dt) < 1) {
    unlink(out.file)
    fwrite(data.frame(), file=out.file, sep="\t")
    q()
}

Y <-
    expr.dt[, 6:ncol(expr.dt)] %>%
    as.matrix() %>%
    t()

gene.info <- expr.dt[, 1:5]
colnames(Y) <- gene.info$hgnc_symbol

.data <- match.with.plink(Y, plink)

message("Successfully parsed data: X, Y")

marginal.dt <- take.marginal.stat(.data$x, .data$y)

message("Done: Marginal QTL calling")

Y <- apply(.data$y, 2, scale)
X <- apply(.data$x, 2, scale)
X[is.na(X)] <- 0

susie.dt <- run.susie(X, Y)

message("Done: SuSie Estimation")

gene.info[, y.col := 1:.N]
snp.info <- plink$BIM
snp.info[, x.col := 1:.N]

out.dt <-
    snp.info %>%
    right_join(marginal.dt) %>%
    left_join(susie.dt) %>%
    left_join(gene.info) %>%
    filter(pip > PIP.CUTOFF) %>%
    select(`#chromosome_name`, `snp.loc`, `plink.a1`, `plink.a2`, `tss`, `tes`, `ensembl_gene_id`, `hgnc_symbol`, `beta`, `se`, `n`, `p.val`, `theta`, `theta.sd`, `pip`, `lfsr`) %>%
    mutate(LD = ld.index) %>% 
    as.data.table() %>% 
    na.omit()

fwrite(out.dt, file=out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
