argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## geno.hdr <- "result/step4/rosmap"
## expr.file.1 <- "result/step3/qc/Mic/Mic_PC50_AD.bed.gz"
## expr.file.2 <- "result/step3/qc/Mic/Mic_PC50_noAD.bed.gz"
## out.file <- "output.txt.gz"

CIS.DIST <- 5e5 # Max distance between SNPs and a gene
PIP.CUTOFF <- 0

if(length(argv) < 6) q()
ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
expr.file.1 <- argv[4]
expr.file.2 <- argv[5]
out.file <- argv[6]

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

    return(ret)
}

run.susie.match <- function(.data.1, .data.2, .suff){

    susie.dt <- data.table()

    xx1 <- apply(.data.1$x, 2, scale)
    xx2 <- apply(.data.2$x, 2, scale)
    xx1[is.na(xx1)] <- 0
    xx2[is.na(xx2)] <- 0

    yy1 <- .data.1$y
    yy2 <- .data.2$y
    m <- ncol(xx1)

    for(k in 1:ncol(yy1)){

        .y1 <- yy1[, k, drop = FALSE]
        .y2 <- yy2[, k, drop = FALSE]
        if(sum(is.finite(.y1)) < 10 || sum(is.finite(.y2)) < 10) next

        .susie.1 <- susie(xx1, .y1, refine=TRUE, na.rm=TRUE)
        .susie.2 <- susie(xx2, .y2, refine=TRUE, na.rm=TRUE)
        .cs1 <- susie_get_cs(.susie.1, coverage=.95)
        .cs2 <- susie_get_cs(.susie.2, coverage=.95)
        .pip1 <- susie_get_pip(.susie.1)
        .pip2 <- susie_get_pip(.susie.2)
        .theta1 <- susie_get_posterior_mean(.susie.1)
        .theta2 <- susie_get_posterior_mean(.susie.2)
        .sd1 <- susie_get_posterior_sd(.susie.1)
        .sd2 <- susie_get_posterior_sd(.susie.2)
        .joint <- rep(0, ncol(xx1))
        max.overlap <- 0
        for(l1 in names(.cs1$cs)){
            for(l2 in names(.cs2$cs)){
                s1 <- .cs1$cs[[l1]]
                s2 <- .cs2$cs[[l2]]
                s12 <- intersect(s1, s2)
                if(length(s12) > 0){
                    p12 <- .pip1[s12] * .pip2[s12]
                    max.overlap <- max(max.overlap, sum(p12))
                    .joint[s12] <- pmax(.joint[s12], p12)
                }
            }
        }

        ## Genetic correlation between two conditions
        X <- apply(rbind(xx1, xx2), 2, scale)
        y1.hat <- X %*% .theta1
        y2.hat <- X %*% .theta2
        .test <- cor.test(y1.hat, y2.hat, use="pairwise.complete.obs")
        if(is.na(.test$estimate)){
            .test$estimate <- 0
            .test$p.value <- 1
        }
        .dt.k <- list()
        .dt.k[[paste0("theta.", .suff[1])]] <- .theta1[1:m]
        .dt.k[[paste0("theta.", .suff[2])]] <- .theta2[1:m]
        .dt.k[[paste0("theta.sd.", .suff[1])]] <- .sd1[1:m]
        .dt.k[[paste0("theta.sd.", .suff[2])]] <- .sd2[1:m]
        .dt.k[[paste0("pip.", .suff[1])]] <- .pip1[1:m]
        .dt.k[[paste0("pip.", .suff[2])]] <- .pip2[1:m]
        .dt.k[["joint"]] <- .joint[1:m]
        .dt.k[["n.overlap"]] <- rep(max.overlap, m)
        .dt.k[["r"]] <- rep(.test$estimate, m)
        .dt.k[["r.pval"]] <- rep(.test$p.value, m)
        .dt.k[["x.col"]] <- 1:m
        .dt.k[["y.col"]] <- rep(k, m)

        susie.dt <- rbind(susie.dt, setDT(.dt.k))
        rm(.susie.1); rm(.susie.2); gc()
        message("Done: ", k)
    }
    return(susie.dt)
}

################################################################

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", pmax(`start` - CIS.DIST, 0), "-", `stop` + CIS.DIST)]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

expr.1 <- fread(cmd = paste0("tabix ", expr.file.1, " ", .query, " -h"))
expr.2 <- fread(cmd = paste0("tabix ", expr.file.2, " ", .query, " -h"))

if(nrow(expr.dt) < 1) {
    unlink(out.file)
    fwrite(data.frame(), file=out.file, sep="\t")
    q()
}

stopifnot(all(expr.1$hgnc_symbol == expr.2$hgnc_symbol))

take.matched.data <- function(expr.dt, plink){

    Y <-
        expr.dt[, 6:ncol(expr.dt)] %>%
        as.matrix() %>%
        t()

    gene.info <- expr.dt[, 1:5]
    colnames(Y) <- gene.info$hgnc_symbol

    match.with.plink(Y, plink)
}

.data.1 <- take.matched.data(expr.1, plink)
.data.2 <- take.matched.data(expr.2, plink)

message("Successfully parsed data: X, Y")

parse.cond <- function(.file.name){
    ret <- tail(strsplit(.file.name, split="[_]")[[1]], 1)
    gsub(".bed.gz","",ret)
}
.suff <- c(parse.cond(expr.file.1),
           parse.cond(expr.file.2))

marginal.1 <- take.marginal.stat(.data.1$x, .data.1$y)
marginal.2 <- take.marginal.stat(.data.2$x, .data.2$y)

marginal.dt <- left_join(marginal.1, marginal.2,
                         by = c("x.col","y.col"),
                         suffix = paste0(".",.suff))

message("Done: Marginal QTL calling")

susie.dt <- run.susie.match(.data.1, .data.2, .suff)

message("Done: SuSie Estimation")

gene.info <- expr.1[, 1:5]
gene.info[, y.col := 1:.N]
snp.info <- plink$BIM
snp.info[, x.col := 1:.N]

out.dt <-
    snp.info %>%
    right_join(marginal.dt) %>%
    left_join(susie.dt) %>%
    left_join(gene.info) %>%
    select(`#chromosome_name`, `snp.loc`, `plink.a1`, `plink.a2`,
           `tss`, `tes`, `hgnc_symbol`,
           starts_with("beta"), starts_with("se"), starts_with("n"),
           starts_with("p.val"), starts_with("theta"), starts_with("pip"),
           `joint`, `r`, `r.pval`) %>%
    mutate(LD = ld.index) %>%
    as.data.table()

fwrite(out.dt, file=out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
