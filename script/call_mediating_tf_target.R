argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

tf.target.ct.file <- "result/step5/tf_target_celltype.txt.gz"
tf.name <- "MEF2A"
geno.hdr <- "result/step4/rosmap"
expr.file <- "result/step5/expr/Mic.bed.gz"
pca.file <- "result/step5/expr/Mic.pca.rds"
out.hdr <- "output"

if(length(argv) != 6) q()

tf.target.ct.file <- argv[1] # e.g., "result/step5/tf_target_celltype.txt.gz"
tf <- argv[2]                # e.g., "MEF2A"
geno.hdr <- argv[3]          # e.g., "result/step4/rosmap"
expr.file <- argv[4]         # e.g., "result/step5/expr/Mic.bed.gz"
pca.file <- argv[5]          # e.g., "result/step5/expr/Mic.pca.rds"
out.hdr <- argv[6]           # e.g., "output.txt.gz"

#############################
## some default parameters ##
#############################

CIS.DIST <- 5e5   # (a) Max distance between SNPs and a gene
PVE.CUTOFF <- 0.8 # (b) proportion of variance explained in PCA
LFSR.CUTOFF <- .1 # (b) To test whether there is a genetic effect
PIP.CUTOFF <- 0   #     confounder vs. collider

dir.create(dirname(out.hdr), recursive=TRUE, showWarnings=FALSE)

temp.dir <- paste0(out.hdr, "_temp")

################################################################

ct.name <- gsub(".bed.gz","",basename(expr.file))
ct.check <- gsub(".pca.rds","",basename(pca.file))
stopifnot(ct.name == ct.check)

message("Cell type: ", ct.name)

out.stat.file <- paste0(out.hdr, ".stat.gz")
out.data.file <- paste0(out.hdr, ".data.gz")

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

run.susie <- function(X, Y, .coverage = .9){

    Y <- apply(Y, 2, scale)
    X <- apply(X, 2, scale)
    X[is.na(X)] <- 0

    susie.dt <- data.table()

    for(k in 1:ncol(Y)){
        yy.k <- Y[,k,drop=FALSE]
        if(sum(is.finite(yy.k)) < 10) next
        .susie.k <- susie(X, yy.k,
                          L = 15,
                          estimate_residual_variance = FALSE,
                          na.rm = TRUE,
                          refine = TRUE)
        .cs <- susie_get_cs(.susie.k, coverage = .coverage)
        m <- ncol(X)
        .factor <- apply(.susie.k$alpha, 2, which.max)[1:m]
        .lfsr <- susie_get_lfsr(.susie.k)
        susie.dt.k <-
            data.table(theta = susie_get_posterior_mean(.susie.k)[1:m],
                       theta.sd = susie_get_posterior_sd(.susie.k)[1:m],
                       k = .factor,
                       pip = susie_get_pip(.susie.k)[1:m],
                       ncs = length(.cs$cs),
                       lfsr = .lfsr[.factor],
                       x.col = 1:m,
                       y.col = k)
        susie.dt <- rbind(susie.dt, susie.dt.k)
        rm(.susie.k); gc()
        ## message("Done: ", k)
    }
    return(susie.dt)
}

###########################
## 0. Take all the lists ##
###########################

tf.target <- fread(tf.target.ct.file)
tf.target <- tf.target[upstream != target]

targets <- tf.target[upstream == tf.name & celltype == ct.name,
                     .(target, tgt.chr, tgt.tss, tgt.tes)]

other.tf <- tf.target[upstream != tf.name & celltype == ct.name,
                      .(upstream, up.chr, up.tss, up.tes)] %>%
    unique()

this.tf <- tf.target[upstream == tf.name & celltype == ct.name,
                     .(upstream, up.chr, up.tss, up.tes)] %>%
    unique()

#####################################################################
## Take genotype data in the cis-regulatory region of the upstream ##
## regulator                                                       ##
#####################################################################

plink.trans <- subset.plink(geno.hdr,
                            chr = this.tf$up.chr,
                            plink.lb = pmax(this.tf$up.tss - CIS.DIST, 0),
                            plink.ub = this.tf$up.tes + CIS.DIST,
                            temp.dir)

unlink(temp.dir, recursive=TRUE)

## 1a. Distinguish between confounder vs. collider factors

.svd <- readRDS(pca.file)
.r <- .svd$d^2
.pve <- cumsum(.r)/sum(.r)
k <- max(1, max(which(.pve <= PVE.CUTOFF)))

svd.uu <- .svd$u[, 1:k, drop = FALSE]
svd.susie <- match.with.plink(svd.uu, plink.trans) %>%
    (function(.data) run.susie(.data$x, .data$y))

collider.factors <-
    unique(svd.susie[lfsr < LFSR.CUTOFF & pip > PIP.CUTOFF, ]$y.col)

svd.conf <- svd.uu

if(length(collider.factors) > 0){
    svd.conf <- svd.uu[, -collider.factors, drop = FALSE]
    message("eliminated ", length(collider.factors), " collider factors")
}

if(ncol(svd.conf) < 1){
    svd.conf <- matrix(0, nrow(svd.conf), 1)
}

rownames(svd.conf) <- rownames(svd.uu)

message("Established confounders from the PCA results")

## Combine other transcription factor genes

other.tf[, query := paste0(`up.chr`, ":", `up.tss`, "-", `up.tes`)]

other.tf.dt <- data.table()

for(k in 1:nrow(other.tf)){
    .tf.k <- other.tf[k, ]$upstream
    .query.k <- other.tf[k, ]$query
    .expr.k <- fread(cmd = paste0("tabix ", expr.file, " ", .query.k, " -h"))
    .expr.k <- .expr.k[hgnc_symbol == .tf.k]
    if(nrow(.expr.k) < 0) next
    other.tf.dt <- rbind(other.tf.dt, .expr.k)
}

## 1b. Distinguish between confounder vs. collider TFs

U <- 
    other.tf.dt[, 6:ncol(other.tf.dt)] %>%
    as.matrix() %>%
    t()

other.info <- other.tf.dt[, 1:5]
colnames(U) <- other.info$hgnc_symbol

other.tf.susie <- match.with.plink(U, plink.trans) %>%
    (function(.data) run.susie(.data$x, .data$y))

collider.genes <-
    unique(other.tf.susie[lfsr < LFSR.CUTOFF & pip > PIP.CUTOFF, ]$y.col)

if(length(collider.genes) > 0){
    conf.dt <- other.tf.dt[- collider.genes, 6:ncol(other.tf.dt)]
    conf.info <- other.tf.dt[- collider.genes, 1:5]
    message("eliminated ", length(collider.genes), " collider genes")
} else {
    conf.dt <- other.tf.dt[, 6:ncol(other.tf.dt)]
    conf.info <- other.tf.dt[, 1:5]
}

if(nrow(conf.dt) < 1) {
    uu <- matrix(0, nrow=ncol(other.tf.dt) - 5, ncol=1)
} else if (nrow(conf.dt) > 10) {
    U <- scale(t(as.matrix(conf.dt)))
    U[is.na(U)] <- 0
    .svd <- rsvd::rsvd(U, k = min(10, ncol(U)))
    uu <- .svd$u
} else {
    uu <- scale(t(as.matrix(conf.dt)))
    uu[is.na(uu)] <- 0
}

rownames(uu) <- colnames(conf.dt)
u.order <- match(rownames(svd.conf), rownames(uu))
uu <- cbind(uu, svd.conf[u.order, , drop=FALSE])

message("Found total ", ncol(uu), " putative confounder variables")

## 1c. Test cis-eQTLs for this TF

this.tf[, query := paste0(`up.chr`, ":", `up.tss`, "-", `up.tes`)]
this.tf.dt <- fread(cmd = paste0("tabix ", expr.file, " ", this.tf$query, " -h"))
this.tf.dt <- this.tf.dt[hgnc_symbol == this.tf$upstream]

M0 <- scale(t(as.matrix(this.tf.dt[, 6:ncol(this.tf.dt)])))
M <- .safe.lm(M0, uu)$residuals
rownames(M) <- rownames(M0)

stat.cis <- match.with.plink(M, plink.trans) %>% 
    (function(.data) {
        take.marginal.stat(.data$x, .data$y) %>% 
            left_join(run.susie(.data$x, .data$y),
                      by = c("x.col","y.col"))
    }) %>%
    select(-y.col)

#########################################################################
## 2. Visit for each target gene with this FT and putative confounders ##
#########################################################################

targets[, query := paste0(`tgt.chr`, ":", `tgt.tss`, "-", `tgt.tes`)]

snp.info <- plink.trans$BIM
snp.info[, x.col := 1:.N]

out.stat <- data.table()
out.data <- data.table()

for(k in 1:nrow(targets)){

    tgt.info <- targets[k, ]

    tgt.query <- tgt.info$query
    tgt.dt <- fread(cmd = paste0("tabix ", expr.file, " ", tgt.query, " -h"))
    tgt.dt <- tgt.dt[hgnc_symbol == tgt.info$target]

    if(nrow(tgt.dt) < 1) next # empty expression

    ## Can we adjust some confounding TFs?
    Y0 <- scale(t(as.matrix(tgt.dt[, 6:ncol(tgt.dt)])))

    ## 2a. Regress out confounder effects
    Y <- .safe.lm(Y0, uu)$residuals
    rownames(Y) <- rownames(Y0)

    .data.k <- data.table(target=as.numeric(Y),
                          tf=as.numeric(M),
                          target0=as.numeric(Y0),
                          tf0=as.numeric(M0),
                          tf.name = tf.name,
                          target.name = tgt.info$target)

    ## 2b. Test a trans-eQTL model
    stat.trans <- match.with.plink(Y, plink.trans) %>% 
        (function(.data) {
            take.marginal.stat(.data$x, .data$y) %>% 
                left_join(run.susie(.data$x, .data$y),
                          by = c("x.col","y.col"))
        }) %>%
        (function(.dt) .dt[, x.col := 1:.N])

    if(max(stat.trans$ncs) < 1) next # empty trans-eQTL discovery

    .stat.k <-
        snp.info[, .(chr,snp.loc,plink.a1,plink.a2,x.col)] %>%
        left_join(stat.trans, by = "x.col") %>%
        select(-y.col) %>%
        mutate(tf.name = tf.name,
               target.name = tgt.info$target)

    out.stat <- rbind(out.stat, as.data.table(.stat.k))
    out.data <- rbind(out.data, as.data.table(.data.k))

    message("Done: " , k , " / ", nrow(targets))
}

out.stat <-
    out.stat %>%
    left_join(stat.cis, by = "x.col", suffix=c(".trans",".cis")) %>%
    mutate(celltype = ct.name) %>% 
    select(-x.col)

out.data <-
    out.data %>% mutate(celltype = ct.name)

fwrite(out.stat, file=out.stat.file, sep="\t", row.names = FALSE)
fwrite(out.data, file=out.data.file, sep="\t", row.names = FALSE)

message("Finished trans-eQTL calling")
unlink(temp.dir, recursive=TRUE)
