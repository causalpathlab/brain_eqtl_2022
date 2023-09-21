argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

svd.file <- argv[1] # e.g., "result/step3/svd.rds"
pheno.file <- argv[2] # e.g., "data/metadata_PFC_all_individuals_092520.tsv.gz"
out.file <- argv[3]

library(data.table)
library(tidyverse)

.svd <- readRDS(svd.file)
uu <- .svd$u
ud <- sweep(uu, 2, .svd$d, `*`)

take.marginal.stat <- function(xx, yy, se.min=1e-8) {
    .xx <- apply(as.matrix(xx), 2, scale)
    .yy <- apply(as.matrix(yy), 2, scale)
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


pheno.dt <- fread(pheno.file, header = TRUE)

.idx <- match(rownames(uu), pheno.dt$projid)

.dt <- data.table(pheno.dt[.idx, .(braaksc, age_death, msex, amyloid, tangles, study)])

.dt[, amyloid := sqrt(amyloid)]
.dt[, amyloid := sqrt(tangles)]
.dt[, study := as.integer(factor(`study`, c("ROS","MAP"), 0:1))]

x.col <- data.table(x.col = 1:ncol(.dt), pheno = colnames(.dt))
y.col <- data.table(PC = 1:ncol(uu), ve = .svd$d^2)
y.col[, pve := cumsum(ve)/sum(ve)]

.stat <-
    take.marginal.stat(.dt, uu) %>%
    left_join(x.col) %>%
    select(-x.col) %>% 
    rename(PC = y.col) %>%
    left_join(y.col) %>% 
    as.data.table()

fwrite(.stat, out.file)
