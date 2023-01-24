options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 89
## geno.hdr <- "result/step4/rosmap"
## gwas.dir <- "result/step4/gwas"
## eqtl.dir <- "result/step4/combined/PC50/all/"
## out.file <- "output.txt.gz"

PIP.CUTOFF <- 1e-8

if(length(argv) < 6) q()
ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
gwas.dir <- argv[4]
eqtl.dir <- argv[5]
out.file <- argv[6]

temp.dir <- paste0(out.file, "_temp")

################################################################

library(data.table)
library(dplyr)
library(reshape2)
library(fqtl)
library(susieR)
`%&%` <- function(a,b) paste0(a,b)

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

fast.cov <- function (x, y) 
{
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0))/n.obs
    return(ret)
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

read.gwas <- function(gwas.file, .query){

    .ret.1 <- fread(cmd="tabix " %&% gwas.file %&% " chr" %&% .query %&% " -h", sep = "\t")
    .ret <- fread(cmd="tabix " %&% gwas.file %&% " " %&% .query %&% " -h", sep = "\t")

    print(gwas.file)

    .ret <- rbind(.ret, .ret.1) %>% 
        as.data.table()

    if("#CHR" %in% colnames(.ret)){
        .ret[, `#chr` := `#CHR`]
    }

    if("chr" %in% colnames(.ret)){
        .ret[, `#chr` := `chr`]
    }

    if("lodds" %in% colnames(.ret)){
        .ret[, beta := lodds]
    }

    .trait <- gsub(".bed.gz$", "", basename(gwas.file))
    .trait <- gsub(".vcf.gz$", "", .trait)
    .ret[, trait := .trait]
    .ret[, `#chr` := sapply(`#chr`, gsub, pattern="chr", replacement="")]
    .ret[, `#chr` := as.integer(`#chr`)]
    return(.ret)
}

read.eqtl <- function(eqtl.file, .query,
                      pval.cutoff = 0.05,
                      pip.cutoff = 0){

    .cmd <- "tabix " %&% eqtl.file %&% " " %&% .query %&% " -h"
    .ret <- fread(cmd=.cmd, sep = "\t")
    if(nrow(.ret) < 1) return(data.table())
    .valid <- .ret[p.val < pval.cutoff &
                   pip > pip.cutoff,
                   .(hgnc_symbol)]

    .celltype <- gsub(".bed.gz$", "", basename(eqtl.file))
    .celltype <- gsub(".vcf.gz$", "", .celltype)

    .ret[hgnc_symbol %in% .valid$hgnc_symbol, ]
    if(nrow(.ret) < 1) return(data.table())
    .ret[, celltype := .celltype]
    .ret
}

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", `start`, "-", `stop`)]

.query <- ld.info[ld.index, ]$query

eqtl.files <- list.files(eqtl.dir,
                         pattern=".vcf.gz$",
                         full.names=TRUE)

eqtl.dt <- do.call(rbind, 
                   lapply(eqtl.files,
                          read.eqtl,
                          .query = .query,
                          pval.cutoff=1e-2,
                          pip.cutoff=PIP.CUTOFF))

if(nrow(eqtl.dt) < 1) {
    fwrite(data.table(), file=out.file, sep="\t")
    q()
}

.plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                       plink.lb = ld.info[ld.index]$`start`,
                       plink.ub = ld.info[ld.index]$`stop`,
                       temp.dir)

unlink(temp.dir, recursive=TRUE)

gwas.files <- list.files(gwas.dir,
                         pattern=".vcf.gz$",
                         full.names=TRUE)

gwas.dt <- do.call(rbind, lapply(gwas.files, read.gwas, .query=.query))

gwas.snps <- gwas.dt %>%
    (function(x)
        .plink$BIM %>%
        mutate(x.col = 1:n()) %>% 
        left_join(x) %>% 
        na.omit()) %>% 
    (function(x) x[, .(snp.loc, x.col, plink.a1, plink.a2)] ) %>%
    unique() %>%
    arrange(snp.loc) %>%
    as.data.table()

.dcast.snps <- function(.dt, f, value.var, fill=0) {
    dcast(.dt, f, value.var=value.var, fun.aggregate = mean, fill = fill) %>%
        (function(x) left_join(gwas.snps[, .(snp.loc)], x)) %>%
        as.data.table()
}

.match.cols <- function(x, y){
    .cols <- colnames(x)
    y[, ...cols]
}

gwas.theta <- .dcast.snps(gwas.dt, snp.loc ~ trait, value.var = "theta")
gwas.se <- .dcast.snps(gwas.dt, snp.loc ~ trait, value.var = "theta.sd", fill=1)
gwas.pip <- .dcast.snps(gwas.dt, snp.loc ~ trait, value.var = "pip")
gwas.se <- .match.cols(gwas.theta, gwas.se)
gwas.pip <- .match.cols(gwas.theta, gwas.pip)

stopifnot(all(colnames(gwas.theta) == colnames(gwas.se)))
stopifnot(all(colnames(gwas.theta) == colnames(gwas.pip)))

eqtl.beta <- .dcast.snps(eqtl.dt, snp.loc ~ hgnc_symbol + celltype,
                         value.var = "beta", fill=0)

eqtl.theta <- .dcast.snps(eqtl.dt, snp.loc ~ hgnc_symbol + celltype,
                          value.var = "theta", fill=0)

eqtl.se <- .dcast.snps(eqtl.dt, snp.loc ~ hgnc_symbol + celltype,
                       value.var = "se", fill=0)

eqtl.pip <- .dcast.snps(eqtl.dt, snp.loc ~ hgnc_symbol + celltype,
                        value.var = "pip", fill=0)

stopifnot(all(gwas.theta$snp.loc == gwas.se$snp.loc))
stopifnot(all(eqtl.theta$snp.loc == gwas.se$snp.loc))
stopifnot(all(eqtl.pip$snp.loc == gwas.se$snp.loc))

take.matrix <- function(x, remove.cols, fill = 0) {
    ret <- as.matrix(x)[, -remove.cols, drop = FALSE]
    ret[is.na(ret)] <- fill
    return(ret)
}

############################
## Predict GWAS and genes ##
############################
X <- scale(.plink$BED[, gwas.snps$x.col, drop = FALSE])

Y.gwas <- scale(X %*% take.matrix(gwas.theta, 1, 0))
X.eqtl <- scale(X %*% take.matrix(eqtl.theta, 1, 0))

.y.name <- data.table(y.col = 1:ncol(Y.gwas), trait = colnames(Y.gwas))
.x.name <- data.table(x.col = 1:ncol(X.eqtl), gene = colnames(X.eqtl))

###################################################
## marginal correlation between genes and traits ##
###################################################

twas.dt <-
    take.marginal.stat(X.eqtl, Y.gwas) %>% 
    left_join(.x.name, by = "x.col") %>% 
    left_join(.y.name, by = "y.col") %>%
    select(trait, gene, x.col, y.col, beta, n, se, p.val)

######################################################
## Perform colocalization between eQTL and GWAS PIP ##
######################################################

dt.1 <- gwas.dt[pip > PIP.CUTOFF, .(snp.loc, pip, trait, k)]
dt.2 <- eqtl.dt[pip > PIP.CUTOFF, .(snp.loc, pip, hgnc_symbol, celltype, k)]

dt.joint <-
    full_join(dt.1, dt.2, by = "snp.loc", suffix = c(".gwas",".eqtl")) %>%
    na.omit() %>%
    as.data.table()

.temp <- dt.joint[, .(pip.overlap = sum(pip.gwas * pip.eqtl)),
                  by = .(hgnc_symbol, celltype, k.eqtl, trait, k.gwas)]

pip.dt <- .temp[order(pip.overlap, decreasing = TRUE),
                head(.SD, 1),
                by = .(trait, hgnc_symbol, celltype)]

out.dt <- twas.dt %>% 
    tidyr::separate(`gene`, c("hgnc_symbol","celltype"), sep="[_]") %>% 
    select(-x.col, -y.col) %>% 
    left_join(pip.dt) %>% 
    filter(`n` > 0) %>% 
    select(-`k.eqtl`, -`k.gwas`) %>% 
    as.data.table()

fwrite(out.dt, out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
