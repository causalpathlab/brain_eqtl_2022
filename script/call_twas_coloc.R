options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## geno.hdr <- "result/step4/rosmap"
## gwas.dir <- "data/gwas"
## gwas.susie.dir <- "result/step4/gwas"
## eqtl.dir <- "result/step4/combined/allIndv_allPC/"
## out.file <- "output.txt.gz"

PIP.CUTOFF <- 1e-8

if(length(argv) < 7) q()
ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
gwas.dir <- argv[4]
gwas.susie.dir <- argv[5]
eqtl.dir <- argv[6]
out.file <- argv[7]

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
        rename(`snp.loc` = `stop`) %>%
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

    if("pval" %in% colnames(.ret)){
        .ret[, p := pval]
    }

    .trait <- gsub(".bed.gz$", "", basename(gwas.file))
    .ret[, trait := .trait]
    .ret[, `#chr` := sapply(`#chr`, gsub, pattern="chr", replacement="")]
    .ret[, `#chr` := as.integer(`#chr`)]
    .ret[, .(`#chr`, `snp.loc`, `a1`, `a2`, `beta`, `se`, `p`, `trait`)]
}

read.gwas.susie <- function(gwas.file, .query){

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

gwas.susie.files <- list.files(gwas.susie.dir,
                               pattern=".vcf.gz$",
                               full.names=TRUE)

gwas.susie.dt <- do.call(rbind, lapply(gwas.susie.files, read.gwas.susie, .query=.query))

gwas.files <- list.files(gwas.dir,
                         pattern=".bed.gz$",
                         full.names=TRUE)

gwas.dt <-
    do.call(rbind, lapply(gwas.files, read.gwas, .query=.query)) %>%
    left_join(gwas.susie.dt) %>% 
    na.omit()

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

eqtl.theta <- .dcast.snps(eqtl.dt, snp.loc ~ hgnc_symbol + celltype,
                          value.var = "theta", fill=0)

stopifnot(all(eqtl.theta$snp.loc == gwas.theta$snp.loc))

take.matrix <- function(x, remove.cols, fill = 0) {
    ret <- as.matrix(x)[, -remove.cols, drop = FALSE]
    ret[is.na(ret)] <- fill
    return(ret)
}

############################
## Predict GWAS and genes ##
############################
X <- apply(.plink$BED[, gwas.snps$x.col, drop = FALSE], 2, scale)
X[is.na(X)] <- 0

Y.gwas <- apply(X %*% take.matrix(gwas.theta, 1, 0), 2, scale)
X.eqtl <- apply(X %*% take.matrix(eqtl.theta, 1, 0), 2, scale)

.y.name <- data.table(y.col = 1:ncol(Y.gwas), trait = colnames(Y.gwas))
.x.name <- data.table(x.col = 1:ncol(X.eqtl), gene = colnames(X.eqtl))

###################################################
## marginal correlation between genes and traits ##
###################################################

twas.dt <-
    take.marginal.stat(X.eqtl, Y.gwas) %>% 
    left_join(.x.name, by = "x.col") %>% 
    left_join(.y.name, by = "y.col") %>%
    dplyr::rename(twas.beta = beta, twas.n = n, twas.se = se, twas.p = p.val) %>% 
    dplyr::select(trait, gene, dplyr::starts_with("twas")) %>%
    as.data.table()

twas.dt[, c("hgnc_symbol", "celltype") := tstrsplit(`gene`, split="_")]
twas.dt[, gene := NULL]

######################################################
## Perform colocalization between eQTL and GWAS PIP ##
######################################################

dt.1 <- gwas.dt[pip > 1e-8, .(snp.loc, beta, se, pip, trait, k, p, a1, a2)]
dt.2 <- eqtl.dt[pip > 1e-8, .(snp.loc, pip, hgnc_symbol, celltype, k)]

dt.joint <-
    left_join(dt.2, dt.1, by = "snp.loc", suffix = c(".gwas",".eqtl")) %>%
    as.data.table()

calc.joint <- function(p1, p2){
    l1 <- log(p1)
    l2 <- log(p2)
    l12 <- l1 + l2
    l12 <- l12 - max(l12)
    l12 - log(sum(exp(l12)))
}

dt.joint[,
         lbf.joint := calc.joint(pip.gwas, pip.eqtl),
         by = .(hgnc_symbol, celltype, trait, k.eqtl, k.gwas)]

argmax.eqtl <-
    eqtl.dt[order(p.val), head(.SD, 1), by = .(hgnc_symbol, celltype)] %>% 
    dplyr::rename(eqtl.beta=beta, eqtl.se=se, eqtl.n=n, eqtl.p=p.val, eqtl.pip = pip, eqtl.lfsr = lfsr, eqtl.snp = snp.loc) %>% 
    dplyr::select(`#chromosome_name`, dplyr::starts_with("eqtl"), celltype, hgnc_symbol)

argmax.gwas <-
    dt.joint[order(p, decreasing = FALSE), head(.SD, 1), by = .(hgnc_symbol, celltype, trait)] %>% 
    dplyr::rename(gwas.beta=beta, gwas.sd=se, gwas.p = p) %>%
    dplyr::mutate(gwas.snp = snp.loc %&% ":" %&% a1 %&% ":" %&% a2) %>% 
    dplyr::select(celltype, hgnc_symbol, trait, dplyr::starts_with("gwas"))

pip.dt <- dt.joint[,
                   .(snp.joint = snp.loc[which.max(lbf.joint)],
                     pip.overlap = ( sum(pip.eqtl * pip.gwas) /
                                     max(sum(pip.eqtl), sum(pip.gwas)) )),
                  by = .(hgnc_symbol, celltype, trait)]

out.dt <-
    argmax.eqtl %>%
    left_join(twas.dt) %>%
    left_join(pip.dt) %>% 
    left_join(argmax.gwas) %>% 
    dplyr::filter(`twas.n` > 0) %>% 
    na.omit() %>% 
    as.data.table()

fwrite(out.dt, out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
