options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 657
## geno.hdr <- "result/step4/rosmap"
## gwas.file <- "data/gwas/bentham_sle.bed.gz"
## out.file <- "output.txt.gz"

if(length(argv) < 5) q()

ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
gwas.file <- argv[4]
out.file <- argv[5]

temp.dir <- paste0(out.file, "_temp")

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

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", `start`, "-", `stop`)]

.query <- ld.info[ld.index, ]$query

.plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                       plink.lb = ld.info[ld.index]$`start`,
                       plink.ub = ld.info[ld.index]$`stop`,
                       temp.dir)

unlink(temp.dir, recursive=TRUE)

gwas.dt <- read.gwas(gwas.file, .query)

gwas.flip <- .plink$BIM %>%
    mutate(x.col = 1:n()) %>%
    left_join(gwas.dt[, .(snp.loc, a1, a2, trait)]) %>%
    na.omit %>%
    as.data.table()

gwas.flip[a1 == plink.a1 & a2 == plink.a2, flip := 1]
gwas.flip[a2 == plink.a1 & a1 == plink.a2, flip := -1]

gwas.dt <- gwas.dt %>%
    left_join(gwas.flip) %>%
    na.omit() %>%
    mutate(beta.qc = beta * flip) %>%
    as.data.table()

## handle very small p-value and SE
p.val.min <- 1e-100
z.abs.max <- qnorm(p.val.min * 2, lower.tail = FALSE)
gwas.dt[`p` < p.val.min & `se` <= 0, `se` := abs(`beta`)/z.abs.max]

gwas.snps <-
    gwas.flip[, .(snp.loc, x.col, plink.a1, plink.a2)] %>%
    unique() %>%
    arrange(snp.loc) %>%
    as.data.table()

if(nrow(gwas.snps) < 1){
    fwrite(data.table(), out.file)
    unlink(temp.dir, recursive=TRUE)
    q()
}

.dcast.snps <- function(.dt, f, value.var, fill=0) {
    dcast(.dt, f, value.var=value.var, fun.aggregate = mean, fill = fill) %>%
    (function(x) left_join(gwas.snps[, .(snp.loc)], x)) %>%
    as.data.table()
}

.match.cols <- function(x, y){
    .cols <- colnames(x)
    y[, ...cols]
}

gwas.beta <- .dcast.snps(gwas.dt, snp.loc ~ trait, value.var = "beta.qc")
gwas.se <- .dcast.snps(gwas.dt, snp.loc ~ trait, value.var = "se", fill=1)
## gwas.pval <- .dcast.snps(gwas.dt, snp.loc ~ trait, value.var = "p")
gwas.se <- .match.cols(gwas.beta, gwas.se)

stopifnot(all(colnames(gwas.beta) == colnames(gwas.se)))
stopifnot(all(gwas.beta$snp.loc == gwas.se$snp.loc))

take.matrix <- function(x, remove.cols, fill = 0) {
    ret <- as.matrix(x)[, -remove.cols, drop = FALSE]
    ret[is.na(ret)] <- fill
    return(ret)
}

X <- .plink$BED[, gwas.snps$x.col]
R <- fast.cov(X, X)

n.gwas <- nrow(.plink$BED) ## Scale down GWAS to eQTL samples...
Z.gwas <- take.matrix(gwas.beta, 1, 0)/take.matrix(gwas.se, 1, 1)

## avoid extreme z-score
z.lb <- max(-20, quantile(Z.gwas, .01, na.rm=TRUE))
z.ub <- min(20, quantile(Z.gwas, .99, na.rm=TRUE))

Z.gwas[Z.gwas < z.lb] <- z.lb
Z.gwas[Z.gwas > z.ub] <- z.ub

.susie <- susie_rss(Z.gwas, R,
                    L=25,
                    max_iter=200,
                    verbose=TRUE,
                    track_fit=FALSE,
                    check_input=FALSE,
                    check_prior=FALSE,
                    refine=FALSE) # refinement was too slow

message("Successfully finished SUSIE")

.cs <- susie_get_cs(.susie, coverage = 0.9)
m <- ncol(X)
.factor <- apply(.susie$alpha, 2, which.max)[1:m]
.lfsr <- susie_get_lfsr(.susie)

susie.dt <-
    data.table(theta = susie_get_posterior_mean(.susie)[1:m],
               theta.sd = susie_get_posterior_sd(.susie)[1:m],
               pip = susie_get_pip(.susie)[1:m],
               ncs = length(.cs$cs),
               coverage = sum(.cs$coverage),
               k = .factor,
               lfsr = .lfsr[.factor],
               x.col = gwas.snps$x.col)

out.dt <- gwas.dt %>%
    select(`#chr`, snp.loc, trait) %>% 
    left_join(gwas.snps) %>% 
    left_join(susie.dt) %>%
    na.omit() %>%
    select(-x.col) %>%
    as.data.table()

fwrite(out.dt, file=out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
