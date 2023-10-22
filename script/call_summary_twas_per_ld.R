options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.index <- 1609
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## herit.dir <- "result/step4/heritability/PC37/"
## qtl.dir <- "result/step4/qtl/PC37/"
## gwas.stat.file <- "data/gwas/AD.vcf.gz"
## gwas.pgs.dir <- "result/step4/gwas/AD/"
## out.file <- "output.txt.gz"


lfsr.cutoff <- .05
pv.cutoff <- 1e-8

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

gwas.dt <- read.gwas(gwas.stat.file, .query)

## match with plink
merged <- plink$map %>%
    mutate(plink.pos = 1:n()) %>%
    as.data.table() %>%
    merge(gwas.dt, by.x = "physical.pos", by.y = "position") %>%
    na.omit() %>%
    as.data.table()

if(nrow(merged) < 1) {
    fwrite(data.table(), file=out.file, sep="\t")
    message("Empty GWAS results")
    q()
}

merged[effect_allele == allele1, beta.flip := beta]
merged[effect_allele == allele2, beta.flip := -beta]
merged[is.na(beta.flip), beta.flip := 0]

zz <- matrix(merged$beta.flip / merged$standard_error, ncol = 1)

message("Read GWAS summary statistics")

max.K <- 100

X <- safe.scale(plink$bed)
x.svd <- rsvd::rsvd(X / sqrt(nrow(X)), k = max.K)
x.svd$v <- x.svd$v[merged$plink.pos, , drop = F]

message("Computed SVD for the LD estimation")

.temp <-
    dcast(qtl.mu.dt, physical.pos ~ celltype + gene, fill = 0, value.var = "mu") %>%
    mutate(physical.pos = physical.pos)

qtl.mu <- left_join(merged[, .(physical.pos)], .temp)

.mu <- as.matrix(qtl.mu[, -1])
.mu[is.na(.mu)] <- 0

twas.out <- compute.twas(x.svd, .mu, zz, .99)

herit.file <- herit.dir %&% "/" %&% ld.index %&% ".txt.gz"
herit.dt <- fread(herit.file)
herit.dt[, col := celltype %&% "_" %&% gene]
herit.dt[, pv := pnorm(rr/pmax(rr.se, 1e-8), lower.tail = F)]

message("Read the heritability files")

out.dt <-
    left_join(herit.dt, twas.out) %>%
    na.omit() %>% 
    mutate(twas.pv = 2 * pnorm(abs(twas.z), lower.tail = F)) %>% 
    as.data.table()
out.dt[, c("celltype", "gene") := tstrsplit(`col`, split="_")]
out.dt[, col:= NULL]
