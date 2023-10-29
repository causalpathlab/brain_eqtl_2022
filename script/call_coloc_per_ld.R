options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.index <- 207
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## expr.file <- "result/step3/log_mean.bed.gz"
## svd.file <- "result/step3/svd.rds"
## max.K <- 37
## gwas.stat.file <- "data/gwas/AD.vcf.gz"
## out.file <- "output.txt.gz"

if(length(argv) < 8) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
expr.file <- argv[4]
svd.file <- argv[5]
max.K <- argv[6]
gwas.stat.file <- argv[7]
out.file <- argv[8]

################################################################

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

read.gwas <- function(gwas.file, .query){

    .ret.1 <- fread(cmd="tabix " %&% gwas.file %&% " chr" %&% .query %&% " -h", sep = "\t")
    .ret <- fread(cmd="tabix " %&% gwas.file %&% " " %&% .query %&% " -h", sep = "\t")

    print(gwas.file)

    .ret <- rbind(.ret, .ret.1) %>%
        as.data.table()

    .trait <- gsub(".vcf.gz$", "", basename(gwas.file))
    .ret[, trait := .trait]

    .ret[order(p_value),
         head(.SD, 1),
         by = .(position)][,
                           .(position, variant_id,
                             p_value,
                             effect_allele,
                             other_allele,
                             beta, standard_error)]
}

match.expr.plink <- function(Y, plink){

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

match.gwas.plink <- function(gwas.dt, plink){
    ## match with plink
    matched <-
        left_join(plink$map, gwas.dt) %>%
        as.data.table()

    if(nrow(matched) < 1) {
        return(data.table())
    }

    matched[effect_allele == allele1, beta.flip := beta]
    matched[effect_allele == allele2, beta.flip := -beta]
    matched[is.na(beta.flip), beta.flip := 0]
    matched[, z := `beta.flip`/`standard_error`]

    return(matched)
}

crop.plink.cis <- function(plink, tss, tes, cis){
    ret <- plink

    .valid <- which(plink$map$physical.pos >= (tss - cis) &
                    plink$map$physical.pos <= (tes + cis))

    ret$map <- ret$map[.valid,]
    ret$bed <- ret$bed[, .valid, drop = F]

    return(ret)
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

svd.pgs <- function(.svd, zz){
    zz[is.na(zz)] <- 0
    lambda <- 1/nrow(.svd$u)
    vdinv <- sweep(.svd$v, 2, 1/(.svd$d + lambda), `*`)
    .svd$u %*% (t(vdinv) %*% zz)
}

################################################################

temp.dir <- paste0(out.file, "_temp")

################################################################

cis.dist <- 5e5

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", `start`, "-", `stop`)]

.query <- ld.info[ld.index, ]$query

gwas.dt <-
    read.gwas(gwas.stat.file, .query) %>%
    mutate(physical.pos = position) %>%
    as.data.table()

message("Read GWAS summary statistics")

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

.svd.covar <- readRDS(svd.file)
.cmd <- "tabix " %&% expr.file %&% " " %&% .query %&% " " %&% " -h"
expr.dt <- fread(cmd=.cmd, sep="\t", header=T)

message("read all the data")

genes <- unique(expr.dt$hgnc_symbol)

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

    covar <- apply(.svd.covar$u[rownames(Y), 1:max.K, drop = FALSE],
                   MARGIN = 2,
                   FUN = scale)

    .lm <- .safe.lm(Y, covar)

    ## Just focus on the cis window
    plink.cis <- crop.plink.cis(plink, tss, tes, cis.dist)

    .data <- match.expr.plink(.lm$residuals, plink.cis)

    x.dt <- as.data.table(plink.cis$map)
    x.dt <- x.dt[, .(`chromosome`,
                     `physical.pos`,
                     `allele1`,
                     `allele2`)] %>%
        cbind(x.col=1:ncol(.data$x))

    y.dt <- data.table(celltype=colnames(.data$y),
                       y.col=1:ncol(.data$y))

    ## slice GWAS
    .gwas <- match.gwas.plink(gwas.dt, plink.cis)

    kk <- min(500, min(dim(.data$x)))

    .svd.geno <- rsvd::rsvd(safe.scale(.data$x) / sqrt(nrow(.data$x)), k=kk)

    y.gwas <- svd.pgs(.svd.geno, as.matrix(.gwas[, .(z)]))

    coloc.g <- data.table()

    for(j in 1:ncol(.data$y)){

        yy.j <- cbind(.data$y[, j, drop = F], y.gwas)

        .joint <- mtSusie::mt_susie(X = .data$x,
                                    Y = yy.j,
                                    L = 10,
                                    clamp = 8,
                                    tol = 1e-4,
                                    prior.var = .01,
                                    coverage = .95,
                                    update.prior = T,
                                    local.residual = T)

        .temp.j <- setDT(.joint$cs)

        .vars.j <-
            .temp.j[lodds > 0, .(n = .N),
                             by = .(variants, levels)] %>%
            dplyr::filter(`n` > 1) %>%
            dplyr::select(-`n`) %>%
            dplyr::left_join(.temp.j) %>% 
            dplyr::filter(traits == 1) %>% ## just keep the eQTL side
            dplyr::mutate(traits = j) %>%  ## 
            dplyr::rename(x.col = variants) %>%
            dplyr::rename(y.col = traits) %>%
            dplyr::left_join(x.dt, by="x.col") %>%
            dplyr::left_join(y.dt, by="y.col") %>%
            dplyr::select(-y.col, -x.col) %>%
            as.data.table()

        coloc.g <- rbind(coloc.g, .vars.j)
    }

    if(nrow(coloc.g) < 1) next

    coloc.g[, gene := g]

    message("Computed: ", g)

    output <- rbind(output, coloc.g)
}


message("Computed all the genes")

if(nrow(output) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
} else {

    out.dt <-
        output %>%
        dplyr::select(`chromosome`, `physical.pos`, `levels`,
                      `gene`, `celltype`,
                      `alpha`, `mean`, `sd`, `lbf`, `z`, `lodds`, `lfsr`) %>%
        dplyr::mutate(`alpha` = round(`alpha`, 4),
                      `mean` = round(`mean`, 4),
                      `lbf` = round(`lbf`, 4),
                      `z` = round(`z`, 4),
                      `lodds` = round(`lodds`, 4),
                      `lfsr` = round(`lfsr`, 4)) %>% 
        arrange(chromosome, physical.pos) %>%
        dplyr::rename(`#chromosome` = `chromosome`) %>%
        as.data.table()

    fwrite(out.dt, file=out.file, sep="\t")

}

unlink(temp.dir, recursive=TRUE)

message("Done")
