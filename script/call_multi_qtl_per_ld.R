argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## geno.hdr <- "result/step4/rosmap"
## expr.file <- "result/step3/log_mean.bed.gz"
## svd.file <- "result/step3/svd.rds"
## max.K <- 37
## out.file <- "output.txt.gz"

if(length(argv) < 7) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
expr.file <- argv[4]
svd.file <- argv[5]
max.K <- argv[6]
out.file <- argv[7]

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

################################################################

temp.dir <- paste0(out.file, "_temp")

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

    plink.cis <- crop.plink.cis(plink, tss, tes, cis.dist)

    if(ncol(plink.cis$bed) < 1 || nrow(plink.cis$map) < 1) next

    .data <- match.with.plink(.lm$residuals, plink.cis)
    .data$x <- apply(.data$x, 2, scale)

    xx <- .data$x
    xx[is.na(xx)] <- 0
    yy <- .data$y

    x.dt <- as.data.table(plink.cis$map)
    x.dt <- x.dt[, .(`chromosome`,
                     `physical.pos`,
                     `allele1`,
                     `allele2`)] %>%
        cbind(x.col=1:ncol(.data$x))

    y.dt <- data.table(celltype=colnames(.data$y),
                       y.col=1:ncol(.data$y))

    ##################
    ## fine-mapping ##
    ##################

    ## Intersection with the multi-trait fine-mapping
    mtsusie <- mtSusie::mt_susie(X = .data$x,
                                 Y = .data$y,
                                 L = 30,
                                 clamp = 4,
                                 tol = 1e-4,
                                 prior.var = .01,
                                 coverage = .95,
                                 update.prior = T,
                                 local.residual = T)

    susie.dt <- setDT(mtsusie$cs)

    if(nrow(susie.dt) < 1) next

    .temp <- 
        as.data.frame(susie.dt) %>% 
        dplyr::rename(x.col = variants) %>%
        dplyr::rename(y.col = traits) %>%
        dplyr::left_join(x.dt, by="x.col") %>%
        dplyr::left_join(y.dt, by="y.col") %>%
        dplyr::select(-y.col, -x.col) %>%
        dplyr::mutate(gene = g) %>%
        na.omit() %>%
        as.data.table()

    .out <- .temp[order(- abs(z)),
                  head(.SD, 1),
                  by = .(physical.pos, celltype)]
    
    message("Computed: ", g)

    output <- rbind(output, .out)
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