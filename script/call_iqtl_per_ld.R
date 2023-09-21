argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## geno.hdr <- "result/step4/rosmap"
## expr.file <- "result/step3/log_mean.bed.gz"
## svd.file <- "result/step3/svd.rds"
## pheno.file <- "data/metadata_PFC_all_individuals_092520.tsv.gz"
## out.file <- "output.txt.gz"

if(length(argv) < 7) q()

ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
expr.file <- argv[4]
svd.file <- argv[5]
pheno.file <- argv[6]
out.file <- argv[7]

temp.dir <- paste0(out.file, "_temp")

################################################################

library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)
library(susieR)

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

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", `start`, "-", `stop`)]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

.svd.covar <- readRDS(svd.file)
.cmd <- "tabix " %&% expr.file %&% " " %&% .query %&% " " %&% " -h"
expr.dt <- fread(cmd=.cmd, sep="\t", header=T)

pheno.dt <- fread(pheno.file, header = T)
pheno.dt[, var_amyloid := scale(sqrt(`amyloid`))]
pheno.dt[, var_nft := scale(sqrt(`nft`))]
pheno.dt[, var_msex := scale(`msex`)]
pheno.dt[, var_age := scale(`age_death`)]
pheno.dt[, var_braaksc := scale(`braaksc`)]

message("read all the data")

genes <- unique(expr.dt$hgnc_symbol)

output <- data.table()

for(gene in genes){

    .temp <- expr.dt[hgnc_symbol == gene]

    y.ct <- .temp$celltype
    Y <- as.matrix(t(.temp[, -(1:6)]))
    colnames(Y) <- y.ct
    Y <- .quantile.norm(Y)

    observed <- apply(!is.na(Y), 2, mean)
    if(sum(observed >= .10) < 1) next

    Y <- Y[, observed >= .10, drop = F]

    W <- pheno.dt[match(rownames(Y), `projid`),
                  .(var_amyloid,
                    var_nft,
                    var_msex,
                    var_age,
                    var_braaksc)] %>% 
        as.matrix()

    rownames(W) <- rownames(Y)

    covar <- apply(.svd.covar$u[rownames(Y), ], 2, scale)
    .lm <- .safe.lm(Y, covar)

    .data <- match.with.plink(.lm$residuals, plink)
    .data$x <- apply(.data$x, 2, scale)
    .data$w <- W[match(rownames(.data$y), rownames(W)), , drop = F]

    x.dt <- as.data.table(plink$map)
    x.dt <- x.dt[, .(`chromosome`,
                     `physical.pos`,
                     `allele1`,
                     `allele2`)] %>%
        cbind(x.col=1:ncol(.data$x))

    w.dt <- data.table(pheno = gsub("var_","",colnames(W))) %>%
        mutate(w.col = 1:n())

    y.dt <- data.table(celltype=colnames(.data$y), y.col=1:ncol(.data$y))

    susie <- mtSusie::mt_susie(X = .data$x,
                               Y = .data$y,
                               W = .data$w,
                               L = 1,
                               prior.var = .01,
                               coverage = .9,
                               output.full.stat = F)

    susie.dt <-
        setDT(susie$cs) %>%
        dplyr::rename(x.col = variants) %>%
        dplyr::rename(y.col = traits) %>%
        dplyr::rename(w.col = interaction) %>%
        na.omit() %>%
        filter(lodds > 0, w.col > 0) %>%
        as.data.table()

   susie.best <-
        susie.dt[order(abs(`z`), decreasing = T),
                 head(.SD, 1),
                 by = .(x.col, y.col)]

    .out <- susie.best %>%
        left_join(x.dt) %>% 
        left_join(y.dt) %>% 
        left_join(w.dt) %>% 
        dplyr::select(-y.col, -x.col, -w.col) %>%
        dplyr::mutate(gene) %>%
        na.omit() %>%
        as.data.table()

    message("Computed: ", gene)

    output <- rbind(output, .out)
}

if(nrow(output) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
} else {

    out.dt <-
        output %>%
        dplyr::select(`chromosome`, `physical.pos`, `pheno`,
                      `gene`, `celltype`,
                      `alpha`, `mean`, `var`, `lbf`, `z`, `lodds`, `lfsr`) %>%
        dplyr::mutate(`alpha` = round(`alpha`, 4),
                      `mean` = round(`mean`, 4),
                      `var` = round(`var`, 4),
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
