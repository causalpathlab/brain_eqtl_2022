argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 207
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

PV <- 1e-4   # p-value
ALPHA <- 0.1 # PIP cutoff
LFSR <- .2  # local false sign range

################################################################

library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)

source("Util.R")

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

    observed <- apply(!is.na(Y), 2, mean)
    if(sum(observed >= .10) < 1) next

    Y <- Y[, observed >= .10, drop = F]

    covar <- apply(.svd.covar$u[rownames(Y), 1:max.K, drop = FALSE],
                   MARGIN = 2,
                   FUN = scale)

    .lm <- safe.lm(Y, covar)

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

    .data$y <- apply(.data$y, 2, scale)
    .data$x <- apply(.data$x, 2, scale)

    .data$x[!is.finite(.data$x)] <- NA
    .data$y[!is.finite(.data$y)] <- NA

    ## Intersection with the multi-trait fine-mapping
    mtsusie <- mtSusie::mt_susie(X = .data$x,
                                 Y = .data$y,
                                 L = 15,
                                 tol = 1e-8,
                                 prior.var = .01,
                                 coverage = .95)

    susie.dt <- setDT(mtsusie$cs)

    if(nrow(susie.dt) < 1) next

    #########################
    ## marginal statistics ##
    #########################
    marg.stat <-
        take.marginal.stat(.data$x, .data$y)

    .out <- 
        as.data.frame(susie.dt) %>% 
        dplyr::rename(x.col = variant) %>%
        dplyr::rename(y.col = trait) %>%
        dplyr::left_join(x.dt, by="x.col") %>%
        dplyr::left_join(y.dt, by="y.col") %>%
        dplyr::left_join(marg.stat) %>% 
        dplyr::select(-y.col, -x.col) %>%
        dplyr::mutate(gene = g) %>%
        na.omit() %>%
        as.data.table()

    ## keeping everything will be too much
    .valid.pos <-
        .out[p.val < PV | (alpha > ALPHA & lfsr < LFSR),
             .(physical.pos)] %>%
        unlist(use.names = F)

    .out.qc <- .out[physical.pos %in% .valid.pos]

    message("Computed: ", g)

    output <- rbind(output, .out.qc)
}

message("Computed all the genes")

if(nrow(output) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
} else {

    out.dt <-
        output %>%
        dplyr::select(`chromosome`, `physical.pos`, `level`,
                      `gene`, `celltype`,
                      `alpha`, `mean`, `sd`, `lbf`, `z`, `lodds`, `lfsr`,
                      `beta`, `se`, `n`, `p.val`) %>%
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
