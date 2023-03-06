argv <- commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## geno.hdr <- "result/step4/rosmap"
## expr.dir <- "result/step3/qc/"
## expr.ext <- "allIndv_allPC.bed.gz"
## out.file <- "output.txt.gz"

CIS.DIST <- 1e6      # Max distance between SNPs and a gene
alpha.cutoff <- 1e-2 #PIP.CUTOFF

if(length(argv) < 6) q()
ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
expr.dir <- argv[4] # "result/step3/qc/"
expr.ext <- argv[5] # "allIndv_allPC.bed.gz"
out.file <- argv[6]

temp.dir <- paste0(out.file, "_temp")

################################################################
library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)

source(here::here("script", "mtSusie.R"))

subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

    .error <- function(e) {
        print(e)
        cat("Failed to read plink!\n", file=stderr())
        return(NULL)
    }

    dir.create(temp.dir, recursive=TRUE, showWarnings=FALSE)

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num <- as.integer(gsub(pattern="chr", replacement="", chr))

        plink.cmd <- sprintf("./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --hwe 1e-6 --chr %d --from-bp %d --to-bp %d --out %s",
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

read.gene <- function(g, expr.dt, plink){
    .info <- expr.dt[g, 1:5]
    .data <- expr.dt[g, 6:ncol(expr.dt)]

    yy <- data.table(sample=colnames(.data), value=unlist(.data,use.names=FALSE))
    yy[, c("projid", "celltype") := tstrsplit(`sample`, split="_")]

    Y.dt <-
        dcast(yy, projid ~ celltype, value=`value`, fun.aggregate=mean) %>%
        mutate(projid=as.integer(projid)) %>%
        as.data.table()

    .plink.fam <- plink$fam %>%
        mutate(projid=as.integer(`sample.ID`)) %>%
        mutate(x.col=1:n())

    .match <-
        left_join(Y.dt[, .(projid)], .plink.fam, by="projid") %>%
        dplyr::mutate(y.col=1:n()) %>% 
        dplyr::select(y.col, x.col) %>%
        na.omit()

    Y <- as.matrix(Y.dt[.match$y.col, -1])
    X <- as.matrix(plink$bed[.match$x.col, , drop=FALSE])
    list(y=Y, x=X, info=.info)
}

################################################################

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", pmax(`start` - CIS.DIST, 0), "-", `stop` + CIS.DIST)]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb=ld.info[ld.index]$`start`,
                      plink.ub=ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

## read all the genes and cell types within the LD block
expr.files <- 
    lapply(list.files(expr.dir, full.names=TRUE), list.files,
           pattern=paste0(expr.ext, "$"),
           full.names=TRUE) %>%
    unlist(use.names=FALSE)

expr.dt <- data.table()

for(.file in expr.files){
    .dt <- fread(cmd=paste0("tabix ", .file, " ", .query, " -h"))
    if(nrow(.dt) > 1) {
        if(nrow(expr.dt) < 1) {
            expr.dt <- .dt
        } else {
            .by <- c("#chromosome_name", "tss", "tes", "ensembl_gene_id", "hgnc_symbol")
            expr.dt <- full_join(expr.dt, .dt, by=.by)
        }
    }
    message("Read ", .file)
}


out.dt <- data.table()

for(g in 1:nrow(expr.dt)){

    data <- read.gene(g, expr.dt, plink)

    y.dt <- data.table(celltype=colnames(data$y), y.col=1:ncol(data$y))
    x.dt <- data.table(physical.pos=plink$map$physical.pos, x.col=1:ncol(data$x))
    message("Successfully parsed data [", g,"]: X, Y")

    susie <- fit_mt_susie(X=data$x, Y=data$y, L=30)
    message("Done: mtSusie Estimation [", g, "]")

    lfsr <-
        reshape2::melt(susie$stat$lfsr, value.name="lfsr") %>%
        dplyr::rename(l=Var1, y.col=Var2)

    susie.dt <-
        susie$cs %>% na.omit() %>% as.data.table() %>%
        dplyr::left_join(x.dt, by="x.col") %>%
        dplyr::left_join(y.dt, by="y.col") %>%
        dplyr::left_join(lfsr, by=c("l", "y.col")) %>%
        dplyr::select(-x.col, -y.col) %>%
        as.data.table()

    susie.dt <- cbind(data$info[, .(tss, tes, ensembl_gene_id, hgnc_symbol)],
                      susie.dt)

    susie.dt <-
        full_join(plink$map, susie.dt[alpha > alpha.cutoff], by = "physical.pos") %>%
        na.omit() %>%
        dplyr::select(-genetic.dist,-marker.ID) %>% 
        dplyr::mutate(z = beta/se, p = 2 * pnorm(abs(z), lower.tail=FALSE)) %>% 
        as.data.table()

    out.dt <- rbind(out.dt, susie.dt)

    message("Done [", g, "]")
}

fwrite(out.dt, file=out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
