argv <- commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)

## ld.file <- "data/LD.info.txt"
## ld.index <- 1609
## geno.hdr <- "result/step4/rosmap"
## expr.dir <- "result/step3/qc/"
## expr.ext <- "allIndv_allPC.bed.gz"
## out.file <- "output.txt.gz"

if(length(argv) < 6) q()
ld.file <- argv[1]
ld.index <- as.integer(argv[2])
geno.hdr <- argv[3]
expr.dir <- argv[4] # "result/step3/qc/"
expr.ext <- argv[5] # "allIndv_allPC.bed.gz"
out.file <- argv[6]

CIS.DIST <- 1e6      # Max distance between SNPs and a gene
pop.pc <- 3          # local population PCs
min.size <- 10       # minimum sample size

temp.dir <- paste0(out.file, "_temp")

################################################################
library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)

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
        mutate(x.row=1:n())

    .match <-
        left_join(Y.dt[, .(projid)], .plink.fam, by="projid") %>%
        dplyr::mutate(y.row=1:n()) %>%
        dplyr::select(y.row, x.row) %>%
        na.omit()

    Y <- as.matrix(Y.dt[.match$y.row, -1])
    X <- as.matrix(plink$bed[.match$x.row, , drop=FALSE])
    list(y=Y, x=X, info=.info, match = .match)
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
    message("Successfully parsed data [", g,"]: Y")

    cis.variants <- which(plink$map$physical.pos > (data$info$tss - CIS.DIST) &
                          plink$map$physical.pos < (data$info$tss + CIS.DIST))

    if(length(cis.variants) < 1) next

    x.pos <- plink$map$physical.pos[cis.variants]
    data$x <- data$x[, cis.variants, drop = FALSE]

    x.dt <- data.table(physical.pos=x.pos, variants=1:ncol(data$x))

    message("Enforce cis-distance around TSS: ", length(cis.variants), " variants")

    x0 <- apply(data$x, 2, scale)
    x0[is.na(x0)] <- 0
    x.svd <- rsvd::rsvd(x0, pop.pc)
    message("Local population structures")
    x.knn <- FNN::get.knn(x.svd$u, 1)

    xx <- apply(data$x, 2, scale)
    yy <- apply(data$y, 2, scale)
    remove <- which(apply(is.finite(yy), 2, sum) <= min.size)

    if(length(remove) == ncol(yy)) next

    yy <- yy[, -remove, drop = FALSE]
    y.dt <- data.table(celltype=colnames(yy), traits=1:ncol(yy))
    y0 <- yy[x.knn$nn.index, , drop = FALSE]

    susie <- mtSusie::mt_susie(xx, yy - y0, L=30, tol=1e-6,
                               output.full.stat=FALSE,
                               min.pip.cutoff = alpha.cutoff)

    susie.dt <- setDT(susie$cs) %>%
        dplyr::left_join(x.dt, by="variants") %>%
        dplyr::left_join(y.dt, by="traits") %>%
        dplyr::select(-traits, -variants) %>%
        na.omit() %>%
        as.data.table()

    susie.dt <-
        cbind(data$info, susie.dt) %>%
        dplyr::select(`#chromosome_name`, `physical.pos`, `levels`,
                      `tss`, `tes`, `hgnc_symbol`, `celltype`,
                      `alpha`, `mean`, `var`, `lbf`, `z`, `lodds`, `lfsr`) %>%
        as.data.table()

    out.dt <- rbind(out.dt, susie.dt)

    message("Done [", g, "]")
}



fwrite(out.dt, file=out.file, sep="\t")
unlink(temp.dir, recursive=TRUE)
