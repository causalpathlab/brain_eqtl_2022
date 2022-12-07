argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

if(length(argv) < 3) q()

expr.file <- argv[1]    # e.g., expr.file = "result/step3/pb/Mic.rds"
feature.file <- argv[2] # e.g., feature.file = "result/step1/features_annotated_GRCh37.txt.gz"
out.file <- argv[3]

library(data.table)
library(dplyr)

.mkdir <- function(...) {
    dir.create(..., recursive = TRUE, showWarnings = FALSE)
}

`%&%` <- function(a,b) paste0(a,b)

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
        ret[.pos.k] <- x.k.valid
    }
    return(ret)
}

.sort.cols <- function(.mat, pheno) {
    .dt <- data.table(col = colnames(.mat))
    .dt[, c("projid", "ct") := tstrsplit(`col`, split="_")]
    .loc <- match(pheno$projid, as.integer(.dt$projid))
    .mat[, .loc, drop = FALSE]
}

.sort.rows <- function(.mat, feature.info){

    rows <- data.table(row = rownames(.mat))
    rows[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`row`, split="_")]
    rows <- rows %>%
        mutate(r = 1:n()) %>%
        left_join(feature.info) %>%
        filter(`chromosome_name` %in% 1:22) %>%
        as.data.table()

    ret <- rows %>%
        select(`chromosome_name`, `tss`, `tes`,
               `ensembl_gene_id`, `hgnc_symbol`) %>%
        cbind(.mat[rows$r, , drop = FALSE]) %>%
        arrange(`chromosome_name`, `tss`) %>%
        rename(`#chromosome_name` = `chromosome_name`) %>% 
        as.data.table()

    return(ret)
}

write.bed.gz <- function(.dt, out.file){
    if(!file.exists(out.file)){
        out.raw.file <- gsub(".gz$", "", out.file)
        fwrite(.dt, out.raw.file, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
        ## Rsamtools::bgzip(out.raw.file, dest=out.file)
        ## Rsamtools::indexTabix(out.file, format="vcf")
        system(paste0("bgzip ", out.raw.file))
        system(paste0("tabix -p bed ", out.file))
        unlink(out.raw.file)
        message("Wrote: ", out.file)
    }
}

## Remove genes with too few expression values
## - we would consider zero expression is missing
filter.mat <- function(.mat, .sum, missing.cutoff = 1, missing.rate = 0.8) {
    .row.missing.rate <- apply(t(.sum) < missing.cutoff, 2, mean)
    .qc <- .mat
    .qc[.sum < missing.cutoff] <- NA
    .qc[.row.missing.rate < missing.rate, , drop = FALSE]
}

.mkdir(dirname(out.file))

feature.info <- fread(feature.file)
expr <- readRDS(expr.file)

message("Read expression data")

mu.qc <- filter.mat(expr$PB$mu, expr$PB$sum)

mu.mat <- 
    .quantile.norm(mu.qc) %>% 
    .sort.cols(pheno = expr$pheno)

mu.dt <- mu.mat %>% 
    .sort.rows(feature.info)

write.bed.gz(mu.dt, out.file)

message("Average expression matrix")
