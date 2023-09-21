argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

feature.file <- argv[1] # e.g., feature.file <- "result/step1/features_annotated_GRCh37.txt.gz"

sum.out.file <- argv[2] # e.g., sum.out.file
mean.out.file <- argv[3]

library(data.table)
library(dplyr)

feature.info <- fread(feature.file)

take.projid <- function(.colnames){
    strsplit(.colnames, split = "_") %>%
        lapply(function(x) x[[1]][1]) %>%
        unlist()
}

rds.files <- list.files("result/step3/pb/", pattern = "rds$", full.names = TRUE)

## 1. Figure out the union of row and column names
projid <- NULL
genes <- NULL

for(.file in rds.files){
    .ct <- gsub(".rds$","",basename(.file))
    .data <- readRDS(.file)
    .projid <- take.projid(colnames(.data$sum))
    .genes <- rownames(.data$sum)
    projid <- sort(unique(c(projid, .projid)))
    genes <- sort(unique(c(genes, .genes)))
    rm(.data)
    message(.ct)
}

## 2. Concatenate them all
sum.out <- NULL
mean.out <- NULL

for(.file in rds.files){
    .ct <- gsub(".rds$","",basename(.file))
    .data <- readRDS(.file)

    ## row-wise concatenation
    .projid <- take.projid(colnames(.data$sum))
    .genes <- rownames(.data$sum)
    .cols <- match(projid, .projid)

    .mat <- .data$sum[.genes, .cols]
    colnames(.mat) <- projid
    rownames(.mat) <- paste0(.genes, "_", .ct)
    sum.out <- rbind(sum.out, .mat)

    ## row-wise concatenation
    .data$ln.mu[.data$sum < 1] <- NA # remove missing spots
    .projid <- take.projid(colnames(.data$ln.mu))
    .genes <- rownames(.data$ln.mu)
    .cols <- match(projid, .projid)

    .mat <- .data$ln.mu[.genes, .cols]
    .mat <- round(.mat, 4)
    colnames(.mat) <- projid
    rownames(.mat) <- paste0(.genes, "_", .ct)
    mean.out <- rbind(mean.out, .mat)

    message(.ct)
    rm(.data)
}

.sort.rows <- function(.mat, feature.info){

    rows <- data.table(row = rownames(.mat))
    rows[, c("ensembl_gene_id", "hgnc_symbol", "celltype") := tstrsplit(`row`, split="_")]
    rows <- rows %>%
        mutate(r = 1:n()) %>%
        left_join(feature.info) %>%
        filter(`chromosome_name` %in% 1:22) %>%
        filter(nchar(`hgnc_symbol`) > 0) %>%
        as.data.table()

    ret <- rows %>%
        select(`chromosome_name`, `tss`, `tes`,
               `ensembl_gene_id`, `hgnc_symbol`, `celltype`) %>%
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

write.bed.gz(.sort.rows(sum.out, feature.info), sum.out.file)
write.bed.gz(.sort.rows(mean.out, feature.info), mean.out.file)
