argv <- commandArgs(trailingOnly = TRUE)
row.file <- argv[1]
out.file <- argv[2]

library(data.table)
library(dplyr)

row.dt <- fread(row.file, header=FALSE, col.names = "V1")
row.dt[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`V1`,split="_")]
row.dt[, V1 := NULL]

ensembl <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
                            host="uswest.ensembl.org",
                            dataset="hsapiens_gene_ensembl",
                            GRCh = 37)

ensembl.hs <- biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)

.attr <- c("ensembl_gene_id",
           "hgnc_symbol",
           "chromosome_name",
           "transcription_start_site",
           "transcript_start",
           "transcript_end",
           "description",
           "percentage_gene_gc_content")

.matched <- biomaRt::getBM(attributes=.attr,
                           filters="ensembl_gene_id",
                           values=row.dt$ensembl_gene_id,
                           mart=ensembl.hs,
                           useCache = FALSE) %>% 
    as.data.table()

.matched <-
    .matched[, .(tss = min(transcript_start),
                 tes = max(transcript_end),
                 gc = mean(percentage_gene_gc_content)),
             by = .(ensembl_gene_id, hgnc_symbol, chromosome_name)]

fwrite(.matched, file = out.file)
