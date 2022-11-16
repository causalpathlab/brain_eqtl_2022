argv <- commandArgs(trailingOnly = TRUE)

annot.file <- argv[1]
out.file <- argv[2]

library(data.table)
library(dplyr)

annot.dt <- fread(annot.file, sep="\t")

.rename.map <- fread("data/inh_merging.tsv", sep = "\t", header=FALSE,
                     col.names = c("celltype.eQTLs","update"))
.rename.map[, update := paste0("Inh-", `update`)]

annot.dt <- annot.dt[, c(2, 3, 11, 14)] %>%
    left_join(.rename.map) %>%
    as.data.table()

annot.dt[!is.na(`update`), celltype.eQTLs := `update`]

annot.dt[, barcode := stringr::str_remove(`X`, paste0(`batch`, "_"))]
annot.dt[, barcode := stringr::str_remove(`barcode`, "-[0-9]+$")]
annot.dt[, batch := stringr::str_remove(`batch`, "SM_")]
annot.dt[, celltype := gsub("[()]","",`celltype.eQTLs`)]
annot.dt[, celltype := gsub("[ /]","-",`celltype`)]

out.dt <- annot.dt[, .(barcode, batch, projid, celltype)]
fwrite(out.dt, out.file, row.names = FALSE, col.names = TRUE, sep="\t")
