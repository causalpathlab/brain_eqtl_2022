argv <- commandArgs(trailingOnly = TRUE)

annot.file <- argv[1]
out.file <- argv[2]

library(data.table)

annot.dt <- fread(annot.file, sep="\t")
annot.dt <- annot.dt[, c(2, 3, 11, 14)]

annot.dt[, barcode := stringr::str_remove(`X`, paste0(`batch`, "_"))]
annot.dt[, barcode := stringr::str_remove(`barcode`, "-[0-9]+$")]
annot.dt[, batch := stringr::str_remove(`batch`, "SM_")]
annot.dt[, celltype := gsub("[()]","",`celltype.eQTLs`)]
annot.dt[, celltype := gsub("[ /]","-",`celltype`)]

out.dt <- annot.dt[, .(barcode, batch, projid, celltype)]
fwrite(out.dt, out.file, row.names = FALSE, col.names = TRUE, sep="\t")
