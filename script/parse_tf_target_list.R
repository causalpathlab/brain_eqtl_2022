## in.file <- "data/markers.eTFset.RData"
## feature.file <- "result/step1/features_annotated_GRCh37.txt.gz"
## out.file <- 

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

in.file <- argv[1]      # "data/markers.eTFset.RData"
feature.file <- argv[2] #
out.file <- argv[3]     #
key.file <- argv[4]     #

load(in.file)

library(data.table)

.dt <- data.table()
for(ct in names(eTFset)){
    markers <- names(eTFset[[ct]])
    for(tf in names(eTFset[[ct]])){
        targets <- eTFset[[ct]][[tf]]
        .dt.tf <- data.table(upstream=tf, target=targets, celltype=ct)
        .dt <- rbind(.dt, .dt.tf)
    }
}

## Annotate upstream and target genes
feat.dt <- fread(feature.file)

up.idx <- match(.dt$upstream, feat.dt$hgnc_symbol)
tgt.idx <- match(.dt$target, feat.dt$hgnc_symbol)
up.loc <- feat.dt[up.idx, .(chromosome_name, tss, tes)]
tgt.loc <- feat.dt[tgt.idx, .(chromosome_name, tss, tes)]

out.dt <- data.table(upstream = .dt$upstream,
                     up.chr = up.loc$chromosome_name,
                     up.tss = up.loc$tss,
                     up.tes = up.loc$tes,
                     target = .dt$target,
                     tgt.chr = tgt.loc$chromosome_name,
                     tgt.tss = tgt.loc$tss,
                     tgt.tes = tgt.loc$tes,
                     celltype = .dt$celltype)

fwrite(out.dt,
       out.file,
       sep="\t",
       row.names=FALSE)

key.dt <- unique(out.dt[,.(upstream, celltype)])
fwrite(key.dt,
       key.file,
       sep="\t",
       row.names=FALSE,
       col.names=FALSE)

