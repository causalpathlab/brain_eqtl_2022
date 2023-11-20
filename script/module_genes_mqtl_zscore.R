argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## mqtl.file <- "result/step4/combined/mqtl_PC37.vcf.gz"
## out.file <- "temp2.txt.gz"

mqtl.file <- argv[1]
out.file <- argv[2]

dir.create(dirname(out.file), recursive = T, showWarnings = F)

library(data.table)
library(dplyr)

`%&%` <- function(a,b) paste0(a,b)

mqtl.dt <- fread(mqtl.file)

## Assign each variant to a unique level to remove redundancy
message("selecting unique gene-level combinations")
.vars <- mqtl.dt[order(p.val, -lbf), head(.SD, 1), by = .(physical.pos, gene, celltype)][, .(physical.pos, levels, gene)] %>% unique()

mqtl.uniq <- left_join(.vars, mqtl.dt) %>% as.data.table()

## Take the lead SNP by p-value and break tie by PIP alpha
message("characterize each gene and level by lead SNP")
lead.snps <- mqtl.uniq[order(p.val, -abs(z)), head(.SD, 1), by = .(celltype, gene, levels)]
                       
.temp <- lead.snps %>% 
    dcast(`#chromosome` +
          `physical.pos` +
          `gene` +
          `levels` ~ `celltype`,
          value.var = "z",
          fun.aggregate = mean,
          fill = 0)

message("dcast by z-scores")

fwrite(.temp, out.file, col.names = T, row.names = F)
