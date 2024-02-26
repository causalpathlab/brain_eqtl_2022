options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.index <- 207
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## qtl.dir <- "result/step4/qtl/PC37/"
## out.file <- "output.txt.gz"

if(length(argv) < 5) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
qtl.dir <- argv[4]
out.file <- argv[5]

################################################################
dir.create(dirname(out.file), recursive = T, showWarnings = F)
temp.dir <- paste0(out.file, "_temp")

library(data.table)
library(dplyr)
library(reshape2)
library(rsvd)


`%&%` <- function(a,b) paste0(a,b)

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

safe.scale <- function(.mat){
    ret <- apply(.mat, 2, scale)
    ret[is.na(ret)] <- 0
    ret
}

ld.info <- fread(ld.file)
ld.info[, chr := as.integer(gsub("chr","",`chr`))]
ld.info[, query := paste0(`chr`, ":", `start`, "-", `stop`)]

.query <- ld.info[ld.index, ]$query

plink <- subset.plink(geno.hdr, ld.info[ld.index]$`chr`,
                      plink.lb = ld.info[ld.index]$`start`,
                      plink.ub = ld.info[ld.index]$`stop`,
                      temp.dir)

unlink(temp.dir, recursive=TRUE)

qtl.file <- qtl.dir %&% "/" %&% ld.index %&% ".txt.gz"
qtl.dt <- fread(qtl.file)

if(nrow(qtl.dt) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
    message("Empty QTL statistics after Q/C")
    q()
}

qtl.mu.dt <- qtl.dt %>%
    mutate(mu = `mean` * `alpha`) %>%
    mutate(beta.z = `beta` / `se`)

if(nrow(qtl.mu.dt) < 1){
    fwrite(data.table(), file=out.file, sep="\t")
    message("Empty QTL statistics after Q/C")
    q()
}

message("Read QTL results")

model.info <- qtl.mu.dt[, .(lfsr = min(lfsr)),
                       by = .(gene, celltype)]

qtl.mu <-
    dcast(qtl.mu.dt,
          physical.pos ~ celltype + gene,
          fill = 0,
          value.var = "mu",
          fun.aggregate = mean)

.models <- colnames(qtl.mu)[-1]

.matched <-
    plink$map %>%
    mutate(r = 1:n()) %>% 
    left_join(qtl.mu) %>%
    na.omit() %>%
    as.data.table()

.effect <- as.matrix(.matched[, ...models])

xx <- safe.scale(plink$bed[, .matched$r, drop = F])
y.hat <- safe.scale(xx %*% .effect)

colnames(y.hat) <- .models
rownames(y.hat) <- plink$fam$sample.ID

out.dt <- reshape2::melt(y.hat) %>%
    rename(projid = Var1, pgs = `value`) %>%
    as.data.table()

out.dt[, c("celltype","gene") := tstrsplit(`Var2`, split="_")]
out.dt[, Var2 := NULL]

out.dt <-
    out.dt %>%
    (function(x) cbind(`ld` = ld.index, x)) %>% 
    left_join(model.info)

fwrite(out.dt, out.file, sep = "\t")
message("done")
