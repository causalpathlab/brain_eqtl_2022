options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

## ld.index <- 207
## ld.file <- "data/LD.info.txt"
## geno.hdr <- "result/step4/rosmap"
## qtl.dir <- "result/step4/iqtl/PC37/"
## gwas.stat.file <- "data/gwas/AD.vcf.gz"
## out.file <- "output.txt.gz"

if(length(argv) < 6) q()

ld.index <- as.integer(argv[1])
ld.file <- argv[2]
geno.hdr <- argv[3]
qtl.dir <- argv[4]
gwas.stat.file <- argv[5]
out.file <- argv[6]

################################################################

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

compute.twas <- function(.svd, .mu, zz, pve.cutoff = .9){
    pve <- cumsum(.svd$d^2) / sum(.svd$d^2)
    k <- min(min(which(pve > pve.cutoff)), length(.svd$d))

    lambda <- 1/nrow(.svd$u)
    vd <- sweep(.svd$v[, 1:k, drop = F], 2, .svd$d[1:k], `*`)

    mu.vd <- t(.mu) %*% vd

    .num <- as.numeric(t(.mu) %*% zz)
    .denom <- apply(mu.vd, 1, function(x) sqrt(sum(x^2) + lambda))
    data.table(col = colnames(.mu),
               stwas.num = .num,
               stwas.denom = .denom,
               stwas.z = .num/.denom)
}

read.gwas <- function(gwas.file, .query){

    .ret.1 <- fread(cmd="tabix " %&% gwas.file %&% " chr" %&% .query %&% " -h", sep = "\t")
    .ret <- fread(cmd="tabix " %&% gwas.file %&% " " %&% .query %&% " -h", sep = "\t")

    print(gwas.file)

    .ret <- rbind(.ret, .ret.1) %>%
        as.data.table()

    .trait <- gsub(".vcf.gz$", "", basename(gwas.file))
    .ret[, trait := .trait]
    return(.ret[, .(position, variant_id, p_value, effect_allele, other_allele,
                    beta, standard_error)])
}

match.gwas.plink <- function(gwas.dt, plink){
    ## match with plink
    matched <-
        left_join(plink$map, gwas.dt) %>%
        na.omit() %>%
        as.data.table()

    if(nrow(matched) < 1) {
        return(data.table())
    }

    matched[effect_allele == allele1, beta.flip := beta]
    matched[effect_allele == allele2, beta.flip := -beta]
    matched[is.na(beta.flip), beta.flip := 0]
    matched[, z := `beta.flip`/`standard_error`]

    return(matched)
}

#####################################
## conventional summary-based TWAS ##
#####################################

create.twas.dt <- function(qtl.mu, gwas.dt, plink){

    .gwas <- match.gwas.plink(gwas.dt, plink)
    zz <- .gwas[order(abs(`z`), decreasing = T),
                head(.SD, 1),
                by = .(physical.pos)]
    zz <- as.data.table(qtl.mu)[, .(physical.pos)] %>%
        left_join(zz)
    zz <- as.matrix(zz[, .(z)])
    zz[is.na(zz)] <- 0

    .match <- match(qtl.mu$physical.pos, plink$map$physical.pos)
    xx <- safe.scale(plink$bed)[, .match, drop = F]
    .mu <- as.matrix(qtl.mu[, -1, drop = F])
    .svd <- rsvd::rsvd(xx/sqrt(nrow(xx)))

    compute.twas(.svd, .mu, zz, pve.cutoff = .95)
}

################################################################

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
                       by = .(gene, cond, W, celltype)]

qtl.mu <-
    dcast(qtl.mu.dt,
          physical.pos ~ cond + W + celltype + gene,
          fill = 0,
          value.var = "mu",
          fun.aggregate = mean)

gwas.dt <-
    read.gwas(gwas.stat.file, .query) %>%
    mutate(physical.pos = position) %>%
    as.data.table()

message("Read GWAS summary statistics")

stwas.stat <- create.twas.dt(qtl.mu, gwas.dt, plink)

stwas.stat[, c("cond", "W", "celltype", "gene") := tstrsplit(`col`, split="_")]
stwas.stat[, W := as.integer(W)]
stwas.stat[, `col` := NULL]
stwas.stat[, stwas.p.val := 2*pnorm(-abs(`stwas.z`), lower.tail=T)]

out.dt <- as.data.frame(stwas.stat) %>%
    dplyr::mutate(ld = ld.index) %>%
    dplyr::select(`ld`, `gene`, `cond`, `W`, `celltype`,
                  dplyr::starts_with("stwas")) %>%
    dplyr::left_join(model.info) %>%
    as.data.table()

fwrite(out.dt, out.file, sep="\t", col.names=T, row.names=F)

message("Done")
