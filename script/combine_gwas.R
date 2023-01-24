argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

dir.name <- argv[1]
ld.file <- argv[2]
out.file <- argv[3]

library(data.table)
library(dplyr)

ld.info <- fread(ld.file)
ld.info[, ld.idx := 1:.N]

read.dt.sorted <- function(dir.name, chr, num.threads=54){
    require(data.table)
    .idx <- which(ld.info$chr == chr | ld.info$chr == paste0("chr",chr))
    .files <- paste0(dir.name, "/", .idx, ".txt.gz")
    .files <- .files[file.exists(.files)]
    cl <- parallel::makeCluster(num.threads)
    parallel::clusterExport(cl, varlist=c())

    .read <- function(x) {
        options(stringsAsFactors = FALSE)
        data.table::setDTthreads(1)
        .sz <- file.info(x)$size
        if(.sz < 1024) return(data.table::data.table())
        ret <- data.table::fread(x,nThread=1,sep="\t")
        ret[, pip := signif(as.numeric(pip),2)]              # remove values too small
        ret[, theta := round(as.numeric(theta),4)]           #
        ret[, theta.sd := round(as.numeric(theta.sd),4)]     #
        ret[, coverage := signif(as.numeric(coverage),2)]    #
        ret[, lfsr := signif(as.numeric(lfsr),3)]            #
        return(ret)
    }

    .dt.list <- parallel::parLapply(cl, .files, .read)
    parallel::stopCluster(cl)
    data.table::setDTthreads(num.threads)
    .dt <- do.call(what=rbind, .dt.list)
    rm(.dt.list)
    gc()
    if(nrow(.dt) < 1) return(data.table())
    .dt.order <- .dt[, do.call(order, .SD), .SDcols=c("#chr", "snp.loc")]
    as.data.table(.dt[.dt.order])
}

out.raw.file <- gsub(".gz$", "", out.file)
unlink(out.raw.file)

for(chr in 1:22){

    .dt <- read.dt.sorted(dir.name, chr)

    fwrite(.dt,
           file = out.raw.file,
           sep = "\t",
           eol = "\n",
           append = (chr!=1),
           row.names = FALSE,
           col.names = (chr==1),
           na = "NA",
           quote = FALSE)
    rm(.dt)

    gc()
    print(chr)
}

## Rsamtools::bgzip(out.raw.file, dest=out.file)
## Rsamtools::indexTabix(out.file, format="vcf")
system(paste0("bgzip ", out.raw.file))
system(paste0("tabix -p vcf ", out.file))
unlink(out.raw.file)
