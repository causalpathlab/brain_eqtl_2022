argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

mct.dir.name <- argv[1]
eqtl.dir.name <- argv[2]
ld.file <- argv[3]
out.file <- argv[4]

LFSR.CUTOFF <- 0.1
P.CUTOFF <- 5e-4

library(data.table)
library(dplyr)

ld.info <- fread(ld.file, header=TRUE)
ld.info[, ld.idx := 1:.N]

eqtl.files <- list.files(eqtl.dir.name, pattern = ".vcf.gz$", full.names=TRUE)
mct.files <- list.files(mct.dir.name, pattern = ".txt.gz$", full.names=TRUE)

## consolidate cell-type-by-cell-type eQTL analysis with multi-cell-type analysis
grab.mct.dt <- function(mct.file, eqtl.files, lfsr.cutoff = LFSR.CUTOFF, p.cutoff = P.CUTOFF){

    require(data.table)
    require(dplyr)
    data.table::setDTthreads(1)

    mct.dt <- fread(mct.file)
    if(nrow(mct.dt) < 1) return(data.table())

    .ld <- gsub(".txt.gz$","",basename(mct.file))

    ## 0. Read multi cell type results
    mct.dt[`var` <= 1e-8, var := 0]  # Deal with numerical zero
    mct.dt[`var` <= 1e-8, lfsr := 1] # This is not a reliable estimate
    mct.dt[`var` <= 1e-8, mean := 0] # Essentially empty
    mct.dt[`var` <= 1e-8, z := 0]    # Essentially empty

    mct.dt <- mct.dt %>%
        dplyr::mutate(se = sqrt(`var`)) %>%
        dplyr::rename(mct.level = `levels`,
                      mct.pip = `alpha`,
                      mct.mean = `mean`,
                      mct.se = `se`,
                      mct.z = `z`,
                      mct.lodds = `lodds`,
                      mct.lfsr = lfsr) %>%
        dplyr::mutate(mct.p.val = 2 * pnorm(abs(mct.z), lower.tail=FALSE)) %>%
        dplyr::select(-`lbf`, -`var`) %>%
        as.data.table

    chr <- unlist(mct.dt[1,1])
    lb <- min(mct.dt$physical.pos)
    ub <- max(mct.dt$physical.pos)
    query <- paste0(chr,":",lb,"-",ub)

    ## 1. Select top snp-gene pairs to investigate
    top.eqtl.dt <- data.table()

    for(eqtl.file in eqtl.files){

        .ct <- gsub(".vcf.gz$", "", basename(eqtl.file))

        .dt <-
            fread(cmd=paste0("tabix ", eqtl.file, " ", query, " -h"), header = TRUE) %>%
            dplyr::rename(physical.pos = snp.loc) %>%
            dplyr::select(physical.pos, hgnc_symbol, plink.a1, plink.a2, beta, se, p.val, pip, lfsr) %>%
            dplyr::mutate(celltype = .ct) %>%
            as.data.table()

        top.eqtl.dt <- rbind(top.eqtl.dt, .dt)
    }

    ## 2. Match with multi celltype SuSiE statistics
    .by <- c("physical.pos","hgnc_symbol","celltype")

    ret <- left_join(top.eqtl.dt, mct.dt, by=.by) %>% na.omit()

    ## 3. Select SNPs that are significant in at least one cell type
    .backbone <-
        ret[mct.lfsr < lfsr.cutoff & mct.p.val < p.cutoff & lfsr < lfsr.cutoff,
                     .(physical.pos, hgnc_symbol, mct.level)] %>%
        unique()

    ret <- left_join(.backbone, ret, by = c("physical.pos", "hgnc_symbol", "mct.level"))

    ret <- ret %>%
        dplyr::select(`#chromosome_name`,
                      `physical.pos`,
                      `hgnc_symbol`,
                      `plink.a1`,
                      `plink.a2`,
                      `celltype`,
                      `mct.level`,
                      `mct.pip`,
                      `mct.mean`,
                      `mct.se`,
                      `mct.z`,
                      `mct.lodds`,
                      `mct.lfsr`,
                      `beta`,
                      `se`,
                      `p.val`,
                      `pip`,
                      `lfsr`) %>%
        dplyr::mutate(LD = .ld) %>%
        as.data.table
}

out.raw.file <- gsub(".gz$", "", out.file)
unlink(out.raw.file)
num.threads <- 54

for(.chr in 1:22){
    .idx <- ld.info[as.character(`chr`) == paste0("chr", .chr)]$ld.idx
    .files <- paste0(mct.dir.name, "/", .idx, ".txt.gz")

    cl <- parallel::makeCluster(num.threads)
    parallel::clusterExport(cl, varlist=c())

    .list <- parallel::parLapply(cl, .files, grab.mct.dt,
                                 eqtl.files = eqtl.files,
                                 lfsr.cutoff = LFSR.CUTOFF,
                                 p.cutoff = P.CUTOFF)
    parallel::stopCluster(cl)
    .dt <- do.call(what = rbind, .list)
    rm(.list); gc()
    .dt.order <- .dt[, do.call(order, .SD), .SDcols=c("#chromosome_name", "physical.pos")]
    .dt <- as.data.table(.dt[.dt.order])

    fwrite(.dt,
           file = out.raw.file,
           sep = "\t",
           eol = "\n",
           append = (.chr!=1),
           row.names = FALSE,
           col.names = (.chr==1),
           na = "NA",
           quote = FALSE)

    rm(.dt); gc()
}

system(paste0("bgzip ", out.raw.file))
system(paste0("tabix -p vcf ", out.file))
unlink(out.raw.file)
