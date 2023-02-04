argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

if(length(argv) != 2) q()

expr.file <- argv[1]    # e.g., expr.file = "result/step3/pb/Mic.rds"
out.file <- argv[2]

library(data.table)
library(dplyr)

.mkdir <- function(...) {
    dir.create(..., recursive = TRUE, showWarnings = FALSE)
}

`%&%` <- function(a,b) paste0(a,b)

.quantile.norm <- function(.mat) {
    stopifnot(is.matrix(.mat))
    ret <- .mat
    for(k in 1:ncol(ret)){
        x.k <- .mat[, k]
        .pos.k <- which(is.finite(x.k))
        x.k.valid <- x.k[.pos.k]
        ngenes <- length(x.k.valid)
        if(ngenes < 1) next
        qq <- qnorm((1:ngenes)/(ngenes + 1))
        x.k.valid[order(x.k.valid)] <- qq
        ret[.pos.k, k] <- x.k.valid
    }
    return(ret)
}

.sort.cols <- function(.mat, pheno) {
    .dt <- data.table(col = colnames(.mat))
    .dt[, c("projid", "ct") := tstrsplit(`col`, split="_")]
    .loc <- match(pheno$projid, as.integer(.dt$projid))
    .mat[, .loc, drop = FALSE]
}

## Remove genes with too few expression values
## - we would consider zero expression is missing
filter.mat <- function(.mat, .sum, missing.cutoff = 1, missing.rate = 0.8) {
    .row.missing.rate <- apply(t(.sum) < missing.cutoff, 2, mean)
    .qc <- .mat
    .qc[.sum < missing.cutoff] <- NA
    .qc[.row.missing.rate < missing.rate, , drop = FALSE]
}

.mkdir(dirname(out.file))

expr <- readRDS(expr.file)

message("Read expression data")

mu.mat <-
    filter.mat(expr$PB$mu, expr$PB$sum) %>% 
    .sort.cols(pheno = expr$pheno)

X <- apply(t(mu.mat), 2, scale)
X[is.na(X)] <- 0
.svd <- rsvd::rsvd(X, k = 100)

uu <- .svd$u
dd <- .svd$d
rownames(uu) <- colnames(mu.mat)
saveRDS(list(u = uu, d = dd), file = out.file)
