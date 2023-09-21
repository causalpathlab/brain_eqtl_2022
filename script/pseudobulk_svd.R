argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()
sum.file <- argv[1] # "result/step3/sum.bed.gz"
out.file <- argv[2]

library(data.table)

X <- fread(sum.file)
X <- as.matrix(X[, -(1:6)])

.mean <- apply(log1p(X), 1, mean, na.rm=T)
ntop <- 10000
top.genes <- order(.mean, decreasing = TRUE)[1:ntop]

xx <- t(X[top.genes, , drop = F])
xx <- apply(xx, 2, scale)
xx[is.na(xx)] <- 0

.svd <- rsvd::rsvd(xx, k = 100)

rownames(.svd$u) <- colnames(X)

saveRDS(.svd, out.file)
