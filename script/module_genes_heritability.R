argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

## herit.file <- "result/step4/combined/ld_heritability_PC37.txt.gz"
## expr.file <- "result/step3/sum.bed.gz"
## out.file <- "temp.txt.gz"

herit.file <- argv[1]
expr.file <- argv[2]
out.file <- argv[3]

dir.create(dirname(out.file), recursive = T, showWarnings = F)

library(data.table)
library(dplyr)
library(mmutilR)

R2.CUTOFF <- .1

qc.z.mat <- function(.feat){
    Z <- as.matrix(.feat[, -1])
    colnames(Z) <- colnames(.feat[, -1])
    rownames(Z) <- unlist(.feat[,1])
    cutoff <- quantile(Z, .99, na.rm = T)
    Z[is.na(Z)] <- 0
    Z[Z > cutoff] <- cutoff
    return(Z)
}

write.temp.data <- function(z_gene_celltype, hdr){
    unlink(unlist(fileset.list(hdr)))
    Zt <- Matrix::Matrix(t(z_gene_celltype), sparse=T)
    write.sparse(Zt, rownames(Zt), colnames(Zt), hdr)
}

degree.cutoff <- function(G, .cutoff = 3) {
    G.sub <- G
    n.remove <- sum(igraph::degree(G.sub) < .cutoff)
    while(n.remove > 0){
        vv <- igraph::V(G.sub)
        .retain <- vv[igraph::degree(G.sub) >= .cutoff]
        G.sub <- igraph::induced_subgraph(G.sub, .retain)
        n.remove <- sum(igraph::degree(G.sub) < .cutoff)
    }
    return(G.sub)
}

comp.cutoff <- function(G, .cutoff = 10){
    .comp <- igraph::components(G)
    .kk <- which(.comp$csize >= .cutoff) # valid component(s)
    .valid <- which(.comp$membership %in% .kk) # valid nodes
    G <- igraph::induced_subgraph(G, igraph::V(G)[.valid])
    return(G)
}

run.knn.umap <- function(.knn, .names,
                         res = 1,
                         nrepeat = 50,
                         min.deg = 0,
                         min.size = 0, ...){

    .dt <-
        data.table(from = c(.knn$src.index, .knn$tgt.index),
                   to = c(.knn$tgt.index, .knn$src.index),
                   weight = c(.knn$weight, .knn$weight))

    message("Sorting out pairs")

    .df <- .dt[`from` != `to`,
               .(weight=mean(weight)),
               by = .(from, to)] %>%
        as.data.frame()

    message("Building a kNN graph object")

    G <- igraph::graph_from_data_frame(.df, directed=FALSE)

    message("Resolving graph clusters")

    C <- list(modularity = 0)
    for(r in 1:nrepeat){
        c.r <- igraph::cluster_louvain(G, resolution = res)
        message(paste("modularity: ", c.r$modularity, "\n"))
        if(r == 1 || max(c.r$modularity) > max(C$modularity)){
            C <- c.r
        }
    }

    .clust <- data.table(gene = .names[as.integer(C$names)],
                         membership = C$membership)

    ## this seems very important for good visualization
    if(min.deg > 0 && min.size > 0){
        G.qc <- degree.cutoff(G, min.deg)
        G.qc <- comp.cutoff(G.qc, min.size)
        A <- igraph::as_adj(G.qc, attr="weight")
    } else {
        A <- igraph::as_adj(G, attr="weight")
    }

    umap.A <- uwot::optimize_graph_layout(A, verbose = TRUE,
                                          n_sgd_threads = "auto",
                                          ...)

    data.table(gene = .names[as.integer(rownames(A))],
               umap1 = umap.A[,1],
               umap2 = umap.A[,2]) %>%
        (function(x) left_join(.clust, x, by = "gene"))
}

################################################################
## temporary data files

temp.dir <- paste0(out.file, "_temp")
dir.create(temp.dir, recursive = T, showWarnings = F)

null.data.hdr <- paste0(temp.dir, "/null")
obs.data.hdr <- paste0(temp.dir, "/obs")

################################################################
## Take cell-type-specific average expressions
expr.dt <- fread(expr.file)

.keys <- c("#chromosome_name",
           "tss",
           "tes",
           "ensembl_gene_id",
           "hgnc_symbol",
           "celltype")

key.cols <- which(colnames(expr.dt) %in% .keys)
.expr.mat <- as.matrix((expr.dt[, .SD, .SDcols = - key.cols]))

.means <- apply(log1p(.expr.mat), 1, mean, na.rm=T)

Z0.dt <- expr.dt[, .(hgnc_symbol, celltype)]
Z0.dt[, v := .means]

.feat0 <- dcast(Z0.dt, hgnc_symbol ~ celltype,
                value.var = "v",
                fill = 0,
                fun.aggregate = max)

Z0 <- qc.z.mat(.feat0)
null.data <- write.temp.data(Z0, null.data.hdr)

#####################################################
## Step 1. Gene modules driven by gene expressions ##
#####################################################
nn <- rcpp_mmutil_info(null.data$mtx)$max.col
mm <- rcpp_mmutil_info(null.data$mtx)$max.row

.bbknn.null <-
    rcpp_mmutil_bbknn_pca(mtx_file = null.data$mtx,
                          r_batches = rep("no.batch", nn),
                          knn = 50,
                          RANK = 25,
                          NUM_THREADS = 16,
                          RECIPROCAL_MATCH = TRUE,
                          TAKE_LN = FALSE,  # already did log1p
                          COL_NORM = 10000, # these are actual genes
                          EM_ITER = 10)

null.umap <- run.knn.umap(.bbknn.null$knn,
                          readLines(null.data$col),
                          nrepeat = 20)

##########################################################
## Step 2. Gene modules driven by heritability patterns ##
##########################################################

.dt <- fread(herit.file)

.feat <-
    .dt[rr > R2.CUTOFF & `gene` %in% null.umap$gene,
        .(gene, celltype, rr)] %>% 
    dcast(gene ~ celltype,
          value.var = "rr",
          fill = 0,
          fun.aggregate = max)

Z <- round(qc.z.mat(.feat) * 100)

obs.data <- write.temp.data(Z, obs.data.hdr)

.batch <-
    fread(obs.data$col, col.names = "gene", header = F) %>%
    left_join(null.umap)

.bbknn <- rcpp_mmutil_bbknn_pca(mtx_file = obs.data$mtx,
                                r_batches = .batch$membership,
                                knn = 100,
                                RANK = 25,
                                NUM_THREADS = 16,
                                TAKE_LN = FALSE,
                                COL_NORM = 10000,
                                EM_ITER = 10)

obs.umap <- run.knn.umap(.bbknn$knn.adj,
                         readLines(obs.data$col),
                         res = 3,
                         nrepeat = 20,
                         min.deg = 3,
                         min.size = 10)

.obs.dt <- fread(herit.file) %>%
    dcast(gene ~ celltype,
          value.var = "rr",
          fill = 0,
          fun.aggregate = max)

out.dt <-
    obs.umap %>%
    left_join(null.umap, by = "gene", suffix = c(".genetic",".expr")) %>%
    left_join(.obs.dt, by = "gene") %>% 
    as.data.table()

fwrite(out.dt, out.file, sep="\t", row.names = F, col.names = T)
unlink(temp.dir, recursive = T)
