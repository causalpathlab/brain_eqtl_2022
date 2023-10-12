options(stringsAsFactors = FALSE)

`%&%` <- function(a,b) paste0(a, b)
`%r%` <- function(mat,rr) mat[rr, , drop = FALSE]
`%c%` <- function(mat,cc) mat[, cc, drop = FALSE]

library(tidyverse)
library(data.table)
library(ggrepel)

num.int <- function(...) format(..., justify="none", big.mark=",", drop0trailing = TRUE)

num.sci <- function(...) format(..., justify="none", digits=2, scientific = TRUE)

row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D = proxy::dist(mat, method <- function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] = 0
    h.out = hclust(D)
    o.out = cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

col.order <- function(pair.tab, .ro, ret.tab = FALSE) {

    M <- pair.tab %>%
        dplyr::select(row, col, weight) %>%
        dplyr::mutate(row = factor(row, .ro)) %>%
        tidyr::spread(key = col, value = weight, fill = 0)

    co <- order(apply(M[, -1], 2, which.max), decreasing = TRUE)
    .co <- colnames(M)[-1][co]
    if(ret.tab) {
        ret <- pair.tab %>%
            dplyr::mutate(row = factor(row, .ro)) %>%
            dplyr::mutate(col = factor(col, .co))
    } else {
        ret <- .co
    }
    return(ret)
}

order.pair <- function(pair.tab, ret.tab=FALSE) {

    require(tidyr)
    require(dplyr)

    .tab <- pair.tab %>% dplyr::select(row, col, weight)

    M <- .tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr <- M[, 1] %>% unlist(use.names = FALSE)
    cc <- colnames(M)[-1] %>% unlist(use.names = FALSE)

    ## Sort rows
    ro <- row.order(M %>% dplyr::select(-row) %>% as.matrix())

    ## Sort columns
    co <- order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    if(ret.tab){
        ret <- pair.tab %>%
            dplyr::mutate(row = factor(row, rr[ro])) %>%
            dplyr::mutate(col = factor(col, cc[co]))
    } else {
        ret <- list(rows = rr[ro], cols = cc[co], M = M)
    }

    return(ret)
}


.sort.matrix <- function(.X) {
    as.matrix(.X) %>%
        reshape2::melt() %>%
        rename(row = Var1, col = Var2, weight = value) %>%
        order.pair(ret.tab=TRUE) %>%
        as.data.table %>%
        dcast(row ~ col, value.var = "weight") %>%
        dplyr::select(-row) %>%
        as.matrix
}

.rnorm <- function(d1, d2) {
    matrix(rnorm(d1 * d2), d1, d2)
}

###############################################################
.matshow <- function(.mat, .lab = 1, .size = .1, .scale=TRUE) {

    .mat <- as.matrix(.mat)
    .cols <- colnames(.mat)
    if(length(.cols) < ncol(.mat)){
        colnames(.mat) <- str_c("c", 1:ncol(.mat))
    }
    .cols <- colnames(.mat)
    .rows <- str_c("r", 1:nrow(.mat))

    .dt <-
        as.data.table(.mat) %>%
        dplyr::mutate(row = str_c("r", 1:dplyr::n())) %>%
        as.data.table %>%
        melt(id.vars = "row", variable.name = "col") %>%
        dplyr::mutate(row = factor(as.character(row), rev(.rows))) %>%
        dplyr::mutate(col = factor(as.character(col), .cols))

    ret <-
        ggplot(.dt, aes(y = row, x = col, fill = value)) +
        theme_void() +
        theme(legend.position = "none") +
        geom_tile(size = .size, colour = "gray")

    if(.scale){
        ret <- ret +
            scale_fill_gradient2(low="blue", high="red", midpoint=0)
    } else {
        ret <- ret +
            scale_fill_distiller(palette = "Greys", direction = 1)
    }

    if(.lab > 0) {
        ret <- ret +
            geom_text(aes(label = round(value,1)), size = .lab)
    }

    return(ret)
}

################################################################
.gg.plot <- function(...) {
    ggplot(...) +
        theme_classic() +
        theme(axis.title = element_text(size=8)) +
        theme(axis.text = element_text(size=6)) +
        theme(legend.spacing = unit(.1, "lines"),
              legend.key.size = unit(.5, "lines"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=6),
              panel.background = element_rect(fill='transparent'),
              plot.background = element_rect(fill='transparent', color=NA),
              legend.background = element_rect(fill='transparent', size=0.05),
              legend.box.background = element_rect(fill='transparent'))
}

as.dt <- function(x, col.names=NULL) {
    .mat <- as.matrix(x)
    if(is.null(col.names)) col.names <- str_c("V",1:ncol(.mat))
    colnames(.mat) <- col.names
    as.data.table(.mat)
}

################################################################
if.needed <- function(.file, .code) {
    if(!all(file.exists(unlist(.file)))){ .code }
    stopifnot(all(file.exists(unlist(.file))))
}

################################################################
take.jaccard <- function(.mat){
    .mat[.mat > 0] <- 1
    num <- crossprod(.mat)
    tot <- apply(.mat, 2, sum)
    denom <- sweep(-num, 1, tot, `+`)
    denom <- sweep(denom, 2, tot, `+`)
    jacc <- num/denom
}

################################################################

build.annoy.graph <- function(Z, knn){
    stopifnot(length(rownames(Z)) == nrow(Z))

    annoy <- new(RcppAnnoy::AnnoyAngular, ncol(Z))
    for(i in 1:nrow(Z)){
        annoy$addItem(i-1, Z[i,])
    }
    annoy$build(50)

    xx <- c()
    yy <- c()
    dd <- c()
    for(i in 1:nrow(Z)){
        xx <- c(xx, rep(i, knn))
        neigh <- annoy$getNNsByItem(i-1, knn)
        neigh.dd <- sapply(neigh, function(j) annoy$getDistance(i-1, j))
        yy <- c(yy, neigh + 1)
        dd <- c(dd, neigh.dd)
    }

    nn <- nrow(Z)
    mm <- length(xx)
    A <- Matrix::spMatrix(nrow=nn, ncol=nn, i=xx, j=yy, x=exp(-dd))
    rownames(A) <- rownames(Z)
    colnames(A) <- rownames(Z)

    return(A)
}

louvain.umap <- function(A){

    A <- as.matrix(A)
    S <- A * t(A) 
    A <- (A + t(A))/2 * (S > 0)    # symetrize
    rm(S)

    G <- igraph::graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE)

    A <- Matrix::Matrix(A, sparse=T)

    res <- 2

    nrepeat <- 20
    C <- list(modularity = 0)
    for(r in 1:nrepeat){
        c.r <- igraph::cluster_louvain(G, resolution = res)
        message(paste("modularity: ", c.r$modularity, "\n"))
        if(r == 1 || max(c.r$modularity) > max(C$modularity)){
            C <- c.r
        }
    }

    umap.A <- uwot::optimize_graph_layout(A, verbose = TRUE, spread=10,
                                          n_sgd_threads = "auto")

    data.table(genes = igraph::V(G)$name, membership = C$membership) %>%
        left_join(cbind(data.table(genes = rownames(A)), as.data.table(umap.A)))

}
