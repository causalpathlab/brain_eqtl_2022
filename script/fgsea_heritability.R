argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

if(length(argv) < 2) q()

in.file <- argv[1]
out.file <- argv[2]

library(data.table)
library(dplyr)

dir.create(dirname(out.file), recursive = T, showWarnings = F)

make.gs.lol <- function(.dt) {
    .dt <- as.data.table(.dt) %>% unique()
    .list <-
        .dt[, .(gene = .(gene_symbol)), by = .(gs_name)] %>%
        as.list()
    .names <- .list$gs_name
    .ret <- .list$gene
    names(.ret) <- .names
    return(.ret)
}

.dt <- fread(in.file)

.feat <-
    .dt[order(-rr),
        head(.SD, 1),
        by = .(gene, celltype)] %>%
    dcast(gene ~ celltype,
          value.var = "rr",
          fill = NA,
          fun.aggregate = max)

Z <- as.matrix(.feat[,-1])
colnames(Z) <- colnames(.feat)[-1]
rownames(Z) <- unlist(.feat[,1])

#########################
## Gene ontology words ##
#########################

gsea.dt <- data.table()

c5.gs.db <-
    msigdbr::msigdbr(species = "human", category = "C5") %>%
    dplyr::filter(`gs_subcat` %in% c("GO:BP","GO:CC","GO:MF")) %>%
    as.data.table()

gs.lol <- make.gs.lol(c5.gs.db[, .(gs_name, gene_symbol)])

for(j in 1:ncol(Z)){

    .ct <- colnames(Z)[j]
    .scores <- unlist(scale(Z[, j]))
    names(.scores) <- rownames(Z)
    .scores <- .scores[!is.na(.scores)]

    .gsea <- fgsea::fgsea(pathways = gs.lol,
                          stats = .scores,
                          minSize = 10,
                          maxSize = 1000,
                          nproc = 16,
                          nPermSimple = 5000,
                          scoreType = "pos")

    .gsea[, gs_name := pathway]
    .gsea[, pathway := NULL]
    .gsea[, celltype := .ct]

    gsea.dt <- rbind(gsea.dt, .gsea)
}

#########################
## Pathway annotations ##
#########################

.sub <- c("CP:KEGG","CP:BIOCARTA","CP:PID","CP:REACTOME","CP:WIKIPATHAYS")

c2.gs.db <-
    msigdbr::msigdbr(species = "human", category = "C2") %>%
    dplyr::filter(gs_subcat %in% .sub) %>%
    as.data.table()

gs.lol <- make.gs.lol(c2.gs.db[, .(gs_name, gene_symbol)])

for(j in 1:ncol(Z)){

    .ct <- colnames(Z)[j]
    .scores <- unlist(scale(Z[, j]))
    names(.scores) <- rownames(Z)
    .scores <- .scores[!is.na(.scores)]

    .gsea <- fgsea::fgsea(pathways = gs.lol,
                          stats = .scores,
                          minSize = 10,
                          maxSize = 1000,
                          nproc = 16,
                          nPermSimple = 5000,
                          scoreType = "pos")

    .gsea[, gs_name := pathway]
    .gsea[, pathway := NULL]
    .gsea[, celltype := .ct]

    gsea.dt <- rbind(gsea.dt, .gsea)
}

.key <- rbind(unique(c2.gs.db[, .(gs_cat, gs_subcat, gs_name)]),
              unique(c5.gs.db[, .(gs_cat, gs_subcat, gs_name)]))

out.dt <- gsea.dt %>% left_join(.key) %>% as.data.table()

out.dt[, padj := p.adjust(pval, method="fdr"), by = .(gs_subcat, celltype)]

fwrite(out.dt, out.file)
