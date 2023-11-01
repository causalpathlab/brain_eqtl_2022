



library(data.table)
library(dplyr)

herit.file <- "result/step4/combined/ld_heritability_PC37.txt.gz"
mqtl.file <- "result/step4/combined/mqtl_PC37.vcf.gz"

mqtl.dt <- fread(mqtl.file)
herit.dt <- fread(herit.file)

## select valid genes (CV R > .1)
## select valid variants and levels (p-value < 1e-4 or lfsr < .05)

valid.genes <-
    herit.dt[rr > .1, .(gene, celltype)] %>%
    unique()

valid.vars <-
    mqtl.dt[gene %in% valid.genes$gene &
            ((lfsr < .05 & abs(`z`) > 3) | `p.val` < 1e-4),
            .(`#chromosome`,
              `physical.pos`,
              `gene`,
              `levels`)] %>%
    unique()

mqtl.marg.feat <-
    valid.vars %>%
    left_join(mqtl.dt) %>%
    as.data.table()

mqtl.marg.feat[, z.marg := `beta`/`se`]

.temp <- mqtl.marg.feat %>% 
    dcast(`#chromosome` +
          `physical.pos` +
          `gene` +
          `levels` ~ `celltype`,
          value.var = "z.marg",
          fun.aggregate = mean,
          fill = 0)




xx <- as.matrix(.temp[gene == "BIN1"][, -(1:4)])
heatmap(xx)




TODO







qtl.dt[order(abs(z), - alpha, decreasing = T),
       head(.SD, 1),
       by = .(`#chromosome`,
              `physical.pos`,
              `celltype`)]


## identify clusters driven by genes



