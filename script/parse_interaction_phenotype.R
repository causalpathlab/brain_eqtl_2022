
pheno.file <- "data/metadata_PFC_all_individuals_092520.tsv.gz"

library(dplyr)
library(data.table)

.vars <- c("projid",            # project ID, unique
           "apoe_genotype",     # factors --> convert them 
           "msex",              # 0/1: 0=female, 1=male
           "smoking",           # 0/1: smoking
           "hypertension_bl",   # 0/1: hypertension
           "cancer_bl",         # 0/1: cancer history
           "headinjrloc_bl",    # 0/1: head injury with loss of consciousness
           "diabetes_sr_rx_bl", # 0/1: diabetes
           "chf_bl",            # 0/1: congestive heart failure
           "claudication_bl",   # 0/1: history of claudication is a marker of peripheral vascular disease
           "heart_bl",          # 0/1: history of heart conditions
           "stroke_bl",         # 0/1: history of stroke
           "thyroid_bl",        # 0/1: history of thyroid disease
           "cogdx",             # quant, 0 - 6
           "educ",              # quant, years
           "alcohol_g_bl",      # quant, grams of alcohol per day at baseline
           "phys5itemsum_bl",   # quant, Physical activity (5 items)
           "gpath",             # quant, global AD pathology burden
           "niareagansc",       # quant, The modified NIA-Reagan diagnosis of Alzheimer's disease
           "amyloid",           # quant, amyloid, immunohistochem + image
           "plaq_d",            # quant, plaque burden across regions
           "plaq_n",            # quant, neuritic plaque burden across regions
           "tangles",           # quant, tau tangles
           "nft",               # quant, neurofibrillary tangles
           "tdp_stage4",        # quant, 1: amygdala; 2: amygdala+limbic; 3: amygdala+limbic+neocortical
           "caa_4gp"            # quant, Cerebral amyloid angiopathy
           )

.dt <-
    fread(pheno.file) %>% 
    select(all_of(.vars)) %>%
    mutate(apoe_genotype = if_else(stringr::str_detect(apoe_genotype, "4"), 1, 0)) %>% 
    select(-ends_with(".NA")) %>% 
    as.data.table()

fwrite(.dt, "data/metadata_selected.tsv.gz", col.names=T, sep="\t")

