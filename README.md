# Single-cell eQTL calling on 29 Brain cell types

## Step 1. Merge data across 14 different batches

Before we merge data sets, we performed gentle Q/C steps to remove "empty" features (genes) in the following criteria:

- the standard deviation across cells $<$ 0.01

- the mean expression $<$ 1e-4

- the number of non-zero cells $<$ 10

We built a quite large data set consisting of gene expression vectors across 2.4M cells measured on 26k+ features

## Step 2. Sort the full data into smaller data sets for each cell type



## Step 3. Expression data preprocessing

### Streamlined data processing by PCA

1. Constructed pseudobulk expression matrix (gene x samples/individuals/`projid`)
2. Perform sample by sample quantile normalization to values generated by standard normal distribution, namely $N(0,1)$
3. For each chromosome, adjust non-genetic effects on the genes located in the chrososome by taking top 50 principal components estimated by genes in the other chromosomes (leaving one chrososome out at a time, LOCO)

`${celltype}/${celltype}_PC50_all.bed.gz`

#### a. AD vs. non-AD

We can partition the samples into two groups to carry out condition-specific eQTL analysis:

`${celltype}/${celltype}_PC50_AD.bed.gz`
`${celltype}/${celltype}_PC50_noAD.bed.gz`

#### b. APOE (e4) vs. APOE WT (without e4)

`${celltype}/${celltype}_PC50_APOE.bed.gz`
`${celltype}/${celltype}_PC50_noAPOE.bed.gz`

#### c. Female vs. male

`${celltype}/${celltype}_PC50_female.bed.gz`
`${celltype}/${celltype}_PC50_male.bed.gz`


### Experimental data processing

#### Confounder adjustment by cell-level kNN matching

We can match AD vs. non-AD cells and selectively remove putative confounding factors present in both types of cells derived from AD and non-AD samples:

`${celltype}/${celltype}_AD_all.bed.gz`

#### Confounder adjustment by individual-level pairwise matching

We can also remove individual-level confounding factors by matching cells between individuals:

`${celltype}/${celltype}_PINE_all.bed.gz`


#### Data sharing

Naming convention:

`${celltype}/${celltype}_${data_processing}_${condition}.bed.gz`

How to extract gene expression vectors:

```
$ tabix Mic/Mic_AD_all.bed.gz 19:45409000-45410000 -h | cut -f 1-10
```

We can retrieve genes by their chromosome name and approximate locations (e.g., LD block):

```
#chromosome_name	tss	tes	ensembl_gene_id	hgnc_symbol	11409232_Mic	11336574_Mic	10260309_Mic	10248033_Mic	20207013_Mic
19	45409011	45412650	ENSG00000130203	APOE	-0.00219808405473058	0.0633145290915309	0.28734561628366	0.123355747962476	0.264699779472434
```

Columns:

1. `#chromosome_name` : chromosome name (between 1 and 22)
2. `tss` : transcription start site (the left most location in `hg19`)
3. `tes` : transcription end site (the right most location in `hg19`)
4. `ensembl_gene_id` : ENSEMBL gene ID
5. `hgnc_symbol` : human-readable gene symbol
6. `${projid}_${celltype}` : many samples

## 2. Expression QTL calling within each cell type 


## 3. 

