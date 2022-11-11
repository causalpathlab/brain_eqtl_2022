from Util import *
import os
import gzip

SAMPLES = list_files("data/raw/", pattern="SM_", full_name = False)
EXT = ["rows.gz","cols.gz","mtx.gz","mtx.gz.index"]
SCR = ["rows","cols"]

##########################################################################
# Step 1. Data integration                                               #
# Step 1a. Select features consistently expressed across all the batches #
# Step 1b. Combine all the data                                          #
##########################################################################

rule step1:
    input:
        "result/step1/valid_features.txt.gz",
        "result/step1/features_annotated_GRCh37.txt.gz",
        expand("result/step1/merged.{ext}", ext=EXT)

rule step1_merge:
    input:
        expand("result/step1/qc/{sample}.{ext}", sample=SAMPLES, ext=EXT),
        "result/step1/valid_features.txt.gz"
    output:
        expand("result/step1/merged.{ext}", ext=EXT)
    shell:
        "Rscript script/merge_data.R result/step1/qc result/step1/merged"

rule step1_row_annotate_GRCh37:
    input: "result/step1/valid_features.txt.gz"
    output: "result/step1/features_annotated_GRCh37.txt.gz"
    shell:
        "mkdir -p result/step1/;"
        "Rscript --vanilla script/annotate_genes_GRCh37.R {input} {output}"

rule step1_row_qc:
    input:
        valid = "result/step1/valid_features.txt.gz",
    output:
        mtx = "result/step1/qc/{sample}.mtx.gz",
        idx = "result/step1/qc/{sample}.mtx.gz.index",
        row = "result/step1/qc/{sample}.rows.gz",
        col = "result/step1/qc/{sample}.cols.gz"
    params:
        data = "result/step1/temp/{sample}",
        out = "result/step1/qc/{sample}"
    shell:
        "Rscript script/qc_row_mtx.R {input.valid} {params.data} {params.out}"

rule step1_glob_row:
    input: expand("result/step1/score/{sample}.rows.gz", sample=SAMPLES)
    output: "result/step1/valid_features.txt.gz"
    shell: "mkdir -p result/step1; Rscript script/test_qc_rows.R result/step1/score {output}"

rule step1_score:
    input:
        mtx = "result/step1/temp/{sample}.mtx.gz",
        col = "result/step1/temp/{sample}.cols.gz",
        row = "result/step1/temp/{sample}.rows.gz"
    params: prefix = "result/step1/temp/{sample}"
    output:
        row = "result/step1/score/{sample}.rows.gz",
        col = "result/step1/score/{sample}.cols.gz"
    shell:
        "mkdir -p result/step1/score;"
        "Rscript --vanilla script/qc_score.R {params.prefix} {output.row} {output.col}"

#############################################################
# consistently map by including ENSEMBL ID and HGNC symbols #
#############################################################

rule step1_row_file:
    input: row = "data/raw/{sample}/filtered_feature_bc_matrix/features.tsv.gz"
    params: dir_ = "result/step1/temp"
    output: row = "result/step1/temp/{sample}.rows.gz"
    run:
        mkdir(params.dir_)
        with gzip.open(input.row) as f:
            rows = map(lambda xx: (xx[0] + '_' + xx[1]).encode(),
                       map(lambda x: x.decode().strip().split(), f))
            with gzip.open(output.row, 'wb') as f_out:
                f_out.write(b'\n'.join(rows) + b'\n')

rule step1_col_file:
    input:
        col = "data/raw/{sample}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        proj = "data/raw/{sample}/filtered_feature_bc_matrix/projid.txt"
    params: dir_ = "result/step1/temp"
    output: col = "result/step1/temp/{sample}.cols.gz"
    run:
        mkdir(params.dir_)
        with open(input.proj) as fh:
            projid = dict(map(lambda x: x.strip().split('\t'), fh))
        with gzip.open(input.col) as f:
            cols = map(lambda xx: (xx[0] + '_' + projid[xx[1]]).encode(),
                       map(lambda x: x.decode().strip().split('-'), f))
            with gzip.open(output.col, 'wb') as f_out:
                f_out.write(b'\n'.join(cols) + b'\n')

rule step1_mtx_file:
    input:
        mtx = "data/raw/{sample}/filtered_feature_bc_matrix/matrix.mtx.gz",
        row = "result/step1/temp/{sample}.rows.gz",
        col = "result/step1/temp/{sample}.cols.gz"
    params: dir_ = "result/step1/temp"
    output: mtx = "result/step1/temp/{sample}.mtx.gz"
    shell: "mkdir -p {params.dir_}; cp {input.mtx} {output.mtx}"

###########################################
# Step 2. Sort cells by known annotations #
###########################################

celltypes = ["Ast", "CAMs", "Endo",
             "Exc-L2-3-CBLN2-LINC02306", "Exc-L3-4-RORB-CUX2", "Exc-L3-5-RORB-PLCH1",
             "Exc-L4-5-RORB-GABRG1", "Exc-L4-5-RORB-IL1RAPL2",
             "Exc-L5-6-RORB-LINC02196", "Exc-L5-ET", "Exc-L5-6-IT-Car3", "Exc-L5-6-NP",
             "Exc-L6-CT", "Exc-L6-THEMIS-NFIA", "Exc-L6b", "Exc-NRGN", "Exc-RELN-CHD7", "Fib",
             "Inh-ALCAM-TRPM3", "Inh-CUX2-MSR1", "Inh-ENOX2-SPHKAP",
             "Inh-FBN2-EPB41L4A", "Inh-GPC5-RIT2",
             "Inh-L1-2-PAX6-SCGN", "Inh-L1-6-LAMP5-CA13", "Inh-L1-PAX6-CA4", "Inh-L3-5-SST-MAFB",
             "Inh-L5-6-PVALB-STON2", "Inh-L5-6-SST-TH", "Inh-L6-SST-NPY", "Inh-LAMP5-NRG1-Rosehip", "Inh-LAMP5-RELN",
             "Inh-PTPRK-FAM19A1", "Inh-PVALB-CA8-Chandelier", "Inh-PVALB-HTR4", "Inh-PVALB-SULF1", "Inh-RYR3-TSHZ2",
             "Inh-SGCD-PDE3A", "Inh-SORCS1-TTN", "Inh-VIP-ABI3BP", "Inh-VIP-CLSTN2", "Inh-VIP-THSD7B", "Inh-VIP-TSHZ2",
             "Mic", "Oli", "Opc", "Per", "SMC", "Tcell"]

rule step2:
    input: expand("result/step2/sorted/{ct}.mtx.gz", ct=celltypes)

rule step2_simplify:
    input: "data/consensus_annot_snPFC.joint_annotation_metadata.eQTL_celltype.tsv.gz"
    output: "result/step2/celltypes.txt.gz"
    shell:
        "mkdir -p result/step2;"
        "Rscript --vanilla script/simplify_celltype_file.R {input} {output}"

rule step2_sort_celltype:
    input:
        mtx = "result/step1/merged.mtx.gz",
        annot = "result/step2/celltypes.txt.gz"
    output:
        mtx = "result/step2/sorted/{ct}.mtx.gz"
    shell:
        "mkdir -p result/step2/sorted;"
        "Rscript --vanilla script/select_celltype_mtx.R {input.mtx} {input.annot} {wildcards.ct} {output.mtx}"
