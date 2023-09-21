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

#################################
# Step 2. Simply sort the cells #
#################################

celltypes = ["Ast", "Oli", "Opc", "Exc-NRGN", "Exc-L3-4-RORB-CUX2", "Exc-RELN-CHD7",
             "Exc-L2-3-CBLN2-LINC02306", "Exc-L5-6-RORB-LINC02196", "Exc-L6-THEMIS-NFIA",
             "Exc-L6b", "Exc-L4-5-RORB-GABRG1", "Exc-L4-5-RORB-IL1RAPL2", "Exc-L6-CT",
             "Exc-L3-5-RORB-PLCH1", "Exc-L5-6-IT-Car3", "Exc-L5-ET", "Exc-L5-6-NP",
             "Inh-PVALB", "Inh-SST", "Inh-VIP", "Inh-LAMP5", "Inh-PAX6", "Mic", "CAMs", "Tcell",
             "Fib", "Endo", "Per", "SMC"]

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

#####################################################
# Step 3. Combine cells and create pseudo-bulk data #
#####################################################

rule step3:
    input:
        expand("result/step3/pb/{ct}.rds", ct=celltypes),
        sum = "result/step3/sum.bed.gz",
        mean = "result/step3/log_mean.bed.gz",
        svd = "result/step3/svd.rds",
        assoc = "result/step3/svd_pheno_assoc.txt.gz"

rule step3_pb_celltype:
    input:
        mtx = "result/step2/sorted/{ct}.mtx.gz",
    output:
        "result/step3/pb/{ct}.rds"
    shell:
        "mkdir -p result/step3/pb/;"
        "Rscript --vanilla script/pseudobulk.R {input.mtx} {output}"

rule step3_pb_concat:
    input:
        expand("result/step3/pb/{ct}.rds", ct=celltypes),
        feat = "result/step1/features_annotated_GRCh37.txt.gz"

    output:
        sum = "result/step3/sum.bed.gz",
        mean = "result/step3/log_mean.bed.gz"

    shell:
        "mkdir -p result/step3/;"
        "Rscript --vanilla script/pseudobulk_concatenate.R {input.feat} {output.sum} {output.mean}"

rule step3_svd:
    input:
        "result/step3/sum.bed.gz"
    output:
        "result/step3/svd.rds"
    shell:
        "mkdir -p result/step3/;"
        "Rscript --vanilla script/pseudobulk_svd.R {input} {output}"

rule step3_svd_assoc:
    input:
        svd = "result/step3/svd.rds",
        pheno = "data/metadata_PFC_all_individuals_092520.tsv.gz"
    output:
        "result/step3/svd_pheno_assoc.txt.gz"
    shell:
        "mkdir -p result/step3/;"
        "Rscript --vanilla script/assoc_svd_pheno.R {input.svd} {input.pheno} {output}"

rule rsync_step3_up:
    shell:
        "rsync -argv ./result/step3/*.* numbers:/home/ypark/work/brain_eqtl_2022/result/step3/ --exclude=\"*temp\" --progress --size-only"

#######################
# Call eQTLs and TWAS #
#######################

rule step4:
    input:
        expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])

###############################
# genotype Q/C and queue jobs #
###############################

rule step4_prepare_genetic_data:
    input:
        pgen = "data/rosmap_wgs_2022/AD_WGS_20210117.pgen",
        pvar = "data/rosmap_wgs_2022/AD_WGS_20210117.pvar",
        psam = "data/rosmap_wgs_2022/AD_WGS_20210117.psam"
    output:
        bed = "result/step4/rosmap.bed",
        bim = "result/step4/rosmap.bim",
        fam = "result/step4/rosmap.fam"
    shell:
        "mkdir -p result/step4/;"
        "plink2 --pgen {input.pgen} "
        "--king-cutoff 0.15 "
        "--pvar {input.pvar} "
        "--psam {input.psam} "
        "--make-bed --out result/step4/rosmap"

rule step4_jobs_heritability:
    input:
        expand("jobs/step4/heritability_{nPC}.sh",
               nPC=list(range(10,101,20)) + [100])

rule rsync_step4_up:
    shell:
        "rsync -argv ./result/step4/rosmap* numbers:/home/ypark/work/brain_eqtl_2022/result/step4/ --exclude=\"*temp\" --progress --size-only"

rule rsync_jobs_up:
    shell:
        "rsync -argv ./jobs numbers:/home/ypark/work/brain_eqtl_2022/ --exclude=\"*temp\" --progress;"

rule _step4_jobs_heritability:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/heritability_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("heritability",
                       "script/call_heritability_per_ld.R",
                       "result/step4/heritability/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap", input.expr, input.svd],
                       mem=2048,
                       maxtime="2:00:00")




# rule step4_run_heritability:
#     input:
#         expand("result/step4/heritability/PC{nPC}/{ld}.txt.gz",
#                nPC=list(range(10,101,20)) + [100],
#                ld=range(1,1704))

# rule _step4_run_heritability_job:
#     input:
#         ldfile="data/LD.info.txt",
#         geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
#     output: "result/step4/heritability/PC{nPC}/{ld}.txt.gz"
#     shell:
#         "mkdir -p result/step4/heritability/PC{wildcards.nPC}/; "
#         "nice Rscript --vanilla script/call_heritability_per_ld.R {input.ldfile} {wildcards.ld} result/step4/rosmap result/step3/log_mean.bed.gz result/step3/svd.rds {wildcards.nPC} {output}"

# rule step4_run_qtl:
#     input:
#         expand("result/step4/qtl/{ld}.txt.gz",
#                ld=range(1,1704))

# rule _step4_run_qtl_job:
#     input:
#         ld="data/LD.info.txt",
#         geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
#     output: "result/step4/qtl/{ld}.txt.gz"
#     shell:
#         "mkdir -p result/step4/qtl/; "
#         "Rscript --vanilla script/call_qtl_per_ld.R {input.ld} {wildcards.ld} result/step4/rosmap result/step3/log_mean.bed.gz result/step3/svd.rds {output}"

# rule _step4_run_iqtl_job:
#     input:
#         ld="data/LD.info.txt",
#         pheno="data/metadata_PFC_all_individuals_092520.tsv.gz",
#         geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
#     output: "result/step4/iqtl/{ld}.txt.gz"
#     shell:
#         "mkdir -p result/step4/iqtl/; "
#         "Rscript --vanilla script/call_iqtl_per_ld.R {input.ld} {wildcards.ld} result/step4/rosmap result/step3/log_mean.bed.gz result/step3/svd.rds {input.pheno} {output}"
