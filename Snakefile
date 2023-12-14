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
        "nice Rscript script/merge_data.R result/step1/qc result/step1/merged"

rule step1_row_annotate_GRCh37:
    input: "result/step1/valid_features.txt.gz"
    output: "result/step1/features_annotated_GRCh37.txt.gz"
    shell:
        "mkdir -p result/step1/;"
        "nice Rscript --vanilla script/annotate_genes_GRCh37.R {input} {output}"

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
        "nice Rscript script/qc_row_mtx.R {input.valid} {params.data} {params.out}"

rule step1_glob_row:
    input: expand("result/step1/score/{sample}.rows.gz", sample=SAMPLES)
    output: "result/step1/valid_features.txt.gz"
    shell: "mkdir -p result/step1; nice Rscript script/test_qc_rows.R result/step1/score {output}"

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
        "nice Rscript --vanilla script/qc_score.R {params.prefix} {output.row} {output.col}"

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
        "nice Rscript --vanilla script/simplify_celltype_file.R {input} {output}"


rule step2_sort_celltype:
    input:
        mtx = "result/step1/merged.mtx.gz",
        annot = "result/step2/celltypes.txt.gz"
    output:
        mtx = "result/step2/sorted/{ct}.mtx.gz"
    shell:
        "mkdir -p result/step2/sorted;"
        "nice Rscript --vanilla script/select_celltype_mtx.R {input.mtx} {input.annot} {wildcards.ct} {output.mtx}"

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
        "nice Rscript --vanilla script/pseudobulk.R {input.mtx} {output}"

rule step3_pb_concat:
    input:
        expand("result/step3/pb/{ct}.rds", ct=celltypes),
        feat = "result/step1/features_annotated_GRCh37.txt.gz"

    output:
        sum = "result/step3/sum.bed.gz",
        mean = "result/step3/log_mean.bed.gz"

    shell:
        "mkdir -p result/step3/;"
        "nice Rscript --vanilla script/pseudobulk_concatenate.R {input.feat} {output.sum} {output.mean}"

rule step3_svd:
    input:
        "result/step3/sum.bed.gz"
    output:
        "result/step3/svd.rds"
    shell:
        "mkdir -p result/step3/;"
        "nice Rscript --vanilla script/pseudobulk_svd.R {input} {output}"

rule step3_svd_assoc:
    input:
        svd = "result/step3/svd.rds",
        pheno = "data/metadata_PFC_all_individuals_092520.tsv.gz"
    output:
        "result/step3/svd_pheno_assoc.txt.gz"
    shell:
        "mkdir -p result/step3/;"
        "nice Rscript --vanilla script/assoc_svd_pheno.R {input.svd} {input.pheno} {output}"

rule rsync_step3_up:
    shell:
        "rsync -argv ./result/step3/*.* numbers:/home/ypark/work/brain_eqtl_2022/result/step3/ --exclude=\"*temp\" --progress --size-only"

#######################
# Call eQTLs and TWAS #
#######################

## list(range(10,50,5)) + list(range(50,101,20)) + [100]
COVAR_PCs = list(range(10,50)) + list(range(50,101,10)) + [100]

rule step4:
    input:
        expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"]),
        expand("result/step4/combined/ld_heritability_PC{pc}.txt.gz", pc=COVAR_PCs),
        expand("result/step4/combined/qtl_PC{pc}.vcf.gz{ext}", pc=[37, 70, 100], ext=["",".tbi"]),
        expand("result/step4/combined/iqtl_PC{pc}.vcf.gz{ext}", pc=[37, 70, 100], ext=["",".tbi"]),
        expand("result/step4/combined/mqtl_PC{pc}.vcf.gz{ext}", pc=[37, 70, 100], ext=["",".tbi"]),
        expand("result/step4/combined/ld_twas_{gwas}_PC{pc}.txt.gz", pc=[37, 70, 100], gwas=["AD"]),
        expand("result/step4/combined/ld_itwas_{gwas}_PC{pc}.txt.gz", pc=[37, 70, 100], gwas=["AD"]),
        expand("result/step4/combined/coloc_{gwas}_PC{pc}.vcf.gz{ext}", pc=[37, 70, 100], gwas=["AD"], ext=["",".tbi"])

rule step4_dropbox:
    shell:
        "rsync -argv result/step4/combined/*heritability* ~/Dropbox/Writing/AD430/1.Results/3.eQTL/heritability/ --progress; "
        "rsync -argv result/step4/combined/*qtl* ~/Dropbox/Writing/AD430/1.Results/3.eQTL/qtl/ --progress --size-only --exclude=\"*.vcf\"; "
        "rsync -argv result/step4/combined/*twas* ~/Dropbox/Writing/AD430/1.Results/3.eQTL/twas/ --progress --size-only; "
        "rsync -argv result/step4/combined/*coloc* ~/Dropbox/Writing/AD430/1.Results/3.eQTL/coloc/ --progress --size-only; "
        "echo \"Done\""

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

rule step4_post_vcf_jobs:
    input:
        "result/step4/combined/{jobname}.vcf"
    output:
        "result/step4/combined/{jobname}.vcf.gz"
    shell:
        "bgzip {input}; "

rule step4_post_vcf_jobs_index:
    input:
        "result/step4/combined/{jobname}.vcf.gz"
    output:
        "result/step4/combined/{jobname}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}; "

rule step4_post_vcf_sort:
    output:
        vcf = "result/step4/combined/{jobname}.vcf"
    params:
        ddir = lambda w: "/".join(w.jobname.split("_")),
        taboo = "#chromosome"
    shell:
        "mkdir -p result/step4/combined/; "
        "[ -d {output.vcf}_temp ] || mkdir -p {output.vcf}_temp; "
        "find result/step4/{params.ddir}/ -name *.gz -size +0 | xargs -I file zcat file | awk -F '\\t' -v TABOO={params.taboo} -v TEMP={output.vcf}_temp/ 'NR == 1 {{ print $0; next }} NR > 1 && $1 != TABOO {{ print $0 | (\"sort  -k1,1 -k2,2n -T \" TEMP) }}' > {output.vcf}; "
        "[ -d {output.vcf}_temp ] && rm -rf {output.vcf}_temp;"

rule step4_post_ld_jobs:
    output: "result/step4/combined/{jobname}.txt.gz"
    params:
        ddir = lambda w: "/".join(w.jobname.split("_")[1:]),
        taboo = lambda w: w.jobname.split("_")[0]
    shell:
        "mkdir -p result/step4/combined/; "
        "cat result/step4/{params.ddir}/*.txt.gz | gzip -cd | awk 'NR == 1 || $1 != \"{params.taboo}\"' | gzip -c > {output}"

rule step4_jobs:
    input:
        expand("jobs/step4/heritability_{nPC}.sh",  nPC=COVAR_PCs),
        expand("jobs/step4/mqtl_{nPC}.sh",  nPC=[37, 70, 100]),
        expand("jobs/step4/qtl_{nPC}.sh",  nPC=[37, 70, 100]),
        expand("jobs/step4/iqtl_{nPC}.sh",  nPC=[37, 70, 100]),
        expand("jobs/step4/coloc_{gwas}_{nPC}.sh",  gwas="AD", nPC=[37, 70, 100]),
        expand("jobs/step4/twas_{gwas}_{nPC}.sh",  gwas="AD", nPC=[37, 70, 100]),
        expand("jobs/step4/itwas_{gwas}_{nPC}.sh",  gwas="AD", nPC=[37, 70, 100])

rule rsync_step4_up:
    shell:
        "rsync -argv ./result/step4/rosmap* numbers:/home/ypark/work/brain_eqtl_2022/result/step4/ --exclude=\"*temp\" --progress --size-only"

rule rsync_step4_down:
    shell:
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/mqtl ./result/step4/ --exclude=\"*temp\" --progress --size-only; "
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/coloc ./result/step4/ --exclude=\"*temp\" --progress --size-only; "
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/twas ./result/step4/ --exclude=\"*temp\" --progress --size-only; "
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/itwas ./result/step4/ --exclude=\"*temp\" --progress --size-only; "
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/qtl ./result/step4/ --exclude=\"*temp\" --progress --size-only; "
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/iqtl ./result/step4/ --exclude=\"*temp\" --progress --size-only;"
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/heritability ./result/step4/ --exclude=\"*temp\" --progress --size-only;"

rule rsync_jobs_up:
    shell:
        "rsync -argv ./jobs numbers:/home/ypark/work/brain_eqtl_2022/ --exclude=\"*temp\" --progress;"

rule step4_run:
    input:
        expand("result/step4/{qtl}/PC{nPC}/{ld}.txt.gz",
               ld=range(1,1704),
               qtl=["qtl","iqtl","mqtl"],
               nPC=[37, 70, 100]),
        expand("result/step4/{twas}/{gwas}/PC{nPC}/{ld}.txt.gz",
               ld=range(1,1704),
               twas=["twas","itwas","coloc"],
               gwas="AD",
               nPC=[37, 70, 100]),
        expand("result/step4/epgs/PC{nPC}/{ld}.txt.gz",
               ld=range(1,1704),
               nPC=[37, 70, 100])

rule _step4_run_epgs:
    input:
        ldfile="data/LD.info.txt",
        qtl_dir="result/step4/qtl/",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/epgs/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/predict_epgs_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.qtl_dir}/PC{wildcards.nPC} {output}"

rule _step4_run_twas:
    input:
        ldfile="data/LD.info.txt",
        qtl_dir="result/step4/qtl/",
        gwas_stat_dir="data/gwas",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/twas/{gwas}/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/call_twas_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.qtl_dir}/PC{wildcards.nPC} {input.gwas_stat_dir}/{wildcards.gwas}.vcf.gz {output}"

rule _step4_run_itwas:
    input:
        ldfile="data/LD.info.txt",
        qtl_dir="result/step4/iqtl/",
        gwas_stat_dir="data/gwas",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/itwas/{gwas}/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/call_twas_conditional_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.qtl_dir}/PC{wildcards.nPC} {input.gwas_stat_dir}/{wildcards.gwas}.vcf.gz {output}"

rule _step4_run_coloc:
    input:
        ldfile="data/LD.info.txt",
        gwas_stat_dir="data/gwas",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/coloc/{gwas}/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/call_coloc_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.expr} {input.svd} {wildcards.nPC} {input.gwas_stat_dir}/{wildcards.gwas}.vcf.gz {output}"

rule _step4_jobs_twas:
    input:
        ldfile="data/LD.info.txt",
        qtl_dir="result/step4/qtl/",
        gwas_stat_dir="data/gwas",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/twas_{gwas}_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("twas",
                       "script/call_twas_per_ld.R",
                       "result/step4/twas/" + wildcards.gwas + "/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap", input.qtl_dir + "/PC" + wildcards.nPC, input.gwas_stat_dir + "/" + wildcards.gwas + ".vcf.gz"],
                       mem=2048,
                       maxtime="4:00:00")

rule _step4_jobs_coloc:
    input:
        ldfile="data/LD.info.txt",
        gwas_stat_dir="data/gwas",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/coloc_{gwas}_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("coloc",
                       "script/call_coloc_per_ld.R",
                       "result/step4/coloc/" + wildcards.gwas + "/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap",
                        input.expr, input.svd, wildcards.nPC,
                        input.gwas_stat_dir + "/" + wildcards.gwas + ".vcf.gz"],
                       mem=2048,
                       maxtime="4:00:00")

rule _step4_jobs_itwas:
    input:
        ldfile="data/LD.info.txt",
        qtl_dir="result/step4/iqtl/",
        gwas_stat_dir="data/gwas",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/itwas_{gwas}_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("twas",
                       "script/call_twas_conditional_per_ld.R",
                       "result/step4/itwas/" + wildcards.gwas + "/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap", input.qtl_dir + "/PC" + wildcards.nPC, input.gwas_stat_dir + "/" + wildcards.gwas + ".vcf.gz"],
                       mem=2048,
                       maxtime="4:00:00")

rule _step4_jobs_iqtl:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        herit = "result/step4/heritability",
        pheno = "data/metadata_PFC_all_individuals_092520.tsv.gz",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/iqtl_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("iqtl",
                       "script/call_iqtl_per_ld.R",
                       "result/step4/iqtl/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap", input.expr, input.herit, input.svd, input.pheno, wildcards.nPC],
                       mem = 2048,
                       maxtime = "4:00:00")

rule _step4_run_iqtl:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        herit = "result/step4/heritability",
        pheno = "data/metadata_PFC_all_individuals_092520.tsv.gz",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/iqtl/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/call_iqtl_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.expr} {input.herit} {input.svd} {input.pheno} {wildcards.nPC} {output}"

rule _step4_jobs_qtl:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/qtl_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("qtl",
                       "script/call_qtl_per_ld.R",
                       "result/step4/qtl/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap", input.expr, input.svd, wildcards.nPC],
                       mem=2048,
                       maxtime="4:00:00")

rule _step4_jobs_mqtl:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        queue="jobs/step4/mqtl_{nPC}.sh"
    run:
        mkdir("jobs/step4")
        with open(output.queue, "w") as fh:
            sys.stdout = fh
            print_Rjob("mqtl",
                       "script/call_multi_qtl_per_ld.R",
                       "result/step4/mqtl/PC" + wildcards.nPC,
                       [input.ldfile, "result/step4/rosmap", input.expr, input.svd, wildcards.nPC],
                       mem=2048,
                       maxtime="4:00:00")

rule _step4_run_qtl:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/qtl/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/call_qtl_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.expr} {input.svd} {wildcards.nPC} {output}"

rule _step4_run_mqtl:
    input:
        ldfile="data/LD.info.txt",
        expr="result/step3/log_mean.bed.gz",
        svd="result/step3/svd.rds",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/mqtl/PC{nPC}/{ld}.txt.gz"
    shell:
        "nice Rscript --vanilla script/call_multi_qtl_per_ld.R {wildcards.ld} {input.ldfile} result/step4/rosmap {input.expr} {input.svd} {wildcards.nPC} {output}"

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
                       [input.ldfile, "result/step4/rosmap", input.expr, input.svd, wildcards.nPC],
                       mem=2048,
                       maxtime="10:00:00")

##############################################
# downstream enrichment, clustering analysis #
##############################################

rule step5:
    input:
        expand("result/step5/module/mqtl_zscore_PC{nPC}.txt.gz", nPC = [37, 70, 100])

rule step5_dropbox:
    shell:
        "rsync -argv result/step5/module/* ~/Dropbox/Writing/AD430/1.Results/4.eQTL_modules/ --progress --exclude=\"*temp*\"; "

rule _step5_module:
    threads: 16
    input: "result/step4/combined/mqtl_PC{nPC}.vcf.gz"
    output: "result/step5/module/mqtl_zscore_PC{nPC}.txt.gz"        
    shell: "nice Rscript --vanilla script/module_genes_mqtl_zscore.R {input} {output}"

rule _step5_fgsea:
    input: "result/step4/combined/ld_mqtl_zscore_PC{nPC}.txt.gz"
    output: "result/step5/fgsea/mqtl_zscore_PC{nPC}.txt.gz"
    shell: "nice Rscript --vanilla script/fgsea_mqtl_zscore.R {input} {output}"
