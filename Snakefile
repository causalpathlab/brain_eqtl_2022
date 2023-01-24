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

data_conditions = ["all","AD","noAD","APOE","noAPOE","male","female"]
data_ext = ["bed.gz","bed.gz.tbi"]
NPC = 50

rule step3:
    input:
        expand("result/step3/pb/{ct}.rds", ct=celltypes),
        expand("result/step3/qc/{ct}/{ct}_PC{pc}_{dt}.{ext}",
               dt=data_conditions, ext=data_ext, pc=NPC, ct=celltypes),
        expand("result/step3/qc/{ct}/{ct}_AD_all.{ext}",
               ext=data_ext, dt=["AD","PINE"], ct=celltypes)

rule step3_dropbox:
    shell:
        "rsync -argv result/step3/qc/* ~/Dropbox/AD430/1.Results/2.scRNA_pseudobulk/ --progress --exclude=\"*temp*\""

rule step3_pb_celltype:
    input:
        mtx = "result/step2/sorted/{ct}.mtx.gz",
        pheno = "data/metadata_PFC_all_individuals_092520.tsv.gz"
    output:
        "result/step3/pb/{ct}.rds"
    shell:
        "mkdir -p result/step3/pb/;"
        "Rscript --vanilla script/pseudobulk.R {input.mtx} {input.pheno} {output}"

rule step3_pb_celltype_qc:
    input:
        rds = "result/step3/pb/{ct}.rds",
        feat_file = "result/step1/features_annotated_GRCh37.txt.gz"
    output:
        ["result/step3/qc/{ct}/{ct}_PC%d_%s.%s"%(NPC, dt, ext)
         for dt in data_conditions
         for ext in data_ext] +
        ["result/step3/qc/{ct}/{ct}_%s_all.%s"%(dt, ext)
         for ext in data_ext
         for dt in ["AD", "PINE"]]
    shell:
        "Rscript --vanilla script/pseudobulk_qc.R {input.rds} {input.feat_file} %d result/step3/qc/{wildcards.ct}"%(NPC)

rule step3_rsync_up:
    shell:
        "rsync -argv ./result/step3 numbers:/home/ypark/work/brain_eqtl_2022/result/ --exclude=\"*temp\" --progress --size-only"

###############################
# genotype Q/C and queue jobs #
###############################

rule step4:
    input:
        expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])

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

TRAITS = [tr.split(".")[0] for tr in
          list_files("data/gwas/",
                     pattern=".bed.gz.tbi",
                     full_name = False)]

##############################
# When everything is done... #
##############################

rule step4_dropbox:
    shell:
        "mkdir -p ~/Dropbox/AD430/1.Results/4.GWAS;"
        "rsync -argv result/step4/gwas/*.vcf.gz* ~/Dropbox/AD430/1.Results/4.GWAS/;"
        "rsync -argv result/step4/combined ~/Dropbox/AD430/1.Results/3.eQTL --progress --exclude=\"*.vcf\""

rule step4_combine_qtl:
    input:
        expand("result/step4/combined/{adj}/{cond}/{ct}.vcf.gz",
               adj=("PC"+str(NPC)), ct=celltypes,
               cond=["all","AD","noAD","APOE","noAPOE","male","female"]),
        expand("result/step4/combined/{adj}/{cond}/{ct}.vcf.gz",
               adj=["AD","PINE"], ct=celltypes, cond="all")

rule _step4_combine_qtl_job:
    input:
        dirname="result/step4/qtl/{adj}/{cond}/{ct}",
        ldfile="data/LD.info.txt"
    output:
        vcf="result/step4/combined/{adj}/{cond}/{ct}.vcf.gz",
        tbi="result/step4/combined/{adj}/{cond}/{ct}.vcf.gz.tbi"
    shell:
        "mkdir -p result/step4/combined/{wildcards.adj}/{wildcards.cond}; "
        "Rscript --vanilla script/combine_qtl.R {input.dirname} {input.ldfile} {output.vcf}"

rule step4_run_qtl:
    input:
        expand("result/step4/qtl/{adj}/{cond}/{ct}/{ld}.txt.gz",
               ct=celltypes, adj=["AD","PINE"],
               cond="all", ld=range(1,1704))

rule _step4_run_qtl_job:
    input:
        expr="result/step3/qc/{ct}/{ct}_{adj}_{cond}.bed.gz",
        ld="data/LD.info.txt",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output: "result/step4/qtl/{adj}/{cond}/{ct}/{ld}.txt.gz"
    shell:
        "mkdir -p result/step4/qtl/{wildcards.adj}/{wildcards.cond}/{wildcards.ct}; "
        "Rscript --vanilla script/call_qtl.R {input.ld} {wildcards.ld} result/step4/rosmap {input.expr} {output}"

rule step4_combine_gwas:
    input:
        expand("result/step4/gwas/{trait}.vcf.gz",trait = TRAITS),
        expand("result/step4/gwas/{trait}.vcf.gz.tbi",trait = TRAITS)

rule _step4_combine_gwas_job:
    input:
        dirname="result/step4/gwas/{trait}",
        ldfile="data/LD.info.txt"
    output:
        vcf="result/step4/gwas/{trait}.vcf.gz",
        tbi="result/step4/gwas/{trait}.vcf.gz.tbi"
    shell:
        "mkdir -p result/step4/gwas/;"
        "Rscript --vanilla script/combine_gwas.R {input.dirname} {input.ldfile} {output.vcf}"

rule step4_run_gwas:
    input:
        expand("result/step4/gwas/{trait}/{ld}.txt.gz",
               trait = TRAITS,
               ld = range(1,1704))

rule _step4_run_gwas_job:
    input:
        ld="data/LD.info.txt",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output: "result/step4/gwas/{trait}/{ld}.txt.gz"
    shell: "Rscript --vanilla script/call_gwas_susie.R {input.ld} {wildcards.ld} result/step4/rosmap data/gwas/{wildcards.trait}.bed.gz {output}"

rule step4_run_twas:
    input:
        expand("result/step4/twas/{cond}/{ld}.txt.gz",
               cond = ["all","AD","noAD","male","female","APOE","noAPOE"],
               ld = range(1,1704))

rule _step4_run_twas_job:
    input:
        ld="data/LD.info.txt",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        "result/step4/twas/{cond}/{ld}.txt.gz"
    shell:
        "mkdir -p result/step4/twas/{wildcards.cond}; "
        "Rscript --vanilla script/call_twas_coloc.R {input.ld} {wildcards.ld} result/step4/rosmap result/step4/gwas result/step4/combined/PC50/{wildcards.cond} {output}"


##################################################
# When we need to run many, many jobs in cluster #
##################################################

rule step4_queue:
    input:
        expand("jobs/step4/{script}_{ct}_{adj}_{cond}.sh",
               script="qtl", ct=celltypes, adj=("PC"+str(NPC)),
               cond=["all","female","male","AD","noAD","APOE","noAPOE"]),
        expand("jobs/step4/{script}_{ct}_{adj}_{cond}.sh",
               script="qtl", ct=celltypes, adj=["AD","PINE"],
               cond="all"),
        expand("jobs/step4_gwas/{trait}.sh", trait=TRAITS)

rule _step4_queue_gwas_job_file:
    input:
        ld="data/LD.info.txt",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"]),
        exe="script/call_gwas_susie.R"
    output:
        script="jobs/step4_gwas/{trait}.sh"
    run:
        mkdir("jobs/step4_gwas/")
        with open(output.script,"w") as fh:
            sys.stdout=fh
            print("""#!/bin/bash -l
#SBATCH -J %(TRAIT)s
#SBATCH -o .log
#SBATCH -e .log
#SBATCH -D ./
#SBATCH -B 1
#SBATCH -t 1:00:00
#SBATCH --mem=4096
#SBATCH --array=1-1703

source /home/${USER}/.bashrc
source activate general

ld_file=%(LD_FILE)s
trait=%(TRAIT)s
script=%(EXE)s
gwas_file=data/gwas/%(TRAIT)s.bed.gz

ld_index=${SLURM_ARRAY_TASK_ID}

logdir=log/${script}_${trait}
mkdir -p ${logdir}/
outdir=result/step4/gwas/${trait}
[ -d $outdir ] || mkdir -p $outdir

outfile=${outdir}/${ld_index}.txt.gz
logfile=${logdir}/$(echo $outfile | awk '{ gsub("/","_"); print }')

[ -f $logfile ] && rm $logfile
if ! [ -f $outfile ]; then
    Rscript --vanilla ${script} ${ld_file} ${ld_index} result/step4/rosmap ${gwas_file} ${outfile}  >> $logfile 2>&1
fi
[ -f $logfile ] && rm $logfile

"""%{"LD_FILE": input.ld, "TRAIT": wildcards.trait, "EXE": input.exe})

rule _step4_queue_qtl_job_file:
    input:
        expr="result/step3/qc/{ct}/{ct}_{adj}_{cond}.bed.gz",
        ld="data/LD.info.txt",
        geno=expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"])
    output:
        script="jobs/step4/{script}_{ct}_{adj}_{cond}.sh"
    run:
        mkdir("jobs/step4/")
        with open(output.script,"w") as fh:
            sys.stdout=fh
            print("""#!/bin/bash -l
#SBATCH -J %(ADJ)s_%(COND)s_%(CELLTYPE)s
#SBATCH -o .log
#SBATCH -e .log
#SBATCH -D ./
#SBATCH -B 1
#SBATCH -t 3:00:00
#SBATCH --mem=16384
#SBATCH --array=1-1703

source /home/${USER}/.bashrc
source activate general

ld_file=%(LD_FILE)s
ct=%(CELLTYPE)s
adj=%(ADJ)s
cond=%(COND)s
script=%(EXE)s
expr_file=%(EXPR_FILE)s

ld_index=${SLURM_ARRAY_TASK_ID}

logdir=log/${script}_${tt}_${ct}
mkdir -p ${logdir}/
outdir=result/step4/${script}/${adj}/${cond}/${ct}
[ -d $outdir ] || mkdir -p $outdir

outfile=${outdir}/${ld_index}.txt.gz
logfile=${logdir}/$(echo $outfile | awk '{ gsub("/","_"); print }')

[ -f $logfile ] && rm $logfile
if ! [ -f $outfile ]; then
    Rscript --vanilla script/call_${script}.R ${ld_file} ${ld_index} result/step4/rosmap ${expr_file} ${outfile}  >> $logfile 2>&1
fi
[ -f $logfile ] && rm $logfile

"""%{"LD_FILE": input.ld, "EXE": wildcards.script, "CELLTYPE": wildcards.ct,
     "EXPR_FILE": input.expr, "ADJ": wildcards.adj, "COND": wildcards.cond})

rule step4_rsync_up:
    shell:
        "rsync -argv ./jobs numbers:/home/ypark/work/brain_eqtl_2022/ --exclude=\"*temp\" --progress;"
        "rsync -argv ./result/step4/rosmap* numbers:/home/ypark/work/brain_eqtl_2022/result/step4/ --exclude=\"*temp\" --progress --size-only"

rule step4_rsync_dn:
    shell:
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/gwas ./result/step4/ --exclude=\"*temp\" --progress;"
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step4/qtl ./result/step4/ --exclude=\"*temp\" --progress"




###################################
# trans-eQTL mediated by cis-eQTL #
###################################

rule step5:
    shell: "echo \"trans-eQTL mediated by cis-eQTL\""

rule step5_:
    input:
        rd = "data/markers.eTFset.RData",
        feat = "result/step1/features_annotated_GRCh37.txt.gz"
    output:
        full = "result/step5/tf_target_celltype.txt.gz",
        key = "result/step5/tf_celltype.txt.gz"
    shell:
        "mkdir -p result/step5/; "
        "Rscript --vanilla script/parse_tf_target_list.R {input.rd} {input.feat} {output.full} {output.key}"

rule step5_pb:
    input:
        expand("result/step5/expr/{ct}.{ext}", ct=celltypes, ext=["bed.gz","pca.rds"])

rule step5_pb_celltype:
    input:
        rds = "result/step3/pb/{ct}.rds",
        feat_file = "result/step1/features_annotated_GRCh37.txt.gz"
    output:
        bed = "result/step5/expr/{ct}.bed.gz",
        pca = "result/step5/expr/{ct}.pca.rds"
    shell:
        "mkdir -p result/step5;"
        "Rscript --vanilla script/pseudobulk_qc_nopca.R {input.rds} {input.feat_file} {output.bed};"
        "Rscript --vanilla script/pseudobulk_pca.R {input.rds} {output.pca};"

rule step5_queue:
    input:
        KEY= "result/step5/tf_celltype.txt.gz",
        TF_TARGET = "result/step5/tf_target_celltype.txt.gz",
        geno = expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"]),
        expr = expand("result/step5/expr/{ct}.bed.gz",ct=celltypes),
        pca = expand("result/step5/expr/{ct}.pca.rds",ct=celltypes)
    output:
        script="jobs/step5/run.sh"
    run:
        mkdir("jobs/step5")
        with open(output.script,"w") as fh:
            sys.stdout=fh
            print("""#!/bin/bash -l
#SBATCH -o .log
#SBATCH -e .log
#SBATCH -D ./
#SBATCH -B 1
#SBATCH -t 4:00:00
#SBATCH --mem=16384
#SBATCH --array=1-296

tf=${SLURM_ARRAY_TASK_ID}
cmd=$(gzip -cd %(KEY)s | head -n ${tf} | tail -n1 | awk -vTF_TARGET=%(TF_TARGET)s '{ tf=$1; ct=$2; data="result/step5/expr/"; output="result/step5/trans/"; print "Rscript --vanilla %(EXE)s %(TF_TARGET)s" FS tf FS "%(GENO)s" FS (data ct ".bed.gz") FS (data ct ".pca.rds") FS (output ct "/" tf); }')
exec $cmd
"""%{"KEY": input.KEY, "TF_TARGET": input.TF_TARGET, "GENO": "result/step4/rosmap", "EXE": "script/call_mediating_tf_target.R"})

rule step5_rsync_up:
    shell:
        "rsync -argv ./jobs numbers:/home/ypark/work/brain_eqtl_2022/ --exclude=\"*temp\" --progress;"
        "rsync -argv ./result/step5 numbers:/home/ypark/work/brain_eqtl_2022/result/ --exclude=\"*temp\" --progress --size-only"

rule step5_rsync_dn:
    shell:
        "rsync -argv numbers:/home/ypark/work/brain_eqtl_2022/result/step5/trans ./result/step5/ --exclude=\"*temp\" --progress --size-only"

rule step5_dropbox:
    shell:
        "rsync -argv result/step5/trans ~/Dropbox/AD430/1.Results/3.eQTL/ --progress"

rule step5_run_:
    input:
        linking = "result/step5/tf_target_celltype.txt.gz",
        geno = expand("result/step4/rosmap.{ext}", ext=["bed","bim","fam"]),
        expr = "result/step5/expr/{ct}.bed.gz",
        pca = "result/step5/expr/{ct}.pca.rds"
    output:
        data = "result/step5/trans/{ct}/{tf}.data.gz",
        stat = "result/step5/trans/{ct}/{tf}.stat.gz"
    shell:
        "mkdir -p result/step5/trans/{wildcards.ct};"
        "Rscript --vanilla script/call_mediating_tf_target.R {input.linking} {wildcards.tf} result/step4/rosmap {input.expr} {input.pca} result/step5/trans/{wildcards.ct}/{wildcards.tf}"
