import os

def mkdir(dir_):
    os.makedirs(dir_, exist_ok=True)

def list_files(dir_, pattern = None, full_name = True):
    ret = []
    if pattern != None: 
        ret = [f for f in os.listdir(dir_) if f.find(pattern) != -1]
    else:
        ret = [f for f in os.listdir(dir_)]

    if full_name:
        ret = [dir_ + "/" + f for f in ret]

    return ret

def print_Rjob(jobname: str, scriptName: str, outDir: str, argv: list,
               mem = 1024, maxtime = "1:00:00",
               minId=1, maxId=1703,
               file_ext="txt.gz"):
    print("""#!/bin/bash -l
#SBATCH -o .log
#SBATCH -e .log
#SBATCH -D ./
#SBATCH -B 1
#SBATCH -t %(TIME)s
#SBATCH --mem=%(MEM)d
#SBATCH --array=%(minID)d-%(maxID)d

#SBATCH -J %(jobname)s

logdir=log/%(EXE)s
mkdir -p ${logdir}/

jobid=${SLURM_ARRAY_TASK_ID}

outfile=%(OUTDIR)s/${jobid}.%(EXT)s
mkdir -p $(dirname $outfile)
logfile=${logdir}/$(echo $outfile | awk '{ gsub("/","_"); print }')
[ -f $logfile ] && rm $logfile

if ! [ -f $outfile ]; then
    Rscript --vanilla %(EXE)s ${jobid} %(ARGV)s ${outfile} >> $logfile 2>&1
fi
[ -f $logfile ] && rm $logfile
"""%{"jobname": jobname,
     "MEM": mem,
     "TIME": maxtime,
     "minID": minId,
     "maxID": maxId,
     "EXE": scriptName,
     "OUTDIR": outDir,
     "EXT": file_ext,
     "ARGV": " ".join(map(str, argv))})
