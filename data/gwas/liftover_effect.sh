#!/bin/bash
if [ $# -lt 5 ]; then
    exit 1
fi

infile=$1
chain=$2
CHR_COL=$3
POS_COL=$4
outfile=$5

# infile=data/gwas/GCST90027158_buildGRCh38.tsv.gz
# chain=result/step4/hg38ToHg19.over.chain.gz
# CHR_COL=3
# POS_COL=4
# outfile=temp.vcf.gz

mkdir -p $(dirname $outfile)
coord=$(dirname $outfile)/$(basename $outfile .gz)_coord.gz
lifted=$(dirname $outfile)/$(basename $outfile .gz)_lifted

tempdir=${outfile}_temp
mkdir -p $tempdir

cat $infile | gzip -cd | tail -n+2 | \
    awk -vCHR=$CHR_COL -vPOS=$POS_COL -F'\t' \
'{
    chr=$(CHR); gsub(/chr/,"",chr);
    pos=int($(POS));
    print "chr" chr FS (pos-1) FS pos FS (chr ":" (pos-1) "-" pos)
}' | sort -T ${tempdir} -k1,1 -k2,2n | gzip > $coord

./bin/liftOver ${coord} ${chain} ${lifted} ${lifted}_remove

echo "Created Liftover"

cat $infile | gzip -cd | head -n1 | \
    awk -vCHR=$CHR_COL -vPOS=$POS_COL -F'\t' \
'{
    printf "#chromosome" FS "position";
    for(j=1; j<=NF; ++j){
        if(j != CHR && j != POS){
            printf FS $j;
        }
    }
    printf "\n";
}' | bgzip -c > $outfile

cat $infile | gzip -cd | tail -n+2 | \
    awk -vLIFTED=${lifted} -vCHR=$CHR_COL -vPOS=$POS_COL -F'\t' \
'BEGIN {
    while((getline line < LIFTED) > 0) {
        split(line,arr,"\t");
        valid[arr[4]] = arr[1] FS arr[3];
    }
}
{
    chr=$(CHR); gsub(/chr/,"",chr);
    pos=int($(POS));
    k =  (chr ":" (pos-1) "-" pos)
    if(k in valid){
       printf "%s", valid[k];
        for(j=1; j<=NF; ++j){
            if(j != CHR && j != POS){
                printf FS $j;
            }
        }
        printf "\n";
    }
}' | sort -T ${tempdir} -k1,1 -k2,2n | bgzip -c >> $outfile

echo "Finished sorting"

[ -f $lifted ] && rm ${lifted} ${lifted}_remove
[ -f $coord ] && rm $coord
[ -d $tempdir ] && rm -r $tempdir

echo "Done"

