all: AD.vcf.gz AD.vcf.gz.tbi hg38ToHg19.over.chain.gz

hg38ToHg19.over.chain.gz:
	wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz" -O $@

AD.vcf.gz: raw/GCST90027158_buildGRCh38.tsv.gz hg38ToHg19.over.chain.gz
	./liftover_effect.sh $^ 3 4 $@

%.vcf.gz.tbi: %.vcf.gz
	tabix -p vcf $<

