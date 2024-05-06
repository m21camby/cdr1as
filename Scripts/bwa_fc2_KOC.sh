#!/bin/bash

# SP081_045_S9_L002_R2_001.fastq.gz
# SP081_046_S10_L002_R2_001.fastq.gz
/data/rajewsky/sequencing/mouse/191031_NB501326_0339_AHCM7WBGXC/fastq/{sample}_R1_001.fastq.gz


for i in KOC1_S11 KOC2_S13 KOC3_S15 KOC4_S17 KOC5_S19 
do
echo $i

	bwa mem -t 40 -k 14 -T 1 -L 3,3 -O 6,6 -E 3,3 /data/rajewsky/indices/mm10_bwa_0.7.17/mm10.fa '/data/rajewsky/sequencing/mouse/191031_NB501326_0339_AHCM7WBGXC/fastq/'$i'_R1_001.fastq.gz' | samtools view -hbuS - > './'$i'_anchors.bam'
 
	samtools view -h $i'_anchors.bam' | python2 /data/rajewsky/home/skim/find_circ2-master/find_circ.py --genome /data/rajewsky/indices/mm10_bwa_0.7.17/mm10.fa --name $i --output $i 2> 'STDERR.fc2_NR_'$i'.txt' 
	echo "Finish $i bwa and fc2"

done






