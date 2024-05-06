#!/bin/bash


# ********************* #
# find top 10 reads from unmapped reads
# ********************* #

./findTopTen.sh WTC1_S1 &

# 1. For unmapped reads

ARG=$1

echo Top ten unmapped 1st read by $i >> $ARG'_TopTenUnmappedReads.txt'
awk '(NR%4==2)' $ARG'Unmapped_R1.fastq' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' |sort -k 2,2 -nr |head -10 >> $ARG'_TopTenUnmappedReads.txt'



# 2. For rRNA mapped reads

for i in {241..248}
	do
		echo Top ten suspected mycoplasma read by $i >> TopTenUnmappedReads.txt 
		awk '{A[$1]++}END{for(i in A)print i,A[i]}' NR_$i'_rRNA_mapped.fasta' |sort -k 2,2 -r|head -10 >> TopTenUnmappedReads.txt 
	
	done



