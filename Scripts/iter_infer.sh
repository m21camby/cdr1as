!#/bin/bash

# infer reads conditions


EXP=("WTC1_S1" "WTC2_S3" "WTC3_S5" "WTC4_S7" "WTC5_S9" "WTm7oex1_S2" "WTm7oex2_S4" "WTm7oex3_S6" "WTm7oex4_S8" "WTm7oex5_S10" "KOC1_S11" "KOC2_S13" "KOC3_S15" "KOC4_S17" "KOC5_S19" "KOm7oex1_S12" "KOm7oex2_S14" "KOm7oex3_S16" "KOm7oex4_S18" "KOm7oex5_S20")



for i in ${EXP[*]}
do

echo $i

infer_experiment.py -r mm10_RefSeq.bed -i $i'sorted_Aligned.out.bam' 

done

