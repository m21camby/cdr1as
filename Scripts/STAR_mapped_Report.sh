#!/bin/bash

#################################
# Description: 
# 22 Nov. 2019 by SJK
# This file is to make easy report of mRNA unmapped from STAR mapping Log.final.out
# Two columns will be crated and easy to load on R as data frame
# How to use: ./STAR_mapped_Report.sh > mRNA_mapped_Report.txt 
#################################

# An array of experiment samples
#Array=("A1_1_S1" "A1_2_S2" "A2_1_S3" "A2_2_S4" "A3_1_S5" "A3_2_S6" "A4_1_S7" "A4_2_S8" "N1_1_S9" "N1_2_S10" "N2_1_S11" "N2_2_S12" "N3_1_S13" "N3_2_S14" "N4_1_S15" "N4_2_S16")
Array=$1


for i in ${Array[*]}; do

#echo $i

ARG1=`awk '$1 == "Number" && $3 == "input"' $i'Log.final.out' | awk '{print $1"_"$2"_"$3"_"$4, $6}'`
echo $i $ARG1
ARG2=`awk '$1 == "Uniquely" && $2 == "mapped" && $4 == "number"' $i'Log.final.out' | awk '{print $1"_"$2"_"$3"_"$4, $6}'`
echo $i $ARG2
ARG3=`awk '$1 == "Uniquely" && $2 == "mapped" && $4 == "%"' $i'Log.final.out' | awk '{sub(/\%/, "", $6); print $1"_"$2"_"$3"("$4")", $6}'`
echo $i $ARG3


done

