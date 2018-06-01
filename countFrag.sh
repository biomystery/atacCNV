#!/bin/bash
cd /home/zhc268/gbm/analysis 
samples=(`cat samples.txt`)

if [ -z $PBS_ARRAYID ]
then
    nth=$1
else
    nth=$PBS_ARRAYID
fi

s=${samples[$nth]}
#echo $s
#sort -k1,1 -k2,2n ${s}.fragments.bed > ${s}.sorted.fragments.bed
s1=${s}.sorted.peak.summit.bed
coverageBed -a ./$s1 -b  ${s}.sorted.fragments.bed -sorted  | awk '{print $7}' > ${s}.summit.frag.count.txt

