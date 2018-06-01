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
# -v bed

coverageBed -a $bed -b  ${s}.pfilter.tn5.srt.tagAlign  -sorted | awk -v OFS='\t' '{print $4,$7}' > ${s}_${bed}.pfilter.coverage.txt

#qsub ./countFrag_genome.sh -t 0-9 -v bed=hg38.1M.windows_a.sorted.bed -q condo -l nodes=1:ppn=1 -N 1M_a




