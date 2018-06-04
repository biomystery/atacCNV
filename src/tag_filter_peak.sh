#!/bin/bash
cd /home/zhc268/gbm/analysis 
libs=($(cat ../../including_libs.txt))

if [ -z $PBS_ARRAYID ]
then
    i=$1
else
    i=$PBS_ARRAYID
fi

s=${libs[2*i+1]};    l=${libs[2*i]}
s1=/home/zhc268/data/outputs/tagAligns/${l}_R1.trim.PE2SE.nodup.tn5.tagAlign.gz
subtractBed -a $s1 -b naive.merged.sorted.bed -A  | sort -k1,1 -k2,2n  >  ${s}.pfilter.tn5.srt.tagAlign


#qsub ./tag_filter_peak.sh -t 0-9 -q condo -l nodes=1:ppn=1 -N tag_f
