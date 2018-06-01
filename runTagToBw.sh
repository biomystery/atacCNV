#! /bin/bash

cd /home/zhc268/gbm/analysis
tags=($(ls -1 *pfilter*srt*tagAlign))
s=${tags[$PBS_ARRAYID]}
g="/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes"
tagTobw.sh $s $g


#qsub -t 0-9 ./runTagToBw.sh -q condo -l nodes=1:ppn=1 
