#! /bin/bash

cd /home/zhc268/gbm/analysis
tags=($(ls -1 *pfilter*srt*tagAlign))
s=${tags[$PBS_ARRAYID]}
g="/projects/ps-epigen/GENOME/hg38/hg38.chrom.sizes"
bedToBdgBw.sh $s $g $window


#qsub -t 0-9 ./runTagToBw_window.sh -q condo -l nodes=1:ppn=1 -v window=500000
