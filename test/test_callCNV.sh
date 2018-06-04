#!/bin/bash

sample="GBM3"
treat="./data/bigwigs/GBM3_seq1_2.pfilter.tn5.srt.avg_cov.w1000K_s500K.bed"

ls -1 ./data/bigwigs/*w1000K_s500K.bed |\
  while read treat
  do 
  echo $treat;   
  sample=$(echo $treat| perl -ne '/(GBM[0-9]+[.0-9]*)/ && print "$1\n"' )
  Rscript callCNV.R -s "$sample" -t $treat 
  done
  
  

  
