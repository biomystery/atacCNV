#!/bin/bash

# check commands
which subtractBed &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }

# check args
if [ $# -lt 2 ];then
    echo "Need 2 parameters! <reads|fragments (bed)> <mask regions(bed)>"
    exit
fi

# receive args
reads=$1 # bed format
region=$2 # bed format 

# output names
out=${reads/.bed/.filt.bed}

subtractBed -a $reads -b $region -A  | sort -k1,1 -k2,2n  > $out

echo -e "input reads:\t$(cat $reads | wc -l)" > ${out/.bed/.log}
echo -e "output reads:\t$(cat $ | wc -l)">> ${out/.bed/.log}

# eg:

# 
