#!/bin/bash 

# check commands:

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bigWigAverageOverBed &>/dev/null || { echo "bigWigAverageOverBed not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }


# end of checking

if [ $# -lt 4 ];then
    echo "Need 4 parameters! <bw|bigwig> <window step size><window size><chromsize file>"
    exit
fi

# receive input 
bwIn=$1
step=$2
swindow=$3
chromsize=$4

# define output 
sample=${bwIn/.bw/};sample=${sample/.bigwig/};

outTab=${sample}.avg_cov.w$[$swindow/1000]K_s$[$step/1000]K.tab
outBed=${sample}.avg_cov.w$[$swindow/1000]K_s$[$step/1000]K.bed
outBw=${sample}.avg_cov.w$[$swindow/1000]K_s$[$step/1000]K.bw

## main
# make window 
bedtools makewindows -g $chromsize -w $swindow -s $step | \
    awk -v OFS='\t' '{print $0,NR}'> tmp.bed

echo -e  "total $(cat tmp.bed|wc -l) sliding windows with step=${step}bp size=${swindow}bp "

# average over window 
bigWigAverageOverBed $bwIn tmp.bed $outTab -bedOut=$outBed

awk  '{printf "%s\t%d\t%d\t%f\n", $1,($2+$3)/2,($2+$3)/2+1,$5}'  $outBed  > tmp.bdg
bedSort tmp.bdg tmp.bdg
bedGraphToBigWig  tmp.bdg  $chromsize $outBw

# delete tmp files 
rm tmp.bed tmp.bdg


# egs:
#smoothBw.sh ../data/hg38_75bp/GRCh38_75.bw  500000 1000000 ~/data/GENOME/hg38/hg38.chrom.sizes
