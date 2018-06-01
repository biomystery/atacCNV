#!/bin/bash

# check commands: slopBed, bedGraphToBigWig and bedClip

which macs2 &>/dev/null || { echo "macs2 not found! Download macs2: <https://github.com/taoliu/MACS/>"; exit 1; }
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 3 ];then
    echo "Need  at least 3 parameters! <tagalign|bed> <chrom size file> <window size>"
    exit
fi


# bedgraph to bigwig 
input=$1 # input 
chrsz=$2 # chrsz
window=$3


prefix=${input%.*}
bdg=${prefix}_wind_${window}.bdg
bw=${prefix}_wind_${window}.bw

macs2 pileup -i $input -B --extsize $[window/2] -o ${bdg}
slopBed -i $bdg -g $chrsz -b 0| bedClip stdin $chrsz ${bdg}_clip # clip
sort -k1,1 -k2,2n ${bdg}_clip > ${bdg} # sort 
bedGraphToBigWig $bdg $chrsz $bw # track 
rm ${bdg}_clip
