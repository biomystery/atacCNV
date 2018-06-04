#!/bin/bash

# check commands: slopBed, bedGraphToBigWig and bedClip

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
    echo "Need 2 parameters! <tagAlign> <chrom size file>"
    exit
fi


# bedgraph to bigwig 
tag=$1 # input
chrsz=$2 # chrsz

#chrsz="/projects/ps-epigen/GENOME/${g}/${g}.chrom.sizes"
prefix=${tag%.tag*} # no >2 tag in the tagalign file name 
bdg=${prefix}.bdg
bdg_srt=${prefix}.srt.bedgraph
bigwig=${prefix}.bw
scale=

bedtools genomecov -i $tag -bg -g $chrsz | slopBed -i stdin -g $chrsz -b 0 | bedClip stdin $chrsz $bdg
sort -k1,1 -k2,2n $bdg > $bdg_srt
bedGraphToBigWig $bdg_srt $chrsz $bigwig
rm -f $bdg $bdg_srt

