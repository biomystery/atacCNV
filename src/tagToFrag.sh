#!/bin/bash

# check commands
which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# check args
if [ $# -lt 2 ];then
    echo "Need 2 parameters! <tagAlign|tagAlign.gz> <chrom size file>"
    exit
fi


# bedgraph to bigwig 
tag=$1 # input
chrsz=$2 # chrsz

#chrsz="/projects/ps-epigen/GENOME/${g}/${g}.chrom.sizes"
prefix=${tag%.tag*} # no >2 tag in the tagalign file name
frag=${perfix}.bed 
bdg=${prefix}.bdg
bdg_srt=${prefix}.srt.bdg
bigwig=${prefix}.bw

# tagAlign to frag

if [ ${tag/*./} == 'gz' ]; then
    zcat $tag | awk 'NR%2{printf "%s ",$0;next;}1' | awk '{print $1,$2,$3,$8,$9}' | awk -v OFS='\t' '{min=$2;max=$2;for(i=2;i<=5;i++) {if ($i>max) max=$i;if($i<min) min=$i;} print $1,min,max}' > $frag
else
    cat $tag | awk 'NR%2{printf "%s ",$0;next;}1' | awk '{print $1,$2,$3,$8,$9}' | awk -v OFS='\t' '{min=$2;max=$2;for(i=2;i<=5;i++) {if ($i>max) max=$i;if($i<min) min=$i;} print $1,min,max}' > $frag
fi

# compute coverage & make track
bedtools genomecov -i $frag -bg -g $chrsz | slopBed -i stdin -g $chrsz -b 0 | bedClip stdin $chrsz $bdg
sort -k1,1 -k2,2n $bdg > $bdg_srt
bedGraphToBigWig $bdg_srt $chrsz $bigwig
rm -f $bdg 

