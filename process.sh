############################################################
# from narrow peak get summit 
############################################################
samples=(`cat samples.txt`)

chromsize=~/data/GENOME/hg38/hg38.chrom.sizes
cd .. 

# extend summit for each lib
for s in ${samples[@]}
do
    echo $s;
    s1=${s}.narrowPeak.gz
    zcat  $s1 | awk -v OFS='\t' '{$2=$2+$10;$3=$2+1}{print $1,$2,$3,$4,$8}' > ${s}.peak.summit.bed
    slopBed -i ${s}.peak.summit.bed -g $chromsize -b 250 | bedClip stdin $chromsize ${s}_1.peak.summit.ext.bed
done


# merge all peaks
cat  ${samples[@]/%/_1.peak.summit.ext.bed} > final.merged.bed
wc -l ./final.merged.bed ;#2466226
sort -k1,1 -k2,2n final.merged.bed | mergeBed -i stdin | sort -k1,1 -k2,2n final.merged.bed | mergeBed -i stdin > final;
mv final final.merged.bed;
wc -l ./final.merged.bed ; # 981621 

awk '{a=$3-$2+1;print a}' final.merged.bed > ./analysis/peak_width.txt 

# binary overlap matrix
for s in ${samples[@]}
do
    echo $s;
    intersectBed -a final.merged.bed  -b ${s}_1.peak.summit.ext.bed -c | awk '{print $4}' >$s.txt
done

echo  ${samples[@]} > ./analysis/binary.overlap.tsv
paste  ${samples[@]/%/.txt} >> ./analysis/binary.overlap.tsv


# peak scores per lib

for s in ${samples[@]}
do
    echo $s;
    awk '{print $5}'  ${s}_1.peak.summit.ext.bed   >${s}_pscore.txt
done

# max pileup/avg  RPM
i=223
for s in ${samples[@]}
do
    a=$(find ../data  -name "JYH_${i}_1_2_R*pileup*.bigwig")
    b=${s}_1.peak.summit.ext.bed
    bigWigAverageOverBed $a $b ${s}.tab 
    awk '{print $6}' ${s}.tab > ${s}_avgPileup.tab &&  rm ${s}.tab
done


############################################################
# from narrow peak all  summit info
############################################################
samples=(`cat samples.txt`)

cd .. 

# extend summit for each lib
for s in ${samples[@]}
do
    echo $s;
    s1=${s}.narrowPeak.gz
    zcat  $s1 | awk -v OFS='\t' '{$2=$2+$10;$3=$2+1}{print $1,$2,$3,$7,$8,$9}' > ${s}.peak.summit.bed
done


for s in ${samples[@]}
do
    echo $s;
    s1=${s}.final.bam
    zcat  ../$s1 | wc -l 
done


# tagalign to fragment bed, need modify NF = '\s

libs=($(cat ../../including_libs.txt))
for i in `seq 0 9`
do
    s=${libs[2*i+1]};    l=${libs[2*i]}
    echo $s;
    s1=/home/zhc268/data/outputs/tagAligns/${l}_R1.trim.PE2SE.nodup.tn5.tagAlign.gz
    zcat $s1 | awk 'NR%2{printf "%s ",$0;next;}1' | awk '{print $1,$2,$3,$8,$9}' | awk -v OFS='\t' '{min=$2;max=$2;for(i=2;i<=5;i++) {if ($i>max) max=$i;if($i<min) min=$i;} print $1,min,max}' > ${s}.fragments.bed#
done



# for each peak sumit, count #fragment
# sorted first
for s in ${samples[@]}
do
    echo $s;
    #s1=${s}.peak.summit.bed
    #sort -k1,1 -k2,2n  ../$s1> ${s}.sorted.peak.summit.bed
    s1=${s}.fragments.bed
    sort -k1,1 -k2,2n  ./$s1> ${s}.sorted.fragments.bed
done

for s in ${samples[@]}
do
    echo $s;
    s1=${s}.sorted.peak.summit.bed
    coverageBed -a ./$s1 -b  ${s}.sorted.fragments.bed -sorted  | awk '{print $7}' > ${s}.summit.frag.count.txt
done


# bined genome into 5k bins
bedtools makewindows -g ~/data/GENOME/hg38/hg38.chrom.sizes  -w 5000 > hg38.5K.windows.bed
bedtools makewindows -g ~/data/GENOME/hg38/hg38.chrom.sizes  -w 1000000 > hg38.1M.windows_a.bed
slopBed -i hg38.1M.windows_a.bed -g ~/data/GENOME/hg38/hg38.chrom.sizes -l -500000 -r 500000 | bedClip stdin $chromsize hg38.1M.windows_b.bed


sort -k1,1 -k2,2n hg38.1M.windows_a.bed > hg38.1M.windows_a.sorted.bed
sort -k1,1 -k2,2n hg38.1M.windows_b.bed > hg38.1M.windows_b.sorted.bed

paste  ${samples[@]/%/_hg38.1M.windows_a.sorted.bed.coverage.txt} > hg38.1M.windows_a.sorted.tsv
paste  ${samples[@]/%/_hg38.1M.windows_b.sorted.bed.coverage.txt} > hg38.1M.windows_b.sorted.tsv
rm ${samples[@]/%/_hg38.1M.windows_a.sorted.bed.coverage.txt} ${samples[@]/%/_hg38.1M.windows_b.sorted.bed.coverage.txt}
sed -i  's/ /\t/g' hg38.1M.windows_b.sorted.tsv 

# filter peaks reads

libs=($(cat ../../including_libs.txt))
for i in `seq 0 9`
do
    s=${libs[2*i+1]};    l=${libs[2*i]}
    echo $s;
    s1=/home/zhc268/data/outputs/tagAligns/${l}_R1.trim.PE2SE.nodup.tn5.tagAlign.gz
    subtractBed -a $s1 -b naive.merged.sorted.bed -A  > ${s}.pfilter.tn5.tagAlign
done




paste  ${samples[@]/%/_hg38.1M.windows_a.sorted.bed.pfilter.coverage.txt} > hg38.1M.windows_a.pfilter.sorted.tsv
paste  ${samples[@]/%/_hg38.1M.windows_b.sorted.bed.pfilter.coverage.txt} > hg38.1M.windows_b.pfilter.sorted.tsv
rm ${samples[@]/%/_hg38.1M.windows_a.sorted.bed.coverage.txt} ${samples[@]/%/_hg38.1M.windows_b.sorted.bed.coverage.txt}



############################################################
# process mappabiilty 
############################################################
bedtools makewindows -g $genome -w $swindow -s $step | awk -v OFS='\t' '{print $0,NR}'> hg38.500k_w_100k_s.bed

mapb="$(pwd)/data/GRCh38_36.bw"
bigWigAverageOverBed $mapb "../hg38.500k_w_100k_s.bed" "${mapb}.hgWindow.tab" -bedOut="${mapb}.hgWindow.bed"

# bedgraphtobigwig
bedSort "${mapb}.hgWindow.bed"  "${mapb}.hgWindow.bed"
awk -v OFS='\t' '{print $1,$2,$3,$5}' "${mapb}.hgWindow.bed" > tmp.bed
bedGraphToBigWig  tmp.bed $chromsize ${mapb}.hgWindow.bw

## make track
bigWigToBedGraph $mapb ${mapb}.bdg #bigwig to bdg 



mapb="$(pwd)/data/hg38_75bp/GRCh38_75.bw"
bigWigAverageOverBed $mapb "../beds/hg38.500k_w_100k_s.bed" "${mapb}.hgWindow.tab" -bedOut="${mapb}.hgWindow.bed"
macs2 pileup -i ${mapb}.hgWindow.bed -o ${mapb}.hgWindow.pileup.bed 


### use bedops
#https://www.biostars.org/p/226051/
awk '{ print $1"\t"$2"\t"$3"\tid-"NR"\t"$4; }' GRCh38_36.bw.bdg |sort-bed - > input.bed

awk  '{printf "%s\t%d\t%d\t%f\n", $1,($2+$3)/2,($2+$3)/2+1,$5}'  GRCh38_36.bw.hgWindow.bed > GRCh38_36.bw.hgWindow.bdg
bedGraphToBigWig  GRCh38_36.bw.hgWindow.bdg  $chromsize GRCh38_36.bw.hgWindow.bw


awk  '{printf "%s\t%d\t%d\t%f\n", $1,($2+$3)/2,($2+$3)/2+1,$5}'  GBM3_seq1_2.avg.cov.bed  > GBM3_seq1_2.avg.cov.bed.bdg
bedSort GBM3_seq1_2.avg.cov.bed.bdg GBM3_seq1_2.avg.cov.bed.bdg
bedGraphToBigWig GBM3_seq1_2.avg.cov.bed.bdg  $chromsize GRCh38_36.bw.hgWindow.bed.bw

############################################################
## denoise by gtak
############################################################
awk -v OFS='\t' 'BEGIN {print "@SQ","SN:HLA-DRB1*15:00:00","LN:11567","AH:*"; print "@RG","ID:GATKCopyNumber","SM:Gbm3";print "CONTIG","START","END","COUNT" } {printf "%s\t%d\t%d\t%d\n", $1,$2,$3,$5}' ./beds/GBM3_seq1_2.avg.cov.bed >./beds/GBM3_seq1_2.avg.cov.bed.tsv

awk -v OFS='\t' 'BEGIN {print "@SQ","SN:HLA-DRB1*15:00:00","LN:11567","AH:*"; print "@RG","ID:GATKCopyNumber","SM:mapb_36";print "CONTIG","START","END","COUNT" } {printf "%s\t%d\t%d\t%d\n", $1,$2,$3,1000*$5}' GRCh38_36.bw.hgWindow.bed > GRCh38_36.bw.hgWindow.tsv

gatk DenoiseReadCounts \
     -I ./beds/GBM3_seq1_2.avg.cov.bed.tsv \
     --count-panel-of-normals GRCh38_36.bw.hgWindow.tsv \
     --standardized-copy-ratios sample.standardizedCR.tsv \
     --denoised-copy-ratios sample.denoisedCR.tsv
 
