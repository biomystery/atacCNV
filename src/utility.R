require(data.table)
require(rtracklayer)
require(tidyverse)
require(ggplot2)
options(scipen=999)
require(DNAcopy)

set.seed(123)

### 
hg38.chrom.sizes <- "../data/hg38.chrom.sizes"
chr.levels <- paste0('chr',c(1:22,'X',"Y"))
chr.levels.noXY <- paste0('chr',c(1:22))

calcCov <- function(gcounts,p=peaks){

}

calcCNV <- function(genome.counts,l=libs){
  # norm to seq depth
  cnts.norm1 <- t(apply(cnts,1,function(x) x/libs$total.reads*1000000))
  
  # no normal sample for second normalization 
  
  # calculate 
}

readBinCounts <- function(fn="../data/genome.5k.counts.tsv",l=libs,
                          bed="../beds/hg38.5k.windows.sorted.bed",
                          p=peaks,
                          g=hg38.chrom.sizes){
  require(data.table)
  require(rtracklayer)
  
  genome.counts <- fread(fn,header = F)
  setDF(genome.counts)
  genome <- list(coverage=genome.counts[,seq(2,20,by = 2)],
                 counts= genome.counts[,seq(1,20,by = 2)])
  names(genome$coverage)<-names(genome$counts) <- l$name
  genome$locs <- import(bed,format = "bed",genome=g)
  genome$counts.cpm <- t(apply(genome$counts,1,function(x) x/l$total.reads*1000000))
  genome$FOM <- apply(genome$counts.cpm,2,function(x) x/median(x))
  #genome$cpm.coverage <- genome$counts.cpm/as.numeric(width(setdiff(genome$locs,p)))
  genome
}

##
write.bed.track <- function(df,fn,ext=250){
  df$start <- as.character(df$start-ext)
  df$end <- as.character(df$end + ext)
  df <- data.frame(df[,1:3],
                   name=paste0("p_",1:nrow(df)),
                   score = df[,5],
                   strand=".")
  write.table(df,file = fn,quote = F,col.names = F,row.names = F,sep = "\t")
}


saveGRtoBed <- function(gr,fn,nm="default"){
  nm=ifelse(nm=="default",c(rep(".",length(gr))),nm)
  df <- data.frame(chr=seqnames(gr),
                   start =start(gr)-1,
                   end = end(gr),
                   names=nm,
                   scores=c(rep(".",length(gr))),
                   strands=strand(gr)
  )
  write.table(df,file=fn,quote = F,sep="\t",row.names = F,col.names = F)
}


DEtest <- function(i,data=peaks.refined.extend.naive.merge.sorted.logcpm){
  #i <- 1
  l <- libs$name[i]
  lref <- libs$name[-i]
  
  m <- apply(data[,-i],1,mean)
  sdev<- apply(data[,-i],1,sd)
  target <- data[,i]
  #smoothScatter(peaks.refined.extend.naive.merge.sorted.logcpm[,i],m)
  
  lfc = target - m
  pval= sapply(1:length(m),function(x) 1-pnorm(target[x],mean = m[x],sd=sdev[x]))
  
  res <- data.frame(base.mean = m, 
                    base.sd=sdev,
                    lfc = lfc,
                    pval=pval,
                    padj.BH=p.adjust(pval,method = "BH"),
                    padj.Bonferroni=p.adjust(pval,method = "bonferroni"))
}


plot_pca <- function(x=peaks.refined.extend.naive.merge.sorted.cpm,...){
  pd.pca <- prcomp(t(x),center =T,...)#scale. = T,
  perct <- as.numeric(round(summary(pd.pca)$importance[2,1:3]*100))
  pd <- data.frame(
    pc1 = pd.pca$x[,1],
    pc2 = pd.pca$x[,2],
    pc3 = pd.pca$x[,3],
    sample=libs$name,
    tot.reads=libs$total.reads,
    tss=libs$tss.enrich
  )
  
  require(scatterD3)
  scatterD3(data=pd,pc1,pc2,lab = sample,col_var = tss,size_var = tot.reads,
            xlab = paste0("PC1: ",perct[1],"%"),
            ylab = paste0("PC2: ",perct[2],"%"),
            point_opacity = 0.5,hover_size = 4, hover_opacity = 1,lasso = T
            
  )
}

filter_gr <- function(gr,lvs=chr.levels.noXY){
  gr <- gr[seqnames(gr) %in% lvs]
  seqlevels(gr) <-lvs
  gr
}

scale_by_chr <- function(gr){
  s <- score(gr)
  for( l in seqlevels(gr)){
    id <- which(seqnames(gr) == l)
    s[id]=s[id]/median(s[id])
  }
  score(gr) <- s
  gr
}

plot_manh <- function(gr,...){
  x.brk<- table(seqnames(gr))
  x.brk <- cumsum(x.brk)
  gr<- gr[order(gr)] 
  require(RColorBrewer)
  col.map <- colorRampPalette(brewer.pal(name = "Paired",n=12))(length(x.brk))
  
  plot(1:length(gr),score(gr),pch=16,type='p',
       col=col.map[as.factor(seqnames(gr))],xaxt="n",
       xlab=NA,...)
  axis(1,at=x.brk,labels = names(x.brk),las=2)
  abline(v=x.brk,col='grey',lwd=1,lty=2)
  abline(h=1,col=1,lwd=1)
  abline(h=c(0.5,2),col='grey',lwd=1,lty=2)
  #points(idx,-log10(res$pval[idx]),cex=.5)
}
write.bed.seg <- function(seg,fn){
  df <- seg[,2:4]
  df[,2] <- df[,2]
  df[,3] <- df[,3]+1000000
  df <- data.frame(df[,1:3],
                   name=paste0("p_",1:nrow(df),"_",seg[1,1]),
                   score = seg$seg.mean,
                   strand=".")
  write.table(df,file = fn,quote = F,col.names = F,row.names = F,sep = "\t")
}
