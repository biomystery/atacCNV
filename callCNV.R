#!/usr/bin/env Rscript
require(optparse)
# parse options  ----------------------------------------------------------
# Inputs: 
option_list = list(
  make_option(c("-c", "--ctrl"), type="character", default="./data/hg38_75bp/GRCh38_75.avg_cov.w1000K_s500K.bed", 
              help="ctrl bed file [default= %default]", metavar="character"),
  make_option(c("-t", "--treat"), type="character", default="out.txt", 
              help="treat bed file", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="../data/hg38.chrom.sizes", 
              help="genome size file [default= %default]", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-u", "--upper"), type="double", default=2, 
              help="upper boundary in trimming log2Ratio [default= %default]", metavar="double"),
  make_option(c("-m", "--mapb_th"), type="double", default=0.75, 
              help="mapbility threshold for regions [default= %default]", metavar="double"),
  make_option(c("-o", "--output_dir"), type="character", default="./data/output/", 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-l", "--lower"), type="double", default=-1.5, 
              help="lower boundary in trimming log2Ratio [default= %default]", metavar="double")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check if treat exist 
if (!file.exists(opt$treat)){
  print_help(opt_parser)
  stop("(treat bed file not exist!).n", call.=FALSE)
}

if (!file.exists(opt$ctrl)){
  print_help(opt_parser)
  stop("(control bed file not exist!).n", call.=FALSE)
}

if (is.null(opt$sample)){
  opt$sample = sub(".bed","",sub(".*/","",opt$treat))
}

# Intro  ------------------------------------------------------------------
suppressMessages(source("./ultility.R"))

# ctrl 
dat.ctrl <- import(opt$ctrl,format = "bed",genome=opt$genome)
dat.ctrl <- filter_gr(dat.ctrl)

# treat
dat.treat <- import(opt$treat,format = "bed",genome=hg38.chrom.sizes)
dat.treat <- filter_gr(dat.treat)

# Check ctrl will be used to weight the coverage at each location 
if (!all.equal(length(match(dat.ctrl,dat.treat)), length(dat.ctrl), length(dat.treat))){
  print_help(opt_parser)
  stop("(contrl and treat are not the same resolution!).n", call.=FALSE)
}

# assign median to low mapb regions
idx.filter <- score(dat.ctrl)>=opt$mapb_th
dat.treat.mod <- dat.treat
score(dat.treat.mod)[!idx.filter] <- median(score(dat.treat.mod)[idx.filter])
raw.t <- score(dat.treat.mod)/median(score(dat.treat.mod))
summary(raw.t)
raw.c <- score(dat.ctrl)
raw.c[!idx.filter] <- median(raw.c[idx.filter]) 

# standardizatin & scaling 
ratio.t <- raw.t/(raw.c+1)*2 # adjusted coverage by mapb 
ratio.t <- ratio.t/median(ratio.t) # ratio of sample's median 
score(dat.treat.mod) <- ratio.t
ratio.t <- score(scale_by_chr(dat.treat.mod)) # scaled by each chr's median
ratio.t.log2 <- log2(ratio.t) # convert to log2Ratio 

# (optional) trim the log2Ratio (upper, lower )
ratio.t.log2[ratio.t.log2> opt$upper] <- opt$upper
ratio.t.log2[ratio.t.log2 < opt$lower] <- opt$lower

# create CNA object 
CNA.object <- CNA(cbind(ratio.t.log2),
                  as.character(seqnames(dat.treat.mod)),
                  start(dat.treat.mod),
                  data.type="logratio",sampleid=opt$sample)

# smooth track by removing outliner and smooth to neighbours 
smoothed.CNA.object <- smooth.CNA(CNA.object)

# detect segments
#segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
sdundo.CNA.object <- segment(smoothed.CNA.object,
                             undo.splits="sdundo",
                             undo.SD=3,verbose=1)
# plots 
if(F){
pdf(file=paste0(opt$output_dir,opt$sample,".pdf"),width = 12,height = 7)
plot(segment.smoothed.CNA.object, plot.type="w",ylim=c(-2,2))
plot(segment.smoothed.CNA.object, plot.type="s",ylim=c(-2,2))
plot(segment.smoothed.CNA.object, plot.type="p",ylim=c(-2,2))
dev.off()}

pdf(file=paste0(opt$output_dir,opt$sample,"_SD.pdf"),width = 12,height = 7)
plot(sdundo.CNA.object, plot.type="w",ylim=c(-2,2))
plot(sdundo.CNA.object, plot.type="s",ylim=c(-2,2))
plot(sdundo.CNA.object, plot.type="p",ylim=c(-2,2))
dev.off()

# export segments by chr and ranking 
write_rds(list(opt=opt,
               smoothed.CNA.object =smoothed.CNA.object,
               sdundo.CNA.object=sdundo.CNA.object),
          path = paste0(opt$output_dir,opt$sample,".rds"))

write_csv(sdundo.CNA.object$output,
          path = paste0(opt$output_dir,opt$sample,"_SD.csv"))

if(sum(abs(sdundo.CNA.object$output$seg.mean)>=1)>0)
  write.bed.seg(sdundo.CNA.object$output%>% filter(abs(seg.mean)>=1),
                fn = paste0(opt$output_dir,opt$sample,".cnv.bed"))
# optional cmds -----------------------------------------------------------

#glFrequency(segment.smoothed.CNA.object)
#plotSample(sdundo.CNA.object,ylim=c(-3,3))
#zoomIntoRegion(sdundo.CNA.object, chrom="chr15",
#               maploc.start = 57851725,maploc.end = 64357541,
#               sampleid="GBM3",
#               pt.pch = 16,pt.cex = .5,ylim=c(-3,3))
#abline(h=0,col='grey')


if(F){
  zoomIntoRegion(sdundo.CNA.object, chrom="chr8",
                 #maploc.start = 67068340,maploc.end = 86000925,
                 sampleid="GBM3",
                 pt.pch = 16,pt.cex = .5,ylim=c(-2,2))
  abline(h=c(-1,1),col='grey')
}


if(F){
  png(file="./atacCNV/figs/GBM3_cmp_smooth.png")
  par(mfrow=c(2,1))
  plot(CNA.object$GBM3,ylim=c(-2.5,8),cex=0.25)
  plot(smoothed.CNA.object$GBM3,ylim=c(-2.5,8),cex=0.25)
  par(mfrow=c(1,1))
  dev.off()
}

### case 1: mapb = 0
if(F){
  hist(score(dat.ctrl),nclass = 50)
  summary(score(dat.ctrl))
  summary(score(dat.treat))
  
  #plot(CNA.object, plot.type="")
  plot(segment.smoothed.CNA.object, plot.type="p")
  
  
  plot(segment.smoothed.CNA.object, plot.type="w")
  plot(sdundo.CNA.object,plot.type="p")
  plot(sdundo.CNA.object,plot.type="w")
  plot(sdundo.CNA.object,plot.type="s")
  # check segmentation
  
  par(mfrow=c(2,1))
  plot(segment.smoothed.CNA.object, plot.type="w")
  plot(sdundo.CNA.object,plot.type="w")
  par(mfrow=c(1,1))  
}
