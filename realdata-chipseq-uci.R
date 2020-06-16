# Real data analysis: Chipseq

# devtools::install_github("tdhock/PeakError")
options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(PeakError)

setwd("~/Documents/gitrep/changepoint-detection/")
source("CLEAN-algorithm.R")

ann.colors = c(noPeaks="#f6f4bf", peakStart="#ffafaf",
               peakEnd="#ff4c4c", peaks="#ff45bc")

outdir = "../drafts/changepoint-method/results-appl/McGill0104-"

# pick a dataset to start with
datadir = "~/Documents/kiwis/changepointdet/chipseq/H3K36me3_AM_immune/samples/monocyte/McGill0104/problems/"
allprobs = paste0(list.dirs(datadir, full.names=F), "/")[-1]  # (first entry is ./)

# Note: bed coordinates are in GRCh37!

# one problem has a very large window which results in very short peaks after downsampling.
# therefore, we split it into two beforehand:
if(!dir.exists(paste0(datadir, "chr10:63000000-125869472/"))){
  labels = read.table(paste0(datadir, "chr10:51448845-125869472/", "labels.bed"))
  labels
  labels1 = filter(labels, V2<=63e6)
  labels2 = filter(labels, V2>63e6)
  dir1 = paste0(datadir, "chr10:51448845-63000000/")
  dir2 = paste0(datadir, "chr10:63000000-125869472/")
  
  # copy the files
  dir.create(dir1)
  dir.create(dir2)
  write.table(labels1, paste0(dir1, "labels.bed"), col.names=F, sep="\t", row.names=F, quote=F)
  write.table(labels2, paste0(dir2, "labels.bed"), col.names=F, sep="\t", row.names=F, quote=F)
  file.copy(paste0(datadir, "chr10:51448845-125869472/coverage.bedGraph.gz"), dir1)
  file.copy(paste0(datadir, "chr10:51448845-125869472/coverage.bedGraph.gz"), dir2)
  
  # update search path list
  allprobs = paste0(list.dirs(datadir, full.names=F), "/")[-1]
}
# exclude the dir that was split
allprobs = allprobs[allprobs!="chr10:51448845-125869472/"]

analyzeProblem <- function(probname){
  print(sprintf("Working on file %s", probname))
  outprefix = paste0(outdir, substr(probname, 1, 7), "-")
  
  chr = regmatches(outprefix, regexpr("chr[0-9]+", outprefix))
  labels = read.table(paste0(datadir, probname, "labels.bed"))
  colnames(labels) = c("chrom", "chromStart", "chromEnd", "annotation")
  print(labels)
  bed = read.table(gzfile(paste0(datadir, probname, "coverage.bedGraph.gz")))
  
  # subset input to a window of 3x range(labels)
  windowsize = max(labels$chromEnd) - min(labels$chromStart)
  windowstart = min(labels$chromStart) - windowsize
  bed = filter(bed, V2>=windowstart,
               V3<=windowstart+3*windowsize)
  print(nrow(bed))

  # assuming no gaps further on, so need to check:
  if(sum(bed$V3-lead(bed$V2), na.rm=T)>0){
    stop("Error: gaps detected in bed file")
  }
  
  # this might differ from the intended position if not all positions were seq'd:
  actualstart = min(bed$V2)
  
  # convert range format into annotations at each bp
  bed = rep(bed$V4, times=bed$V3-bed$V2)
  print(length(bed))
  
  # downsample to get ~1k samples (by windowed mean)
  DSfactor = round(length(bed)/1000, -2)
  print(sprintf("Downsampling %d x", DSfactor))
  bedDS = data.frame(pos = round((seq_along(bed)+actualstart)/DSfactor)*DSfactor, cov=bed) %>%
    group_by(pos) %>%
    summarize(m = mean(cov))
  print(nrow(bedDS))

 # fixed max lookback of 50 kbp
  maxpeaklen = round(50000/DSfactor)
  if(maxpeaklen<3){
    stop("Error: max segment length seems very short")
  } else if(maxpeaklen>500){
    stop("Error: max segment length close to data size")
  }
  print(sprintf("Using max segment length: %d points", maxpeaklen))
  
  # plot(bedDS$pos/1e3, bedDS$m)
  
  # get confirmation before continuing
  confirm = readline("Continue?")
  
  ## Apply Algorithm 1:
  ts = bedDS$m
  pen = autoset_penalty(ts)
  alg1 = shortsgd(ts, MAXLOOKBACK=maxpeaklen, PEN=pen, SD=sd(ts))
  print(alg1$segs)
  
  # final theta estimate:
  finaltheta = alg1$wt[length(alg1$wt)]
  print(finaltheta)
  
  # we are interested only in segments of increased mean:
  alg1incr = alg1$segs[alg1$segs[,3]>finaltheta,,drop=F]
  
  # empty dataframe with a "segment" at 0-1 to allow plotting
  alg1res = data.frame(chrom=chr, chromStart=actualstart, chromEnd=actualstart+1)
  if(nrow(alg1incr)>0){
    alg1res = data.frame(chrom=chr,
                         chromStart=alg1incr[,1]*DSfactor + actualstart,
                         chromEnd=alg1incr[,2]*DSfactor + actualstart+1, 
                         row.names = NULL)  
  }
  print(alg1res)
  print(labels)
  
  errs = PeakError(alg1res, labels)
  
  p1 = ggplot()+
    geom_rect(aes(xmin=chromStart/1000, xmax=chromEnd/1000,
                  ymin=-2, ymax=0, fill=annotation,
                  linetype=fn.status, size=fp.status),
              data=errs, color="black") +
    geom_point(aes(x=pos/1000, y=m), size=0.5, data=bedDS) +
    scale_linetype_manual(values=c("false negative"="dotted", correct="solid"))+
    scale_size_manual(values=c("false positive"=2, correct=0.7))+
    scale_fill_manual(values=ann.colors, breaks=names(ann.colors))+
    theme_bw()+
    guides(fill=guide_legend(order=1),
           linetype=guide_legend(order=2, override.aes=list(fill="white")),
           size=guide_legend(order=3, override.aes=list(fill="white")))+
    geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3),
                 data=alg1res, color="deepskyblue", size=3)+
    xlab("position on chromosome, kbp") + ylab("coverage")
  print(p1)
  
  # save output
  save(alg1, file=paste0(outprefix, "alg1.RData"))
  write.table(errs, paste0(outprefix, "alg1.tsv"), quote=F, row.names=F, sep="\t")
  ggsave(paste0(outprefix, "alg1.png"), width=14, height=8, units="cm")
  
  
  ## Apply Algorithm 2:
  ts = bedDS$m
  pen = autoset_penalty(ts)
  alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen, PEN=pen, PEN2=pen, BURNIN=2, SD=sd(ts), prune=2)
  print(alg2$segs)
  
  # we are interested only in segments of increased mean:
  alg2incr = alg2$segs[alg2$segs[,3]>finaltheta,, drop=F]
  
  alg2res = data.frame(chrom=chr,
                       chromStart=alg2incr[,1]*DSfactor + actualstart,
                       chromEnd=alg2incr[,2]*DSfactor + actualstart+1, 
                       segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                       row.names = NULL)
  print(alg2res)
  print(labels)
  
  errs = PeakError(alg2res[alg2res$segtype=="S",], labels)
  
  p2 = ggplot()+
    geom_rect(aes(xmin=chromStart/1000, xmax=chromEnd/1000,
                  ymin=-2, ymax=0, fill=annotation,
                  linetype=fn.status, size=fp.status),
              data=errs, color="black") +
    geom_point(aes(x=pos/1000, y=m), size=0.5, data=bedDS) +
    scale_linetype_manual(values=c("false negative"="dotted", correct="solid"))+
    scale_size_manual(values=c("false positive"=2, correct=0.7))+
    scale_fill_manual(values=ann.colors, breaks=names(ann.colors))+
    theme_bw()+
    guides(fill=guide_legend(order=1),
           linetype=guide_legend(order=2, override.aes=list(fill="white")),
           size=guide_legend(order=3, override.aes=list(fill="white")))+
    xlab("position on chromosome, kbp") + ylab("coverage")
  if("S" %in% alg2res$segtype){
    p2 = p2 + geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3),
                           data=alg2res[alg2res$segtype=="S",], color="deepskyblue", size=3)
  }
  if("N" %in% alg2res$segtype){
    p2 = p2 + geom_segment(aes(chromStart/1000, -4, xend=chromEnd/1000, yend=-4),
                           data=alg2res[alg2res$segtype=="N",], color="purple", size=3)
  }
  print(p2)
  
  # save output
  write.table(errs, paste0(outprefix, "alg2.tsv"), quote=F, row.names=F, sep="\t")
  save(alg2, file=paste0(outprefix, "alg2.RData"))
  ggsave(paste0(outprefix, "alg2.png"), width=14, height=8, units="cm")
}

for(problem in allprobs){
  try(analyzeProblem(problem), outFile=paste0(substr(problem, 1, 6), "analysisErrors.log"))
}

# Calculate accuracy metrics
datadir = "~/Documents/kiwis/changepointdet/chipseq/H3K36me3_AM_immune/samples/monocyte/McGill0104/problems/"
allprobs = paste0(list.dirs(datadir, full.names=F), "/")[-1]  # (first entry is ./)
# exclude the dir that was split
allprobs = allprobs[allprobs!="chr10:51448845-125869472/"]

getErr <- function(probname){
  print(sprintf("Working on file %s", probname))
  outprefix = paste0(outdir, substr(probname, 1, 7), "-")
  
  chr = regmatches(outprefix, regexpr("chr[0-9]+", outprefix))
  labels = read.table(paste0(datadir, probname, "labels.bed"))
  colnames(labels) = c("chrom", "chromStart", "chromEnd", "annotation")
  print(labels)
  
  # get Bed properties
  bed = read.table(gzfile(paste0(datadir, probname, "coverage.bedGraph.gz")))
  
  # subset input to a window of 3x range(labels)
  windowsize = max(labels$chromEnd) - min(labels$chromStart)
  windowstart = min(labels$chromStart) - windowsize
  bed = filter(bed, V2>=windowstart,
               V3<=windowstart+3*windowsize)
  print(nrow(bed))
  
  # this might differ from the intended position if not all positions were seq'd:
  actualstart = min(bed$V2)
  
  # convert range format into annotations at each bp
  bed = rep(bed$V4, times=bed$V3-bed$V2)
  print(length(bed))
  
  # downsample to get ~1k samples (by windowed mean)
  DSfactor = round(length(bed)/1000, -2)
  print(sprintf("Downsampling %d x", DSfactor))
  
  # load Alg 2 results
  alg2res = NULL
  if(file.exists(paste0(outprefix, "alg2.RData"))){
    load(paste0(outprefix, "alg1.RData"))
    # final theta estimate:
    finaltheta = alg1$wt[length(alg1$wt)]
    print(finaltheta)
    
    load(paste0(outprefix, "alg2.RData"))
    
    # we are interested only in segments of increased mean:
    alg2incr = alg2$segs[alg2$segs[,3]>finaltheta,, drop=F]
    alg2res = data.frame(chrom=chr,
                         chromStart=alg2incr[,1]*DSfactor + actualstart,
                         chromEnd=alg2incr[,2]*DSfactor + actualstart+1, 
                         segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                         row.names = NULL)
  } else {
    # or create an "empty" df if Alg2 didn't detect anything
    alg2res = data.frame("chrom"=chr, "chromStart"=1, "chromEnd"=2)
  }
  print(alg2res)
  
  # without counting N as detections:
  errs = PeakError(alg2res[alg2res$segtype=="S",], labels)
  # including N as detections:
  # (first do only N and then correct for S)
  errs2 = PeakError(alg2res[alg2res$segtype=="N",], labels)
  # i.e.: if either S or N produced an FP -> FP
  errs2$fp.status = ifelse(errs2$fp.status=="correct" & errs$fp.status=="correct", "correct", "false positive")
  # i.e.: if either S or N, or both, detected -> not FN
  errs2$fn.status = ifelse(errs2$fn.status=="correct" | errs$fn.status=="correct", "correct", "false negative")
  errs2$status = ifelse(errs2$fp.status=="correct" & errs2$fn.status=="correct", "correct",
                        ifelse(errs2$fp.status=="false positive", "false positive", "false negative"))
  return(list("onlyS" = errs, "withN" = errs2))
}

allerrs = data.frame()
allerrsN = data.frame()
for(prob in allprobs){
  err = getErr(prob)
  allerrs = bind_rows(allerrs, err$onlyS)
  allerrsN = bind_rows(allerrsN, err$withN)
}

# analyse errors without N 
TP = sum(allerrs$tp)
FP = sum(allerrs$fp)
TN = sum(allerrs$possible.fp) - FP
FN = sum(allerrs$possible.tp) - TP
sens = TP / (TP+FN)
spec = TN / (TN+FP)
acc = (TP+TN)/(TP+TN+FP+FN)
sens;spec;acc

# analyse errors including N
TP = sum(allerrsN$tp)
FP = sum(allerrsN$fp)
TN = sum(allerrsN$possible.fp) - FP
FN = sum(allerrsN$possible.tp) - TP
sens = TP / (TP+FN)
spec = TN / (TN+FP)
acc = (TP+TN)/(TP+TN+FP+FN)
sens;spec;acc

###
group_by(err, prob.dir) %>%
  summarize(merr = mean(errors), min(errors), max(errors)) %>%
  summarize(sum(merr))
(111-25.5)/111

# Alternative targets:
# H3K4me3 (ENCODE or immune) does not seem useful for our method: peaks are all narrow and background is even
# H3K27ac (some) sort of intermediate, could use
# ATAC has very flat bg (0) and lots of small peaks which aren't marked as peaks, so not easy to interpret
