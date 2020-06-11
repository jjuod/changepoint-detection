# Real data analysis: Chipseq

# devtools::install_github("tdhock/PeakError")
options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(PeakError)

ann.colors = c(noPeaks="#f6f4bf", peakStart="#ffafaf",
               peakEnd="#ff4c4c", peaks="#ff45cc")

# pick a dataset to start with
setwd("~/Documents/kiwis/changepointdet/chipseq/")

# problem = "H3K36me3_AM_immune/samples/monocyte/McGill0104/problems/chr9:92678796-133073060/"
# outprefix = "../drafts/changepoint-method/results-appl/McGill0104-chr9-"

# problem = "H3K36me3_AM_immune/samples/monocyte/McGill0104/problems/chr2:21178113-87668206/"
# outprefix = "../drafts/changepoint-method/results-appl/McGill0104-chr2-"

problem = "H3K36me3_AM_immune/samples/monocyte/McGill0104/problems/chr10:18024675-38818835/"
outprefix = "../drafts/changepoint-method/results-appl/McGill0104-chr10-"

chr = regmatches(outprefix, regexpr("chr[0-9]+", outprefix))
labels = read.table(paste0(problem, "labels.bed"))
colnames(labels) = c("chrom", "chromStart", "chromEnd", "annotation")
labels
bed = read.table(gzfile(paste0(problem, "coverage.bedGraph.gz")))

setwd("~/Documents/gitrep/changepoint-detection/")
source("CLEAN-algorithm.R")

# subset input to a window of 3x range(labels)
windowsize = max(labels$chromEnd) - min(labels$chromStart)
windowstart = min(labels$chromStart) - windowsize
bed = filter(bed, V2>=windowstart,
                V3<=windowstart+3*windowsize)
nrow(bed)
summary(bed)
max(bed$V3) - min(bed$V2)

# assuming no gaps further on, so need to check:
if(sum(bed$V3-lead(bed$V2), na.rm=T)>0){
  stop("Error: gaps detected in bed file")
}

# this might differ from the intended position if not all positions were seq'd:
actualstart = min(bed$V2)

# convert range format into annotations at each bp
bed = rep(bed$V4, times=bed$V3-bed$V2)
length(bed)

# downsample to get ~1k samples (by windowed mean)
DSfactor = round(length(bed)/1000, -3)
DSfactor
bedDS = data.frame(pos = round((seq_along(bed)+actualstart)/DSfactor)*DSfactor, cov=bed) %>%
  group_by(pos) %>%
  summarize(m = mean(cov))
nrow(bedDS)

# Note: bed coordinates seem to be in GRCh37!

# max lookback = max diff between reported peakStart and peakEnd
if("peakStart" %in% labels$annotation & "peakEnd" %in% labels$annotation){
  peakpos = which(labels$annotation=="peakStart")
  maxpeaklen = max(labels$chromEnd[peakpos+1] - labels$chromStart[peakpos])
  maxpeaklen = ceiling(maxpeaklen/DSfactor)
} else if("peaks" %in% labels$annotation){
  peakpos = which(labels$annotation=="peaks")
  maxpeaklen = max(labels$chromEnd[peakpos] - labels$chromStart[peakpos])
  maxpeaklen = ceiling(maxpeaklen/DSfactor)
} else {
  stop("Error: no peak annotations found")
}
maxpeaklen

plot(bedDS$pos/1e3, bedDS$m)


## Apply Algorithm 1:
ts = bedDS$m
pen = autoset_penalty(ts)
alg1 = shortsgd(ts, MAXLOOKBACK=maxpeaklen, PEN=pen, SD=sd(ts))
alg1$segs

# final theta estimate:
finaltheta = alg1$wt[length(alg1$wt)]
finaltheta

# we are interested only in segments of increased mean:
alg1incr = alg1$segs[alg1$segs[,3]>finaltheta,,drop=F]

alg1res = data.frame(chrom=chr,
                     chromStart=alg1incr[,1]*DSfactor + actualstart,
                     chromEnd=alg1incr[,2]*DSfactor + actualstart, 
                     row.names = NULL)
alg1res
labels

errs = PeakError(alg1res, labels)

ggplot()+
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

# save output
save(alg1, file=paste0(outprefix, "alg1.RData"))
write.table(errs, paste0(outprefix, "alg1.tsv"), quote=F, row.names=F, sep="\t")
ggsave(paste0(outprefix, "alg1.png"), width=14, height=8, units="cm")


## Apply Algorithm 2:
ts = bedDS$m
pen = autoset_penalty(ts)
alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen, PEN=pen, PEN2=pen, BURNIN=2, SD=sd(ts), prune=2)
alg2$segs

# we are interested only in segments of increased mean:
alg2incr = alg2$segs[alg2$segs[,3]>finaltheta,, drop=F]

alg2res = data.frame(chrom=chr,
                     chromStart=alg2incr[,1]*DSfactor + actualstart,
                     chromEnd=alg2incr[,2]*DSfactor + actualstart, 
                     segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                     row.names = NULL)
alg2res
labels

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
  geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3),
               data=alg2res[alg2res$segtype=="S",], color="deepskyblue", size=3)+
  xlab("position on chromosome, kbp") + ylab("coverage")
if("N" %in% alg2res$segtype){
  p2 = p2 + geom_segment(aes(chromStart/1000, -4, xend=chromEnd/1000, yend=-4),
               data=alg2res[alg2res$segtype=="N",], color="purple", size=3)
}
p2

# save output
write.table(errs, paste0(outprefix, "alg2.tsv"), quote=F, row.names=F, sep="\t")
save(alg2, file=paste0(outprefix, "alg2.RData"))
ggsave(paste0(outprefix, "alg2.png"), width=14, height=8, units="cm")
