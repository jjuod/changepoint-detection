# Real data analysis: Chipseq UCI

options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(anomaly)
library(not)
library(PeakSegDisk)

setwd("~/Documents/gitrep/changepoint-detection/")
source("MAIN-algorithm.R")

# Setup input-output dirs for this window
datadir = "~/Documents/kiwis/changepointpaper/chipseq/H3K36me3_AM_immune/samples/monocyte/McGill0104/problems/"
probname = "chr12:37856694-109373470/"
outprefix = paste0("../drafts/changepoint-method/results-appl/McGill0104-", substr(probname, 1, 7), "-")

ann.colors = c(noPeaks="#fffeef", peakStart="#ffafaf",
               peakEnd="#ff4c4c", peaks="#baa3ff")


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
if(maxpeaklen<2){
  stop("Error: max segment length seems very short")
} else if(maxpeaklen>500){
  stop("Error: max segment length close to data size")
}
print(sprintf("Using max segment length: %d points", maxpeaklen))

plot(bedDS$pos/1e3, bedDS$m)


## DETECTION RUNS:
## Shared parameters
ts = bedDS$m
finaltheta = median(ts)
pen = autoset_penalty(ts)
robust_sd = mad(diff(bedDS$m))
ts_norm = (ts-finaltheta)/robust_sd
tsdf = data.frame(chrom=chr, chromStart=as.integer(bedDS$pos),
                  chromEnd=as.integer(bedDS$pos+DSfactor),
                  count=as.integer(round(bedDS$m)))

## Apply standard epidemic chp detection (anomaly)
res.anom = capa.uv(ts_norm, beta=10*pen, beta_tilde=10*pen, type="mean",
                   min_seg_len=2, max_seg_len=maxpeaklen, transform=identity)
# Extract segments, and convert effect sizes to raw means (similar to Alg2)
res.anom.c = collective_anomalies(res.anom)[,1:3]
res.anom.c$mean.change = apply(res.anom.c, 1, function(x) mean(ts[x[1]:x[2]]))
res.anom.p = point_anomalies(res.anom)
if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
  res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location,
                          mean.change=ts[res.anom.p$location])
  res.anom.c = rbind(res.anom.c, res.anom.p)
}
# we are interested only in segments of increased mean:
anomincr = res.anom.c[res.anom.c[,3]>finaltheta,, drop=F]
anomres = data.frame(chrom=chr,
                     chromStart=anomincr[,1]*DSfactor + actualstart,
                     chromEnd=anomincr[,2]*DSfactor + actualstart+DSfactor,
                     segtype="S",
                     row.names = NULL)
print(anomres)
# save output
save(res.anom, file=paste0(outprefix, "anom.RData"))

## Apply NOT (Narrowest-over-threshold detector)
res.not = not(ts_norm, method="not", contrast="pcwsConstMean")
res.not = features(res.not, method="threshold", th=0.4*pen)$cpt
res.not = cbind(c(1, res.not),
                c(res.not-1, length(ts)))
res.not = cbind(res.not, apply(res.not, 1, function(x){
  if(any(is.na(x))) { NA }
  else { mean(ts[x[1]:x[2]]) }
}))
notincr = res.not[res.not[,3]>finaltheta,, drop=F]
notres = data.frame(chrom=chr,
                    chromStart=notincr[,1]*DSfactor + actualstart,
                    chromEnd=notincr[,2]*DSfactor + actualstart+DSfactor, 
                    segtype=ifelse(notincr[,3]>finaltheta+robust_sd, "S", "S2"),
                    row.names = NULL)
print(notres)
# save output
save(res.not, file=paste0(outprefix, "not.RData"))

## Apply PeakSegFPOP
# use the penalty chosen in Broad data
res.fpop = PeakSegFPOP_df(tsdf, 10e3, "/tmp/")$segments
fpopincr = filter(data.frame(res.fpop), status=="peak", mean>finaltheta) %>%
  mutate(segtype = "S")
print(fpopincr)
# save output
save(res.fpop, file=paste0(outprefix, "fpop.RData"))


## Apply Algorithm 2:
alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen,
                          PEN=10*pen, PEN2=10*pen, SD=robust_sd, prune=2)
print(alg2$segs)

# we are interested only in segments of increased mean:
alg2incr = alg2$segs[alg2$segs[,3]>finaltheta,, drop=F]
alg2res = data.frame(chrom=chr,
                     chromStart=alg2incr[,1]*DSfactor + actualstart,
                     chromEnd=alg2incr[,2]*DSfactor + actualstart+DSfactor, 
                     segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                     row.names = NULL)
print(alg2res)
print(labels)
# save output
save(alg2, file=paste0(outprefix, "alg2.RData"))

## Plot all results:
toplotL = data.frame(label=c("anomaly", "not", "PeakSeg", "proposed"),
                     # x=max(alg2res$chromEnd)/1000+5,
                     x = 63260,
                     y=c(-4, -7, -10, -13))
p2 = ggplot()+
  geom_rect(aes(xmin=chromStart/1000, xmax=chromEnd/1000,
                ymin=-2.5, ymax=-0.5, fill=annotation),
            # , linetype=fn.status, size=fp.status),
            data=labels, color="black") +
  # geom_point(aes(x=pos/1000, y=m), alpha=0.5, data=bedDS) +
  geom_point(aes(x=pos/1000, y=m), shape=21, col="black", fill="grey50", stroke=0.5, data=bedDS) +
  geom_segment(aes(chromStart/1000, -4, xend=chromEnd/1000, yend=-4, color=segtype, size=segtype),
               data=anomres)
if(nrow(notres)>0){
  p2 = p2 +
    geom_segment(aes(chromStart/1000, -7, xend=chromEnd/1000, yend=-7, color=segtype, size=segtype),
                 data=notres) 
}
if(nrow(fpopincr)>0){
  p2 = p2 +
    geom_segment(aes(chromStart/1000, -10, xend=chromEnd/1000, yend=-10, color=segtype, size=segtype),
                 data=fpopincr)
}
if(nrow(alg2res)>0){
  p2 = p2 +
    geom_segment(aes(chromStart/1000, -13, xend=chromEnd/1000, yend=-13, color=segtype, size=segtype),
                 data=alg2res)
}
p2 = p2 +  
  geom_text(aes(x=x, y=y, label=label), data=toplotL, hjust=0) +
  # scale_linetype_manual(values=c("false negative"="dotted", correct="solid"))+
  scale_size_manual(values=c("false positive"=2, correct=0.7, "N"=7, "S"=6, "S2"=6))+
  scale_fill_manual(values=ann.colors, breaks=names(ann.colors), name="annotations")+
  # scale_color_manual(values=c("S"="blue", "N"="gold"), name="detections", labels=c("nuisance", "signal")) +
  scale_color_manual(values=c("S"="blue", "N"="gold", "S2"="#D8EFFF"), name="detections",
                     labels=c("nuisance", "signal", "sign. low")) +
  theme_bw()+ theme(legend.key = element_rect(), text=element_text(size=12))+
  guides(fill=guide_legend(order=1),
         # linetype=guide_legend(order=3, override.aes=list(fill="white")),
         color=guide_legend(order=2, override.aes=list(size=7)),
         size="none")+
  # xlim(c(62460, 62720)) + 
  xlab(paste0("position on chromosome ", substring(chr, 4), ", kbp")) + ylab("coverage")

print(p2)

# ggsave(paste0(outprefix, "alg2.png"), width=18, height=10, units="cm")
ggsave(paste0(outprefix, "alg2.eps"), width=18, height=10, units="cm", device=cairo_ps)


## Calculate fit statistics (BIC/SIC)
# BIC
# Calculates likelihood assuming a normal model with current finaltheta and SD.
# For nuisance segments, uses only the points N\S.
#  - df: a dataframe of only the segments detected.
#  - Requires "start" and "end" columns in ts positions, and "segtype" with S or N
getBIC = function(df, ts, mu0, sigma0){
  nparams = 0
  logl = bglogl = 0
  if(!"segtype" %in% colnames(df)){
    df$segtype = "S"
  }
  
  # first, mask any S segment points so that they would be excluded
  # from nuisance and bg likelihoods later
  # and mask N points for excluding from bg
  isBg = rep(T, length(ts))
  isNuis = rep(T, length(ts))
  for(j in 1:nrow(df)){
    isBg[df$start[j]:df$end[j]] = F
    if(df$segtype[j]=="S"){
      isNuis[df$start[j]:df$end[j]] = F
    }
  }

  for(i in 1:nrow(df)){
    start = df$start[i]
    end = df$end[i]
    nparams = nparams+1
    
    # now, treat S segments as usual
    if(df$segtype[i]=="S"){
      segx = ts[start:end]
    } else {
      # for N segments, first exclude any overlapping S points
      nuismask = logical(length(ts))
      nuismask[start:end] = T
      segx = ts[nuismask & isNuis]
    }
    logl = logl - 2*sum(dnorm(segx, mean(segx), sigma0, log=T))
  }
  
  # deal with bg points
  bgpoints = ts[isBg]
  bglogl = bglogl - 2*sum(dnorm(bgpoints, mu0, sigma0, log=T))

  # num parameters = num segments x3 (for mean and two endpoints in each) + 2
  klogn = (3*nrow(df)+2)*log(length(ts))
  cat(sprintf("Log-likelihood of background points: %.1f\n", bglogl))
  cat(sprintf("Log-likelihood of non-bg points: %.1f\n", logl))
  cat(sprintf("penalty k log(n): %.1f\n", klogn))
  return(logl+bglogl+klogn)
}

bic.anom = getBIC(res.anom.c, ts, finaltheta, robust_sd)

# define "background" segments
res.not.segs = data.frame(res.not[res.not[,3]-finaltheta>0, ])
colnames(res.not.segs) = c("start", "end", "mean")
bic.not = getBIC(res.not.segs, ts, finaltheta, robust_sd)

# extract background and convert back to ts coordinates
res.fpop.segs = filter(res.fpop, status=="peak") %>% arrange(chromStart)
res.fpop.segs$start = match(res.fpop.segs$chromStart, bedDS$pos)
res.fpop.segs$end = match(res.fpop.segs$chromEnd, bedDS$pos)
bic.fpop = getBIC(res.fpop.segs, ts, finaltheta, robust_sd)

res.alg2.segs = data.frame(alg2$segs)
colnames(res.alg2.segs) = c("start", "end", "mean", "segtype")
res.alg2.segs$segtype = ifelse(res.alg2.segs$segtype==1, "S", "N")
bic.alg2 = getBIC(res.alg2.segs, ts, finaltheta, robust_sd)

bic.anom
bic.not
bic.fpop
bic.alg2
