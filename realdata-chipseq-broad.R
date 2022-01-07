# Real data analysis: Chipseq Broad

options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(anomaly)
library(not)
library(PeakSegDisk)

setwd("~/Documents/gitrep/changepoint-detection/")
source("MAIN-algorithm.R")

# Note: first use prep-broad-chipseq.sh to convert UCSC raw data into bedGraph 

indir = "~/Documents/kiwis/changepointpaper/chipseq/BROAD/"
bedC = read.table(paste0(indir, "filtered_control.bedGraph"))
bedH = read.table(paste0(indir, "filtered_h3k27ac.bedGraph"))

plot(bedC$V2, bedC$V4)
plot(bedH$V2, bedH$V4)

bedCf = filter(bedC, V2>120100e3, V3<=120700e3)
bedHf = filter(bedH, V2>120100e3, V3<=120700e3)

chr = "chr1"

# align different start positions
actualstart = max(min(bedCf$V2), min(bedHf$V2))
actualend = min(max(bedCf$V2), max(bedHf$V2))
bedCf = filter(bedCf, V2>=actualstart, V3<=actualend)
bedHf = filter(bedHf, V2>=actualstart, V3<=actualend)

# sections with no reads appear as gaps - need to fill those in with 0 counts
toappend = data.frame()
for (i in 2:nrow(bedCf)){
  if(bedCf$V2[i]-bedCf$V3[i-1]>0){
    neededpos = seq(bedCf$V3[i-1], bedCf$V2[i]-1, 25)
    toappend = bind_rows(toappend, data.frame(V1=chr, V2=neededpos, V3=neededpos+25, V4=0))
  }
}
bedCf = bind_rows(bedCf, toappend) %>% arrange(V2)
# same for treated
toappend = data.frame()
for (i in 2:nrow(bedHf)){
  if(bedHf$V2[i]-bedHf$V3[i-1]>0){
    neededpos = seq(bedHf$V3[i-1], bedHf$V2[i]-1, 25)
    toappend = bind_rows(toappend, data.frame(V1=chr, V2=neededpos, V3=neededpos+25, V4=0))
  }
}
bedHf = bind_rows(bedHf, toappend) %>% arrange(V2)

range(lead(bedCf$V2)-bedCf$V3, na.rm = T)
range(lead(bedHf$V2)-bedHf$V3, na.rm = T)
head(bedCf)
head(bedHf)

# convert to long vectors (for each bp)
bedCv = rep(bedCf$V4, bedCf$V3-bedCf$V2)
bedHv = rep(bedHf$V4, bedHf$V3-bedHf$V2)
length(bedCv); length(bedHv)

# downsample to a manageable number
DSfactor = 500
bedCds = data.frame(pos = round((seq_along(bedCv)+actualstart)/DSfactor)*DSfactor, m=bedCv) %>%
  group_by(pos) %>%
  summarize(m = mean(m))
bedHds = data.frame(pos = round((seq_along(bedHv)+actualstart)/DSfactor)*DSfactor, m=bedHv) %>%
  group_by(pos) %>%
  summarize(m = mean(m))

ggplot(bind_rows(bedCds,bedHds, .id="trt")) + geom_point(aes(pos, m, col=trt), alpha=0.5)


## WRAPPERS

capa_wrapper = function(ts, theta, est_sd, pen, len, DSfactor){
    data_norm = (ts-theta)/est_sd
    res.anom = capa.uv(data_norm, beta=pen, beta_tilde=pen, type="mean",
                       min_seg_len=2, max_seg_len=len, transform=identity)
    # Extract segments, and convert effect sizes to raw means (similar to Alg2)
    res.anom.c = collective_anomalies(res.anom)[,1:3]
    res.anom.c$mean.change = apply(res.anom.c, 1, function(x) mean(ts[x[1]:x[2]]))
    # attach point anomalies, if any
    res.anom.p = point_anomalies(res.anom)
    if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
        res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location,
                                mean.change=ts[res.anom.p$location])
        res.anom.c = rbind(res.anom.c, res.anom.p)
    }
    # we are interested only in segments of increased mean:
    if(nrow(res.anom.c)==0){
        return(data.frame())
    }
    anomincr = res.anom.c[res.anom.c[,3]>theta,, drop=F]
    anomres = data.frame(chrom=chr,
                         chromStart=anomincr[,1]*DSfactor + actualstart,
                         chromEnd=anomincr[,2]*DSfactor + actualstart+DSfactor, 
                         segtype="S",
                         row.names = NULL)
}

not_wrapper = function(ts, theta, est_sd, pen, DSfactor){
    res.not = not(ts, method="not", contrast="pcwsConstMean")
    # res.not = features(res.not)$cpt
    res.not = features(res.not, method="threshold", th=pen)$cpt
    
    res.not = cbind(c(1, res.not),
                    c(res.not-1, length(ts)))
    res.not = cbind(res.not, apply(res.not, 1, function(x){
        if(any(is.na(x))) { NA }
        else { mean(ts[x[1]:x[2]]) }
    }))
    notincr = res.not[res.not[,3]>theta,, drop=F]
    notres = data.frame(chrom=chr,
                        chromStart=notincr[,1]*DSfactor + actualstart,
                        chromEnd=notincr[,2]*DSfactor + actualstart+DSfactor, 
                        segtype=ifelse(notincr[,3]>theta+est_sd, "S", "S2"),
                        row.names = NULL)
}

alg2_wrapper = function(ts, theta, est_sd, pen, len, DSfactor){
    alg2 = fulldetector_prune(ts, theta0=theta, MAXLOOKBACK=len,
                              PEN=pen, PEN2=pen, SD=est_sd, prune=2)

    # we are interested only in segments of increased mean:
    alg2incr = alg2$segs[alg2$segs[,3]>theta,, drop=F]
    
    alg2res = data.frame(chrom=chr,
                         chromStart=alg2incr[,1]*DSfactor + actualstart,
                         chromEnd=alg2incr[,2]*DSfactor + actualstart+DSfactor, 
                         segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                         row.names = NULL)
}


## ---- PENALTY TRAINING ----
# Shared training data:
ts = bedCds$m
tsdf = data.frame(chrom="chr1", chromStart=as.integer(bedCds$pos),
                  chromEnd=as.integer(bedCds$pos+DSfactor),
                  count=as.integer(round(bedCds$m)))
finaltheta = median(bedCds$m)
maxpeaklen = round(50e3/DSfactor)
pen = autoset_penalty(ts)
robust_sd = mad(diff(ts))

## Apply standard epidemic chp detection (anomaly)
for(penmult in c(1, 3, 5, 10, 20)){
    res = capa_wrapper(ts, finaltheta, robust_sd, penmult*pen, maxpeaklen, DSfactor)
    print(res)
    cat(sprintf("Penalty %.2f / Number of detections %d\n", penmult*pen, nrow(res)))
}
# Could use x5 or x10 for capa (3/2 dets)

## Narrowest-over-threshold detector
for(penmult in c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8)){
    res = not_wrapper(ts, finaltheta, robust_sd, penmult*pen, DSfactor)
    print(res)
    cat(sprintf("Penalty %.2f / Number of detections %d\n", penmult*pen, nrow(res)))
}
# Could use 0.2 (3 dets) or 0.3 (2 dets)

## PeakSegFPOP
for(lambda in c(5e2, 1e3, 3e3, 5e3, 7e3, 1e4)){
    res.fpop = PeakSegFPOP_df(tsdf, lambda, "/tmp/")$segments
    fpopincr = filter(data.frame(res.fpop), status=="peak", mean>finaltheta)
    print(fpopincr)
    cat(sprintf("Penalty %d / Number of detections %d\n", lambda, nrow(fpopincr)))
}
# 5e3 gives 3 dets

## Proposed Alg2
for(penmult in c(3, 5, 10)){
    sink("/tmp/alg2.log")
    res = alg2_wrapper(ts, finaltheta, robust_sd, penmult*pen, maxpeaklen, DSfactor)
    sink()
    print(res)
    cat(sprintf("Penalty %.2f / Number of detections %d\n", penmult*pen, nrow(res)))
}
# 1 gives 5, 3 and 5 give same 2, 10 gives 1. x5 matches CAPA so use that


## ---- DETECTION RUNS -----

# Shared detection parameters
ts = bedHds$m
finaltheta = median(bedHds$m)
maxpeaklen = round(50e3/DSfactor)
pen = autoset_penalty(ts)
tsdf = data.frame(chrom="chr1", chromStart=as.integer(bedHds$pos),
                  chromEnd=as.integer(bedHds$pos+DSfactor),
                  count=as.integer(round(bedHds$m)))
robust_sd = mad(diff(bedCds$m))


## Apply standard epidemic chp detection (anomaly)
data_norm = (ts-finaltheta)/robust_sd
res.anom = capa.uv(data_norm, beta=5*pen, beta_tilde=5*pen, type="mean",
                   min_seg_len=2, max_seg_len=maxpeaklen, transform=identity)
# Extract segments, and convert effect sizes to raw means (similar to Alg2)
res.anom.c = collective_anomalies(res.anom)[,1:3]
res.anom.c$mean.change = apply(res.anom.c, 1, function(x) mean(ts[x[1]:x[2]]))
# attach point anomalies, if any
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

## Apply NOT, Narrowest-over-threshold detector
res.not = not(ts, method="not", contrast="pcwsConstMean")
# res.not = features(res.not)$cpt
res.not = features(res.not, method="threshold", th=0.2*pen)$cpt

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

## Apply Algorithm 2
alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen,
                          PEN=5*pen, PEN2=5*pen, SD=robust_sd, prune=2)
# we are interested only in segments of increased mean:
alg2incr = alg2$segs[alg2$segs[,3]>finaltheta,, drop=F]
alg2res = data.frame(chrom=chr,
                     chromStart=alg2incr[,1]*DSfactor + actualstart,
                     chromEnd=alg2incr[,2]*DSfactor + actualstart+DSfactor, 
                     segtype=ifelse(alg2incr[,4]==1, "S", "N"),
                     row.names = NULL)

## Apply PeakSegFPOP
res.fpop = PeakSegFPOP_df(tsdf, 5e3, "/tmp/")$segments
fpopincr = filter(data.frame(res.fpop), status=="peak", mean>finaltheta)
fpopincr$segtype = "S"
print(fpopincr)


## Plot all results
toplotAn = mutate(anomres, treatment=2)
toplotNot = mutate(notres, treatment=2)
toplotFP = mutate(fpopincr, treatment=2)
toplot2 = mutate(alg2res, treatment=2)
toplotL = data.frame(label=c("anomaly", "not", "PeakSeg", "proposed"), x=120700, y=c(-3, -6, -9, -12), treatment=2)
p3 = bind_rows(bedCds,bedHds, .id="treatment") %>%
  ggplot() +
  # geom_point(aes(pos/1000, m), alpha=0.5) +
  geom_point(aes(pos/1000, m), shape=21, col="black", fill="grey50", stroke=0.5) +  # to allow nice exporting to eps
  geom_segment(aes(chromStart/1000, -3, xend=chromEnd/1000, yend=-3, color=segtype),
               data=toplotAn, size=6)+
  geom_segment(aes(chromStart/1000, -6, xend=chromEnd/1000, yend=-6, color=segtype),
               data=toplotNot, size=6)+
  geom_segment(aes(chromStart/1000, -9, xend=chromEnd/1000, yend=-9, color=segtype),
               data=toplotFP, size=6)+
  geom_segment(aes(chromStart/1000, -12, xend=chromEnd/1000, yend=-12, color=segtype, size=segtype),
               data=toplot2) +
  geom_point(aes(pos/1000, m), alpha=0, data=bedHds[which.max(bedHds$m),]) +  # invisible point to adjust scales nicer
  geom_text(aes(x=x, y=y, label=label), data=toplotL, hjust=0) + 
  facet_grid(treatment~., labeller = labeller(treatment=c("1"="control", "2"="H3K27ac")),
             scales = "free_y", space="free_y") + 
  xlab("position on chromosome 1, kbp") + ylab("coverage") +
  xlim(c(120100, 120780)) +
  scale_color_manual(values=c("S"="blue", "N"="gold", "S2"="#D8EFFF"),
                     name="detections", labels=c("nuisance", "signal", "sign. low")) +
  scale_size_manual(values=c("N"=7, "S"=6), guide="none")+
  theme_bw() + theme(legend.position = "right")

print(p3)

# note that no segments where removed due to the < median constraint.

# save output
outprefix = "../drafts/changepoint-method/results-appl/broad-chr1-120mbp-ds500-"
write.table(anomres, paste0(outprefix, "anom.tsv"), quote=F, row.names=F, sep="\t")
save(res.anom, file=paste0(outprefix, "anom.RData"))
write.table(notres, paste0(outprefix, "not.tsv"), quote=F, row.names=F, sep="\t")
save(res.not, file=paste0(outprefix, "not.RData"))
write.table(fpopincr, paste0(outprefix, "fpop.tsv"), quote=F, row.names=F, sep="\t")
save(res.fpop, file=paste0(outprefix, "fpop.RData"))
write.table(alg2res, paste0(outprefix, "alg2.tsv"), quote=F, row.names=F, sep="\t")
save(alg2, file=paste0(outprefix, "alg2.RData"))
# ggsave(paste0(outprefix, "merged.png"), plot=p3, width=18, height=15, units="cm")
ggsave(paste0(outprefix, "merged.eps"), plot=p3, width=18, height=15, units="cm", device=cairo_ps)


## ---- Calculate fit statistics (BIC/SIC) -----

# Calculates likelihood assuming a normal model with current finaltheta and SD.
# For nuisance segments, uses only the points N\S.
#  - df: a dataframe of only the segments detected.
#  - Requires "start" and "end" columns in ts positions, and "segtype" with S or N
getBIC = function(df, ts, mu0, sigma0){
  nparams = 0
  lastsegend = 0
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
res.fpop.segs$start = match(res.fpop.segs$chromStart, bedHds$pos)
res.fpop.segs$end = match(res.fpop.segs$chromEnd, bedHds$pos)
bic.fpop = getBIC(res.fpop.segs, ts, finaltheta, robust_sd)

res.alg2.segs = data.frame(alg2$segs)
colnames(res.alg2.segs) = c("start", "end", "mean", "segtype")
res.alg2.segs$segtype = ifelse(res.alg2.segs$segtype==1, "S", "N")
bic.alg2 = getBIC(res.alg2.segs, ts, finaltheta, robust_sd)

bic.anom
bic.not
bic.fpop
bic.alg2
