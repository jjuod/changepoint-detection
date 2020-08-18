# Real data analysis: Spanish mortality data

options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

library(anomaly)
library(not)
library(PeakSegDisk)
setwd("~/Documents/gitrep/changepoint-detection/")
source("MAIN-algorithm.R")

setwd("~/Documents/kiwis/changepointdet/covid/")

dt = read.csv("eurostat/demo_r_mweek3_1_Data.csv")

dt = mutate(dt, ValueNum = sub(":", "", Value)) %>%
  mutate(ValueNum = as.numeric(sub(",", "", ValueNum)))
dt

tsdf = filter(dt, GEO=="Spain", AGE=="Y60-64") %>%
  group_by(TIME) %>%
  summarize(ValueNum = mean(ValueNum))
  
tsdf = separate(tsdf, "TIME", c("YR", "WK"), sep="W", remove=F) %>%
  mutate(YR=as.numeric(YR), WK=as.numeric(WK))

# two missing values at the end:
tail(tsdf)
which(is.na(tsdf$ValueNum))
# also filter to an exact 4 yr period (2017W22 to 2020W21)
tsdf = filter(tsdf, !is.na(ValueNum), YR>2017 | WK>21)
table(tsdf$YR)


## DETECTION RUNS

## Shared parameters
ts = tsdf$ValueNum
plot(ts)
pen = autoset_penalty(ts)
maxpeaklen = 10
SD = sd(ts[1:52])
finaltheta = median(ts[1:52])

## Apply standard epidemic chp detection (anomaly)
res.anom = capa.uv((ts-finaltheta)/sd(ts), beta=pen, beta_tilde=pen, type="mean", min_seg_len=2, max_seg_len=maxpeaklen, transform=identity)

# Extract segments, and convert effect sizes to raw means (similar to Alg2)
res.anom.c = collective_anomalies(res.anom)[,1:3]
res.anom.c$mean.change = apply(res.anom.c, 1, function(x) mean(ts[x[1]:x[2]]))
res.anom.p = point_anomalies(res.anom)
if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
  res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location,
                          mean.change=ts[res.anom.p$location])
  res.anom.c = rbind(res.anom.c, res.anom.p)
}
colnames(res.anom.c)[3] = "theta"
print(res.anom.c)

## Apply NOT
# Narrowest-over-threshold detector
res.not = not(ts, method="not", contrast="pcwsConstMean")
res.not = features(res.not)$cpt
res.not = cbind(c(1, res.not),
                c(res.not-1, length(ts)))
res.not = cbind(res.not, apply(res.not, 1, function(x){
  if(any(is.na(x))) { NA }
  else { mean(ts[x[1]:x[2]]) }
}))
res.not = data.frame(res.not)
res.not = filter(res.not, abs(X3-finaltheta)>SD)
colnames(res.not) = c("start", "end", "theta")
print(res.not)

## Apply PeakSegFPOP
fpopdf = data.frame(chrom="chr1", chromStart=as.integer(seq_along(ts)), chromEnd=as.integer(seq_along(ts)+1),
                  count=as.integer(round(ts)))
for(lambda in c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6)){
  res.fpop = PeakSegFPOP_df(fpopdf, lambda, "/tmp/")$segments
  fpopincr = filter(data.frame(res.fpop), status=="peak", mean>finaltheta)
  cat(sprintf("Penalty %d / Number of detections %d\n", lambda, nrow(fpopincr)))
}
res.fpop = PeakSegFPOP_df(fpopdf, 1e2, "/tmp/")$segments
fpopincr = filter(data.frame(res.fpop), status=="peak")
fpopincr = mutate(fpopincr[,c("chromStart", "chromEnd", "mean")], segtype="S")
colnames(fpopincr) = c("start", "end", "theta", "segtype")
print(fpopincr)

## Apply Algorithm 2:
alg2 = fulldetector_prune(ts, theta0=finaltheta, MAXLOOKBACK=maxpeaklen, PEN=pen, PEN2=pen, SD=SD, prune=2)
print(alg2$segs)

alg2res = data.frame(alg2$segs)
colnames(alg2res) = c("start", "end", "theta", "segtype")
rownames(alg2res) = NULL
alg2res$segtype = ifelse(alg2res$segtype==1, "S", "N")

# fix time for plotting
tsdf$date = parse_date_time(paste(tsdf$YR, tsdf$WK, "1"), "%Y %W %w")

toplot = bind_rows(anomaly=mutate(res.anom.c, segtype="S"),
          PeakSeg=fpopincr,
          not=mutate(res.not, segtype="S"),
          proposed=alg2res,
          .id="method")

toplot$stdate = min(tsdf$date) + dweeks(toplot$start-1)
toplot$edate = min(tsdf$date) + dweeks(toplot$end-1)
toplot$y = ifelse(toplot$method=="anomaly", 260, ifelse(toplot$method=="PeakSeg", 200, ifelse(toplot$method=="not", 230, 165)))
toplotL = group_by(toplot, method) %>% top_n(1, start)

p2 = ggplot(tsdf) +
  geom_point(aes(x=as.Date(date), y=ValueNum), size=0.5) +
  geom_segment(aes(as.Date(stdate), y, xend=as.Date(edate), yend=y, color=segtype, size=segtype),
               data=toplot) +
  geom_text(aes(x=as.Date(max(tsdf$date)+dweeks(1)), y=y, label=method), hjust=0, data=toplotL) +
  scale_color_manual(values=c("S"="blue", "N"="gold"), name="detections", labels=c("nuisance", "signal")) +
  ylab("weekly deaths") + xlab(NULL) + # ggtitle("Weekly death count in Spain over 2017-2020, ages 60-64") +
  scale_x_date(date_breaks = "3 months", date_labels="%Y-%m", minor_breaks = NULL,
               limits = c(as.Date(min(tsdf$date)), as.Date(max(tsdf$date) + dweeks(16)))) +
  scale_size_manual(values=c("false positive"=2, correct=0.7, "N"=5.5, "S"=4.8)) +
  guides(size="none", color=guide_legend(order=2, override.aes=list(size=7))) +
  theme_bw()+ theme(text=element_text(size=12), axis.text.x=element_text(size=12, angle=30, vjust=0.8, hjust=1),
                    legend.position = c(0.5, 0.8), legend.box.background = element_rect(color="grey50"),
                    legend.box.margin = margin(4,8,4,8))
print(p2)

outprefix = "../drafts/changepoint-method/results-appl/det-mortality-"
save(res.anom, file=paste0(outprefix, "anom.RData"))
save(res.not, file=paste0(outprefix, "not.RData"))
save(res.fpop, file=paste0(outprefix, "fpop.RData"))
save(alg2, file=paste0(outprefix, "alg2.RData"))
ggsave(paste0("../drafts/changepoint-method/results-appl/det-mortality", ".png"), width=17, height=10, units="cm")


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

bic.anom = getBIC(res.anom.c, ts, finaltheta, SD)

# define "background" segments
res.not.segs = data.frame(res.not[res.not[,3]-finaltheta>0, ])
colnames(res.not.segs) = c("start", "end", "mean")
bic.not = getBIC(res.not.segs, ts, finaltheta, SD)

# extract background and convert back to ts coordinates
res.fpop.segs = filter(res.fpop, status=="peak") %>% arrange(chromStart)
res.fpop.segs$start = res.fpop.segs$chromStart
res.fpop.segs$end = res.fpop.segs$chromEnd
bic.fpop = getBIC(res.fpop.segs, ts, finaltheta, SD)

res.alg2.segs = data.frame(alg2$segs)
colnames(res.alg2.segs) = c("start", "end", "mean", "segtype")
res.alg2.segs$segtype = ifelse(res.alg2.segs$segtype==1, "S", "N")
bic.alg2 = getBIC(res.alg2.segs, ts, finaltheta, SD)

bic.anom
bic.not
bic.fpop
bic.alg2

