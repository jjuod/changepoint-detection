# Saddleback analysis script
# for ecoacoustics congress 2020 abstract

options(stringsAsFactors = F)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rjson)
library(cowplot)
library(RcppCNPy)

## 0. DIRECTORY SETUP / READ IN DATA
setwd("~/Documents/audiodata/Saddleback/Train/")

# export JSON in AviaNZ format
exportAviaNZ <- function(dets, outfile){
  specieslabel = list(filter="D", species="Saddleback (Nth Is)", certainty=100)
  annotexp = mutate(dets, fl=0, fh=0, sp=rep(list(list(specieslabel)), nrow(dets)))
  annotexp = t(annotexp)
  rownames(annotexp) = NULL
  annotexp = as.list(as.data.frame(annotexp))
  names(annotexp) = NULL
  annotexp = c(list(list(Operator="Auto", Reviewer="Auto", Duration=300.0)), annotexp)
  annotexp = toJSON(annotexp)
  annotexp
  writeChar(annotexp, paste0("~/Documents/audiodata/Saddleback/Train/", outfile, ".data"),
            eos=NULL)
}

# import JSON from AviaNZ format
importAviaNZ <- function(infile){
  a = fromJSON(file=paste0("../VLannots/", infile, ".data"))
  
  if(length(a)<=1){
    next
  }
  
  a = a[-1] # drop metadata
  a = data.frame(t(sapply(a, c))) # to dataframe
  # a$time = parse_date_time(substr(wavfile, 4, 18), "Ymd_HMS") # time from filename
  # a$start = a$time + seconds(a$X1) # call start time to absolute time
  # a$end = a$time + seconds(a$X2)
  colnames(a) = c("start", "end", "flow", "fhigh", "species")
  
  # just convert these from lists
  a = unnest(a, start, end, flow, fhigh)
  
  # deal with new format nested species labels
  a = unnest(a, species)
  
  # separate out species names (drop filter name, certainty)
  a = unnest(a, species)
  a = a[seq(2, nrow(a), 3),]
  
  # remove 1-5 labels which were used to mark quality here
  a = unnest(a, species)
  a = filter(a, !species %in% 1:5)
  
  a$filename = infile
  return(a)
}

# parse annotations
wavfiles = list.files(pattern="c7.*.wav$")

annot = data.frame()
for(wavfile in wavfiles){
  a = importAviaNZ(wavfile)
  annot = rbind(annot, a)
}

# check results (all annotations):
nrow(annot)
table(annot$species)

# tieke positive times
tieke = filter(annot, species=="Saddleback (Nth Is)")

# produce 0/1 ground truth for file 1
# 500 wcs per second, 301 second for rounding errors
times = rep(F, 301*500)
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_130504.wav"){
    start = round(tieke[s,"start"]*500)
    end = round(tieke[s,"end"]*500)
    times[start:end] = T
  }
}

# load the wavelet coefs for all nodes, all timepoints
energy1 = npyLoad(paste0("c7_20181224_130504.wav", ".energies.npy"))
energy1 = t(energy1)

# we should actually downsample b/c the leaves are stored non-downsampled
energy1 = energy1[seq(1, nrow(energy1), 2),]
dim(energy1)

colnames(energy1) = paste0("n", 31:62)


## 1. standard "wavelet filter", no smoothing
## (max log E over 10 WCs)
logE1 = log(abs(energy1[,"n40"]))
logE2 = log(abs(energy1[,"n46"]))

# Convert pres/abs into intervals
getintervals = function(ts){
  dfint = data.frame(t=which(ts))
  dfint = mutate(dfint, contbefo=t-lag(t)==1, contpast=t-lead(t)==-1)
  # fix ends
  dfint$contbefo[is.na(dfint$contbefo)] = F
  dfint$contpast[is.na(dfint$contpast)] = F
  # things that have no continuity before and have continuity past are starts
  # or those that have no continuity either way
  starts = dfint$t[!dfint$contbefo]
  ends = dfint$t[!dfint$contpast]
  
  # convert time to true seconds
  samplerate = round(length(ts)/300)
  starts = starts/samplerate - 1/samplerate
  ends = ends/samplerate
  dfint = data.frame(starts, ends)
  
  return(dfint)
}

bind_cols(e=logE2, t=times[1:nrow(energy1)]) %>%
  ggplot() + geom_density(aes(fill=t, x=e), alpha=0.5)

df = data.frame(times=1:nrow(energy1), gt = times[1:nrow(energy1)], e1=logE1, e2=logE2)
df = mutate(df, pred1 = logE1 > mean(logE1)+sd(logE1)*1.5,
            pred2 = logE2 > mean(logE2)+sd(logE2)*1.5,
            pred12 = pred1 | pred2)
table(df$pred1)
table(df$pred2)
table(df$pred12)

dfw = group_by(df, t = round(times/10)) %>%
  summarize(pred12 = any(pred12), gt = any(gt))

dfw %>% filter(t>170*50, t<200*50) %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred12), pch="|", size=8) + theme_bw()

dfw %>% filter(t>185*50, t<200*50) %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred12), pch="|", size=8) + theme_bw()

# total number of intervals found
dets1 = getintervals(dfw$pred12)
dets1
nrow(dets1); sum(dets1$ends-dets1$starts); table(dfw$gt, dfw$pred12, dnn=c("GT", "PRED"))


## 2. standard "wavelet filter", smoothing for 0.5 s
## (smooth = windowed mean on the original scale)
df = data.frame(times=1:nrow(energy1), gt = times[1:nrow(energy1)],
                e1=abs(energy1[,"n40"]), e2=abs(energy1[,"n46"]))

dfw = group_by(df, t = round(times/250)) %>%
  summarize(e1 = mean(e1), e2 = mean(e2), gt = any(gt))

dfw = mutate(dfw, pred1 = log(e1) > mean(log(e1))+sd(log(e1))*0.4,
            pred2 = log(e2) > mean(log(e2))+sd(log(e2))*0.4,
            pred12 = pred1 | pred2)
table(dfw$pred1)
table(dfw$pred2)
table(dfw$pred12)

p1 = dfw %>%
  ggplot() + geom_point(aes(x=t/2, y=gt, color=pred12), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="wavelet detector output") + theme_bw() +
  theme(legend.position = "top")
p1

# total number of intervals found
dets2 = getintervals(dfw$pred12)
dets2
nrow(dets2); sum(dets2$ends-dets2$starts); table(dfw$gt, dfw$pred12, dnn=c("GT", "PRED"))

# the two bands of interest
# ggplot(dfw) + geom_line(aes(x=t, y=log(e1)-mean(log(e1)))) +
  # geom_line(aes(x=t, y=log(e2)-mean(log(e2))), color="darkolivegreen") + theme_bw()


## 3. Cross correlation ??????
# load the audio xcorrs for all timepoints
xcorr1 = npyLoad(paste0("c7_20181224_130504.wav", ".xcorr.npy"))

# downsample to 1 measure per 10 ms
xcorr1 = data.frame(times = seq_along(xcorr1), r = xcorr1)
xcorr1 = group_by(xcorr1, t=round(times/160)) %>%
  summarize(r = max(abs(r)))
ggplot(xcorr1) + geom_line(aes(x=t, y=r))


## 4. Changepoint detection, assuming known theta0
## (standard OP algorithm)

cost0 <- function(xs, theta0){
  liks = dnorm(xs, 0, theta0)
  -2 * sum(log(liks))
}
cost1 <- function(xs){
  liks = dnorm(xs, 0, sd(xs))
  -2 * sum(log(liks))
}

detectchp <- function(ts, MAXLOOKBACK, PEN, THETA0){
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0, 0)
  # for each t, a matrix of segment starts-ends
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  segs[[2]] = matrix(0, nrow=1, ncol=3)
  
  for(t in 3:length(ts)){
    # if this is from background:
    bgcost = bestcost[t-1] + cost0(ts[t], THETA0)
    
    # if this is from segment, min over k:
    # (seg length >= 2)
    # (maxlookback limited if t short)
    segcost = rep(Inf, min(MAXLOOKBACK, t-2))
    for(k in seq_along(segcost)){
      segcost[k] = bestcost[t-k-1] + cost1(ts[(t-k):t])
    }
    
    # first x of the proposed segment:
    bestsegstart = which.min(segcost)  # in k = negative coord relative to t
    bestsegcost = segcost[bestsegstart] + PEN
    bestsegstart = t - bestsegstart  # in absolute position
    # print(sprintf("Best segment was %d-%d: %.1f = %.1f + F(t-k-1)",
    # bestsegstart, t, segcost[t-bestsegstart], cost1(ts[bestsegstart:t])))
    
    # is background better than seg?
    if(bgcost < bestsegcost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
      # print(sprintf("It was worse than C0 at %.1f", bgcost))
    } else {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, sd(ts[bestsegstart:t]))
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
    }
  }
  
  segs[[length(segs)]]
}

# check at most this many ks:
MAXLOOKBACK = 30
# penalty for starting a segment
PEN = 30

# assumed background sd (not var!)
# THETA0 = 1

# test time series
# ts = c(rnorm(100, 0, 1), rnorm(50, 0, 3), rnorm(50, 0, 1))
# detectchp(ts, MAXLOOKBACK, PEN, THETA0)

# REAL data, downsampled to 10/s
tsE1 = energy1[seq(1, nrow(energy1), 50),"n40"]
tsE2 = energy1[seq(1, nrow(energy1), 50),"n46"]

segs1 = detectchp(tsE1, MAXLOOKBACK, PEN, THETA0=sd(tsE1[1:100]))
segs2 = detectchp(tsE2, MAXLOOKBACK, PEN, THETA0=sd(tsE2[1:100]))
segs1
segs2

dets4a = data.frame(starts=segs1[-1,1], ends=segs1[-1,2])
dets4b = data.frame(starts=segs2[-1,1], ends=segs2[-1,2])
# convert to actual seconds
samplerate = round(length(tsE1)/300)
dets4a = mutate(dets4a, starts = starts/samplerate, ends = ends/samplerate)
dets4b = mutate(dets4b, starts = starts/samplerate, ends = ends/samplerate)

# convert segment time boundaries to pres/abs
df = data.frame(t = seq_along(tsE1)/samplerate, gt=F, pred1=F, pred2=F)
for(s in 1:nrow(segs1)){
  df$pred1[segs1[s,1]:segs1[s,2]] = T
}
for(s in 1:nrow(segs2)){
  df$pred2[segs2[s,1]:segs2[s,2]] = T
}
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_130504.wav"){
    start = round(tieke[s,"start"]*samplerate)
    end = round(tieke[s,"end"]*samplerate)
    df$gt[start:end] = T
  }
}

df = mutate(df, pred12 = pred1 | pred2)

p2 = df %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred12), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="changepoint det. output") + theme_bw() +
  theme(legend.position = "top")
p2

# export JSON in AviaNZ format
exportAviaNZ(dets4a, "c7_20181224_130504.wav")


## 5. Changepoint detection, OUR ALGORITHM (unknown theta0)
## (one iteration of stochastic fitting, and one batch. No formal converging yet)

detectchpU <- function(ts, MAXLOOKBACK, PEN){
  # bestcost[t] := F(all x[1:t])
  theta0 = rep(0, length(ts))
  bestcost = c(0, 0)
  # for each t, a matrix of segment starts-ends
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  segs[[2]] = matrix(0, nrow=1, ncol=3)
  
  # initialize theta0
  theta0[2] = sd(ts[1:2])
  for(t in 3:length(ts)){
    # if this is from background:
    bgcost = bestcost[t-1] + cost0(ts[t], theta0[t-1])
    
    # if this is from segment, min over k:
    # (seg length >= 2)
    # (maxlookback limited if t short)
    segcost = rep(Inf, min(MAXLOOKBACK, t-2))
    for(k in seq_along(segcost)){
      segcost[k] = bestcost[t-k-1] + cost1(ts[(t-k):t])
    }
    
    # first x of the proposed segment:
    bestsegstart = which.min(segcost)  # in k = negative coord relative to t
    bestsegcost = segcost[bestsegstart] + PEN
    bestsegstart = t - bestsegstart  # in absolute position
    # print(sprintf("Best segment was %d-%d: %.1f = %.1f + F(t-k-1)
    #               vs cost0: %.1f | theta0: %.2f",
    #               bestsegstart, t, segcost[t-bestsegstart], cost1(ts[bestsegstart:t]),
    #               bgcost, theta0[t-1]))
    # 
    # is background better than seg?
    if(bgcost < bestsegcost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
      # find current background set
      bgpoints = 1:t
      if(nrow(segs[[t]])>1){
        segpoints = c()
        for(s in 2:nrow(segs[[t]])){
          segpoints = c(segpoints, (segs[[t]][s,1]) : (segs[[t]][s,2]))
        }
        bgpoints = bgpoints[-segpoints]
      }
      # update theta0
      theta0[t] = sd(ts[bgpoints])
    } else {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, sd(ts[bestsegstart:t]))
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
      theta0[t] = theta0[bestsegstart-1]
    }
  }
  # one final iteration with the right theta0
  print(sprintf("Estimated theta0: %.2f", theta0[t]))
  res = detectchp(ts, MAXLOOKBACK, PEN, theta0[t])
  res
}

# check at most this many ks:
MAXLOOKBACK = 30
# penalty for starting a segment
PEN = 30

# test time series
#ts = c(rnorm(100, 0, 1), rnorm(50, 0, 3), rnorm(50, 0, 1))
#detectchpU(ts, MAXLOOKBACK, PEN)

# REAL data, downsampled to 10/s
tsE1 = energy1[seq(1, nrow(energy1), 50),"n45"]
tsE2 = energy1[seq(1, nrow(energy1), 50),"n46"]

segs1 = detectchpU(tsE1, MAXLOOKBACK, PEN)
segs2 = detectchpU(tsE2, MAXLOOKBACK, PEN)
segs1
segs2

dets5a = data.frame(starts=segs1[-1,1], ends=segs1[-1,2])
dets5b = data.frame(starts=segs2[-1,1], ends=segs2[-1,2])
# convert to actual seconds
samplerate = round(length(tsE1)/300)
dets5a = mutate(dets5a, starts = starts/samplerate, ends = ends/samplerate)
dets5b = mutate(dets5b, starts = starts/samplerate, ends = ends/samplerate)


# convert segment time boundaries to pres/abs
df = data.frame(t = seq_along(tsE1)/5, gt=F, pred1=F, pred2=F)
for(s in 1:nrow(segs1)){
  df$pred1[segs1[s,1]:segs1[s,2]] = T
}
for(s in 1:nrow(segs2)){
  df$pred2[segs2[s,1]:segs2[s,2]] = T
}
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_130504.wav"){
    start = round(tieke[s,"start"]*samplerate)
    end = round(tieke[s,"end"]*samplerate)
    df$gt[start:end] = T
  }
}

df = mutate(df, pred12 = pred1 | pred2)

p3 = df %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred12), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="new algorithm output") + theme_bw() +
  theme(legend.position = "top")
p3

plot_grid(p1+xlab(NULL), p2+xlab(NULL), p3, nrow=3)

# export JSON in AviaNZ format
exportAviaNZ(dets5a, "c7_20181224_130504.wav")



###### TEST ON FILE 2

# produce 0/1 ground truth for file 2
# 500 wcs per second, 301 second for rounding errors
times = rep(F, 301*500)
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_134004.wav"){
    start = round(tieke[s,"start"]*500)
    end = round(tieke[s,"end"]*500)
    times[start:end] = T
  }
}

# load the wavelet coefs for all nodes, all timepoints
energy1 = npyLoad(paste0("c7_20181224_134004.wav", ".energies.npy"))
energy1 = t(energy1)

# we should actually downsample b/c the leaves are stored non-downsampled
energy1 = energy1[seq(1, nrow(energy1), 2),]
dim(energy1)

colnames(energy1) = paste0("n", 31:62)


## 1. standard "wavelet filter", no smoothing
## (max log E over 10 WCs)
logE1 = log(abs(energy1[,"n40"]))
logE2 = log(abs(energy1[,"n46"]))

bind_cols(e=logE2, t=times[1:nrow(energy1)]) %>%
  ggplot() + geom_density(aes(fill=t, x=e), alpha=0.5)

df = data.frame(times=1:nrow(energy1), gt = times[1:nrow(energy1)], e1=logE1, e2=logE2)
df = mutate(df, pred1 = logE1 > mean(logE1)+sd(logE1)*1.5,
            pred2 = logE2 > mean(logE2)+sd(logE2)*1.5,
            pred12 = pred1 | pred2)
table(df$pred1)
table(df$pred2)
table(df$pred12)

dfw = group_by(df, t = round(times/10)) %>%
  summarize(pred12 = any(pred12), gt = any(gt))

# total number of intervals found
dets1 = getintervals(dfw$pred12)
dets1
nrow(dets1); sum(dets1$ends-dets1$starts); table(dfw$gt, dfw$pred12, dnn=c("GT", "PRED"))

## 2. standard "wavelet filter", smoothing for 0.5 s
## (smooth = windowed mean on the original scale)
df = data.frame(times=1:nrow(energy1), gt = times[1:nrow(energy1)],
                e1=abs(energy1[,"n40"]), e2=abs(energy1[,"n46"]))

dfw = group_by(df, t = round(times/250)) %>%
  summarize(e1 = mean(e1), e2 = mean(e2), gt = any(gt))

dfw = mutate(dfw, pred1 = log(e1) > mean(log(e1))+sd(log(e1))*0.4,
             pred2 = log(e2) > mean(log(e2))+sd(log(e2))*0.4,
             pred12 = pred1 | pred2)
table(dfw$pred1)
table(dfw$pred2)
table(dfw$pred12)

p1 = dfw %>%
  ggplot() + geom_point(aes(x=t/2, y=gt, color=pred12), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="wavelet detector output") + theme_bw() +
  theme(legend.position = "top")
p1

# total number of intervals found
dets2 = getintervals(dfw$pred12)
dets2
nrow(dets2); sum(dets2$ends-dets2$starts); table(dfw$gt, dfw$pred12, dnn=c("GT", "PRED"))

# the two bands of interest
# ggplot(dfw) + geom_line(aes(x=t, y=log(e1)-mean(log(e1)))) +
# geom_line(aes(x=t, y=log(e2)-mean(log(e2))), color="darkolivegreen") + theme_bw()


## 3. Cross correlation ??????
# load the audio xcorrs for all timepoints
xcorr1 = npyLoad(paste0("c7_20181224_134004.wav", ".xcorr.npy"))

# downsample to 1 measure per 10 ms
xcorr1 = data.frame(times = seq_along(xcorr1), r = xcorr1)
xcorr1 = group_by(xcorr1, t=round(times/160)) %>%
  summarize(r = max(abs(r)))
xcorr1 = xcorr1[125:30125,]

samplerate = round(nrow(xcorr1)/300)
xcorr1$gt = F
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_134004.wav"){
    start = round(tieke[s,"start"]*samplerate)
    end = round(tieke[s,"end"]*samplerate)
    xcorr1$gt[start:end] = T
  }
}
xcorr1 = mutate(xcorr1, pred = r > 1.01*median(r))

gt3 = getintervals(xcorr1$gt)
ggplot(xcorr1, aes(x=t/samplerate)) + geom_line(aes(y=r)) +
  geom_rect(data=gt3, ymin=1e9, ymax=1.6e9, aes(xmin=starts, xmax=ends), fill="blue2", alpha=0.3) +
  geom_point(aes(y=1e9, col=pred), pch="|", size=8) + xlab("t, seconds") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="xcorr output") + theme_bw()

# total number of intervals found
dets3 = getintervals(xcorr1$pred)
dets3
nrow(dets3); sum(dets3$ends - dets3$starts); table(xcorr1$gt, xcorr1$pred, dnn=c("GT", "PRED"))

pxc = xcorr1 %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="cross-corr. output") + theme_bw() +
  theme(legend.position = "top")
pxc

## 4. Changepoint detection, assuming known theta0
## (standard OP algorithm)

# REAL data, downsampled to 10/s
tsE1 = energy1[seq(1, nrow(energy1), 50),"n40"]
tsE2 = energy1[seq(1, nrow(energy1), 50),"n46"]

segs1 = detectchp(tsE1, MAXLOOKBACK, PEN, THETA0=sd(tsE1[1:100]))
segs2 = detectchp(tsE2, MAXLOOKBACK, PEN, THETA0=sd(tsE2[1:100]))
segs1
segs2

dets4a = data.frame(starts=segs1[-1,1], ends=segs1[-1,2])
dets4b = data.frame(starts=segs2[-1,1], ends=segs2[-1,2])
# convert to actual seconds
samplerate = round(length(tsE1)/300)
dets4a = mutate(dets4a, starts = starts/samplerate, ends = ends/samplerate)
dets4b = mutate(dets4b, starts = starts/samplerate, ends = ends/samplerate)

# convert segment time boundaries to pres/abs
df = data.frame(t = seq_along(tsE1)/samplerate, gt=F, pred1=F, pred2=F)
for(s in 1:nrow(segs1)){
  df$pred1[segs1[s,1]:segs1[s,2]] = T
}
for(s in 1:nrow(segs2)){
  df$pred2[segs2[s,1]:segs2[s,2]] = T
}
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_134004.wav"){
    start = round(tieke[s,"start"]*samplerate)
    end = round(tieke[s,"end"]*samplerate)
    df$gt[start:end] = T
  }
}

df = mutate(df, pred12 = pred1 | pred2)

p2 = df %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred12), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="changepoint det. output") + theme_bw() +
  theme(legend.position = "top")
p2

# export JSON in AviaNZ format
exportAviaNZ(dets4a, "c7_20181224_134004.wav")


## 5. Changepoint detection, OUR ALGORITHM (unknown theta0)
## (one iteration of stochastic fitting, and one batch. No formal converging yet)

# REAL data, downsampled to 10/s
tsE1 = energy1[seq(1, nrow(energy1), 50),"n40"]
tsE2 = energy1[seq(1, nrow(energy1), 50),"n46"]

segs1 = detectchpU(tsE1, MAXLOOKBACK, PEN)
segs2 = detectchpU(tsE2, MAXLOOKBACK, PEN)
segs1
segs2

dets5a = data.frame(starts=segs1[-1,1], ends=segs1[-1,2])
dets5b = data.frame(starts=segs2[-1,1], ends=segs2[-1,2])
# convert to actual seconds
samplerate = round(length(tsE1)/300)
dets5a = mutate(dets5a, starts = starts/samplerate, ends = ends/samplerate)
dets5b = mutate(dets5b, starts = starts/samplerate, ends = ends/samplerate)


# convert segment time boundaries to pres/abs
df = data.frame(t = seq_along(tsE1)/5, gt=F, pred1=F, pred2=F)
for(s in 1:nrow(segs1)){
  df$pred1[segs1[s,1]:segs1[s,2]] = T
}
for(s in 1:nrow(segs2)){
  df$pred2[segs2[s,1]:segs2[s,2]] = T
}
for(s in 1:nrow(tieke)){
  if(tieke$filename[s]=="c7_20181224_134004.wav"){
    start = round(tieke[s,"start"]*samplerate)
    end = round(tieke[s,"end"]*samplerate)
    df$gt[start:end] = T
  }
}

df = mutate(df, pred12 = pred1 | pred2)

p3 = df %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pred12), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  scale_color_manual(values=c("chocolate1", "blue3"), name="new algorithm output") + theme_bw() +
  theme(legend.position = "top")
p3

plot_grid(p1+xlab(NULL), pxc+xlab(NULL), p3, nrow=3)

# export JSON in AviaNZ format
exportAviaNZ(dets5a, "c7_20181224_134004.wav")


############## ABSTRACT V2 ###############


# read in median clipping and Harma results
det = data.frame()
for (wavfile in wavfiles){
  det1 = read.table(paste0("~/Documents/gitrep/paper2/results-seg/", wavfile, ".segsh"), header=F)
  det1 = mutate(det1, start=gsub("[\\[,]", "", V1), end=gsub(",", "", V2),
                thr=gsub(",", "", V3), par2=gsub("\\]", "", V4), file=wavfile)
  det1 = det1[,-c(1:4)]
  
  det2 = read.table(paste0("~/Documents/gitrep/paper2/results-seg/", wavfile, ".segsm"), header=F)
  det2 = mutate(det2, start=gsub("[\\[,]", "", V1), end=gsub(",", "", V2),
                thr=gsub(",", "", V3), par2=gsub("\\]", "", V4), file=wavfile)
  det2 = det2[,-c(1:4)]

  det3 = read.table(paste0("~/Documents/gitrep/paper2/results-seg/", wavfile, ".segslm"), header=F)
  det3 = mutate(det3, start=gsub("[\\[,]", "", V1), end=gsub(",", "", V2),
                thr=gsub(",", "", V3), par2=gsub("\\]", "", V4), file=wavfile)
  det3 = det3[,-c(1:4)]
  
  det1 = bind_rows(harma=det1, med=det2, medlog=det3, .id="method")
  det = bind_rows(det, det1)
}
nrow(det)

detgroups = group_by(det, file, method, thr, par2) %>%
  summarise()
detgroups

detpresabs = data.frame()
for(gr in 1:nrow(detgroups)){
  # create an "Absent" df for one combo of file x params
  detonegr = semi_join(det, detgroups[gr,], by=c("method", "file", "thr", "par2"))
  detpresabs1 = detonegr[1, c("method", "file", "thr", "par2")]
  detpresabs1 = bind_rows(replicate(300, detpresabs1, simplify = F))
  detpresabs1$gt = F
  detpresabs1$pres = F
  
  # fill it out
  for(s in 1:nrow(detonegr)){
    start = floor(as.numeric(detonegr[s,"start"]))
    end = ceiling(as.numeric(detonegr[s,"end"]))
    detpresabs1$pres[start:end] = T
  }
  
  # add ground truth
  for(s in 1:nrow(tieke)){
    if(tieke$filename[s]==detonegr$file[1]){
      start = floor(tieke[s,"start"])
      end = ceiling(tieke[s,"end"])
      detpresabs1$gt[start:end] = T
    }
  }
  
  detpresabs = bind_rows(detpresabs, detpresabs1)
}

nrow(detpresabs)/300
table(detpresabs$gt)
table(detpresabs$pres)

# obtain sens, spec etc
detstats = group_by(detpresabs, method, thr, par2) %>%
  summarize(n=n(), TP=mean(gt & pres), TN=mean(!gt & !pres),
            FP = mean(!gt & pres), FN = mean(gt & !pres)) %>%
  mutate(sens = TP/(TP+FN), spec=TN/(TN+FP), ppv=TP/(TP+FP), npv=TN/(TN+FN))

## OPTIMIZE PARAMETERS:
# all methods
ggplot(detstats) + geom_point(aes(x=ppv, y=sens, col=par2)) +
  xlim(c(0,1)) + ylim(c(0,1)) + theme_bw() + facet_wrap(~method)

# pick max sens given >25 % prec

# Harma
filter(detstats, method=="harma") %>%
  ggplot() + geom_point(aes(x=ppv, y=sens, col=par2)) +
  xlim(c(0,1)) + ylim(c(0,1)) + theme_bw()
# best option: 79 % sens, 26 % ppv @ thr=10, par2=0.7
filter(detstats, method=="harma", ppv>0.25) %>% arrange(desc(sens))

# Median clipping
filter(detstats, method=="med") %>%
  ggplot() + geom_point(aes(x=ppv, y=sens, col=par2)) +
  xlim(c(0,1)) + ylim(c(0,1)) + theme_bw()
# best options: 91 % sens, 31 % ppv @ thr=4.0, par2=60
filter(detstats, method=="med", ppv>0.25) %>% arrange(desc(sens))

# Median clipping (log scale)
filter(detstats, method=="medlog") %>%
  ggplot() + geom_point(aes(x=ppv, y=sens, col=par2)) +
  xlim(c(0,1)) + ylim(c(0,1)) + theme_bw() 
# best options: 88 % sens, 28 % ppv @ thr=1.3, par2=70
filter(detstats, method=="medlog", ppv>0.25) %>% arrange(desc(sens))


## Plot results with optimal parameters
detgood = filter(detpresabs, method=="harma" & par2==0.7 & thr==10 | 
                   method=="med" & par2==60 & thr=="4.0" |
                   method=="medlog" & par2==70 & thr==1.3)

group_by(detgood, method) %>% mutate(t = row_number()) %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pres), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  facet_grid(method~.) +
  scale_x_continuous(breaks = seq(0, 1200, 300), minor_breaks = seq(0, 1200, 60)) +
  scale_color_manual(values=c("chocolate1", "blue3"), name="changepoint det. output") +
  theme_bw() + theme(legend.position = "top")


## Changepoint detection, assuming known theta0
## (standard OP algorithm)

cost0 <- function(xs, theta0){
  liks = dnorm(xs, 0, theta0)
  -2 * sum(log(liks))
}
cost1 <- function(xs){
  liks = dnorm(xs, 0, sd(xs))
  -2 * sum(log(liks))
}

detectchp <- function(ts, MAXLOOKBACK, PEN, THETA0){
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0, 0)
  # for each t, a matrix of segment starts-ends
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  segs[[2]] = matrix(0, nrow=1, ncol=3)
  
  for(t in 3:length(ts)){
    # get cost if this is from background:
    bgcost = bestcost[t-1] + cost0(ts[t], THETA0)
    
    # get costs if this is from segment of length k:
    # (seg length >= 2)
    # (maxlookback limited if t short)
    segcost = rep(Inf, min(MAXLOOKBACK, t-2))
    for(k in seq_along(segcost)){
      segcost[k] = bestcost[t-k-1] + cost1(ts[(t-k):t])
    }
    
    # first x of the proposed segment:
    bestsegstart = which.min(segcost)  # in k = negative coord relative to t
    bestsegcost = segcost[bestsegstart] + PEN
    bestsegstart = t - bestsegstart  # in absolute position
    # print(sprintf("Best segment was %d-%d: %.1f = %.1f + F(t-k-1)",
    # bestsegstart, t, segcost[t-bestsegstart], cost1(ts[bestsegstart:t])))
    
    # is background better than seg?f
    if(bgcost < bestsegcost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
      # print(sprintf("It was worse than C0 at %.1f", bgcost))
    } else {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, sd(ts[bestsegstart:t]))
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
    }
  }
  
  list(segs=segs[[length(segs)]], cost=bestcost[[length(bestcost)]])
}

analyzeFile <- function(filename, MAXLOOKBACK, PEN, plotme=F){
  # load the spec for all timepoints
  energy1 = npyLoad(paste0(filename, ".spec.npy"))
  
  # transform
  energy1 = apply(log(energy1), 1, max)
  
  # downsample to 10/s
  energy1DS = tapply(energy1, rep(1:3000, each=25, length.out=length(energy1)), mean)
  minmean = mean(c(median(energy1DS[1:100]), median(energy1DS[601:700]), median(energy1DS[1201:1300])))
  
  energy1DS = energy1DS - minmean
  dim(energy1DS)
  
  # Q: how to combine the multiple bands into one?
  # Let's say we average over freqs:
  tsE1 = energy1DS * rep(c(-1, 1), length.out=length(energy1DS))
  
  theta1 = min( c(sd(tsE1[1:50]), sd(tsE1[601:651]),
                  sd(tsE1[1201:1251]), sd(tsE1[1501:1551]),
                  sd(tsE1[1801:1851]), sd(tsE1[2101:2151]),
                  sd(tsE1[2401:2451]), sd(tsE1[2701:2751])) )

  # theta1 = sd(tsE1[1:50])
  message("using sd: ", theta1)
  res1 = detectchp(tsE1, MAXLOOKBACK, PEN, THETA0=theta1)
  
  # theta1 = sd(tsE1[601:651])
  # message("using sd: ", theta1)
  # res2 = detectchp(tsE1, MAXLOOKBACK, PEN, THETA0=theta1)
  # 
  # theta1 = sd(tsE1[1201:1251])
  # message("using sd: ", theta1)
  # res3 = detectchp(tsE1, MAXLOOKBACK, PEN, THETA0=theta1)

  # pick the objectively best theta
  # bestrun = which.min(c(res1$cost, res2$cost, res3$cost))
  # segs = list(res1$segs, res2$segs, res3$segs)[[bestrun]]
  segs = res1$segs
  if(plotme){
    plot(tsE1)
    for(s in 2:nrow(segs)){
      abline(v=segs[s,1], col="blue")
      abline(v=segs[s,2], col="blue")
      points(x=mean(segs[s,1:2]), y=6, col="blue", pch=16)
      text(x=mean(segs[s,1:2]), y=s %% 3+1, labels=round(segs[s,3], 2))
    }
    for(s in 1:nrow(tieke)){
      if(tieke$filename[s]==filename){
        abline(v=tieke[s,1]*10, col="green")
        abline(v=tieke[s,2]*10, col="green")
      }
    }
  }
  
  segs = data.frame(starts=segs[-1,1], ends=segs[-1,2])
  
  # convert to actual seconds
  samplerate = round(length(tsE1)/300)
  segs = mutate(segs, starts = starts/samplerate, ends = ends/samplerate)
  
  return(segs)
}

# check at most this many ks:
MAXLOOKBACK = 50
# penalty for starting a segment
PEN = 2100

# Prepare output
detnew = data.frame(method=rep("chp", 1200), file=rep(wavfiles, each=300),
                    thr="1", par2="1", gt=F, pres=F)
detnewsegs = data.frame()
# add ground truth 
for(s in 1:nrow(tieke)){
  filestart = which(detnew$file==tieke$filename[s])[1]
  start = floor(tieke[s,"start"]) + filestart
  end = ceiling(tieke[s,"end"]) + filestart
  detnew$gt[start:end] = T
}

# RUN DETECTION
for(wavfile in wavfiles){
  message("analyzing ", wavfile)
  dets4a = analyzeFile(wavfile, MAXLOOKBACK, PEN)
  detnewsegs = bind_rows(detnewsegs,
                data.frame(method="chp", start=dets4a$starts, end=dets4a$ends,
                thr=as.character(MAXLOOKBACK), par2=as.character(PEN),
                file = wavfile))
  message("obtained ", nrow(dets4a), " segments")
  filestart = which(detnew$file==wavfile)[1]
  detnew$pres[filestart:(filestart+299)] = F # just to quickly reset while testing
  
  # convert segment time boundaries to pres/abs
  for(s in 1:nrow(dets4a)){
    start = min(floor(as.numeric(dets4a[s,"starts"])), 299) + filestart
    end = min(ceiling(as.numeric(dets4a[s,"ends"])), 299) + filestart
    detnew$pres[start:end] = T
  }
}
# best combo: LB 50, PEN 800 = 81 % sens, 30 % prec, 14-11-10-7 segs (for Mean)
# best combo: LB 50, PEN 2100 = 87 % sens, 28 % prec, 15-12-12-11 segs (for Max)

bind_rows(detnew, detgood) %>%
  group_by(method) %>% mutate(t = row_number()) %>%
  ggplot() + geom_point(aes(x=t, y=gt, color=pres), pch="|", size=8) +
  ylab("true presence") + xlab("time, s") +
  facet_grid(method~.) +
  scale_x_continuous(breaks = seq(0, 1200, 300), minor_breaks = seq(0, 1200, 60)) +
  scale_color_manual(values=c("chocolate1", "blue3"), name="changepoint det. output") +
  theme_bw() + theme(legend.position = "top")

# export JSON in AviaNZ format
# exportAviaNZ(dets4a, "c7_20181224_130504.wav")

# obtain sens, spec etc
bind_rows(detnew, detgood) %>%
  group_by(method) %>%
  summarize(n=n(), TP=mean(gt & pres), TN=mean(!gt & !pres),
            FP = mean(!gt & pres), FN = mean(gt & !pres), P=sum(pres)) %>%
  mutate(sens = TP/(TP+FN), spec=TN/(TN+FP), ppv=TP/(TP+FP), npv=TN/(TN+FN))


# number of segments produced by each method ("to review")
filter(det, method=="harma" & par2==0.7 & thr==10 | 
         method=="med" & par2==60 & thr=="4.0" |
         method=="medlog" & par2==70 & thr==1.3) %>%
  bind_rows(detnewsegs) %>%
  group_by(method, file) %>%
  summarize(n(), sum(end-start))

# number of segments that would remain after review
det$start = as.numeric(det$start)
det$end = as.numeric(det$end)
detoverl = filter(det, method=="harma" & par2==0.7 & thr==10 | 
         method=="med" & par2==60 & thr=="4.0" |
         method=="medlog" & par2==70 & thr==1.3) %>%
  bind_rows(detnewsegs) %>%
  left_join(tieke, by=c("file" ="filename")) %>%
  filter(end.x>start.y & start.x<end.y)
detoverl = detoverl[,1:6]
detoverl = detoverl[!duplicated(detoverl),]

group_by(detoverl, method) %>%
  summarize(n=n())
