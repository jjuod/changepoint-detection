options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("~/Documents/gitrep/changepoint-detection/")
source("MAIN-algorithm.R")

#### ------ SIMULATIONS 1 --------
## For checking consistency and convergence (FIGURE 1 and TABLE 1):
# generate time series from 3 different scenarios,
# store the final estimates of w_t from our alg,
# median of the entire dataset for comparison w/ CAPA,
# and segment positions estimated w/ or w/o the last step.

## RUN THE SIMULATIONS:
# Recommend running that outside Rstudio, as older versions of it sometimes
# don't manage to do garbage collection efficiently
#system("Rscript simulations1-separate.R")

ntotest = c(30, 60, 90, 130, 180, 240, 320, 440, 550, 750)
niter = 500

# read output:
detsalln = read.table("../drafts/changepoint-method/results-sim-new/sim1-estimates.tsv", h=T)
allsegs = read.table("../drafts/changepoint-method/results-sim-new/sim1-detections.tsv", h=T)

## PLOT RESULTS
## 1. test if the final estimate theta0 converges (FIGURE 1)

head(detsalln)
temp1 = group_by(detsalln, n, scenario) %>%
  summarize(estimator="alg1", q2.5 = quantile(theta0, 0.05), q25=quantile(theta0, 0.25),
            q75=quantile(theta0, 0.75), q97.5=quantile(theta0, 0.95))
temp2 = group_by(detsalln, n, scenario) %>%
  summarize(estimator="anomaly", q2.5 = quantile(median, 0.05), q25=quantile(median, 0.25),
            q75=quantile(median, 0.75), q97.5=quantile(median, 0.95))
temp3 = group_by(detsalln, n, scenario) %>%
    summarize(estimator="aPELT", q2.5 = quantile(theta.apelt, 0.05), q25=quantile(theta.apelt, 0.25),
              q75=quantile(theta.apelt, 0.75), q97.5=quantile(theta.apelt, 0.95))
temp4 = detsalln[,c("n", "scenario", "fBG")] %>%
  unique() %>%
  mutate(se = ifelse(scenario==3, sqrt(3), 1) / sqrt(fBG*n) ) %>%
  mutate(estimator="oracle", q2.5 = qnorm(0.05, sd=se), q25=qnorm(0.25, sd=se),
            q75=qnorm(0.75, sd=se), q97.5=qnorm(0.95, sd=se))

dets_sum = bind_rows(temp1, temp2, temp3, temp4)

ungroup(dets_sum) %>%
  filter(n<700) %>%
  mutate(scenario=factor(scenario, labels=c("one segm.", "multiple", "heavy tail"))) %>%
  ggplot() +  geom_hline(yintercept=0, col="grey40") + 
  geom_linerange(aes(x=n, col=estimator, ymin=q2.5, ymax=q97.5), alpha=0.5, size=2.5, position=position_dodge(20)) +
  geom_linerange(aes(x=n, col=estimator, ymin=q25, ymax=q75), alpha=0.9, size=2.5, position=position_dodge(20)) +
  facet_grid(scenario~., scales = "free_y") + 
  scale_x_continuous(breaks=ntotest, minor_breaks = NULL) +
  scale_color_brewer(palette="Set1") +
  ylab("quantile value") + theme_bw(base_size = 15) +
  theme(panel.grid.major.x=element_blank(), legend.text=element_text(size=13),
          legend.position="bottom")
# ggsave("../drafts/changepoint-method/results-sim/fig1.png", width=23, height=12, units="cm", dpi=150)
ggsave("../drafts/changepoint-method/v21dec/fig1.eps", width=21, height=12, units="cm", device=cairo_ps, fallback_resolution=600)


## 2. test if changepoints are estimated consistently (TABLE 1)

# add the number of segments detected in each iteration
allsegs = group_by(allsegs, NPOINTS, run, scenario, i) %>%
  mutate(nsegs=n())

# extracts distances to true chps for TPR
dist = function(pos1, pos2, scen, n){
  n = n[1]
  scen = scen[1]
  if(scen==1){
    truepos = floor(c(0.3*n+1, 0.5*n))
  } else if(scen==2){
    truepos = floor(c(0.2*n+1, 0.3*n, 0.5*n+1, 0.6*n, 0.7*n+1, 0.8*n))
  } else if(scen==3){
    truepos = floor(c(0.2*n+1, 0.6*n))
  }
  # for each true segment, calculate distance to chp, in points:
  d = rep(Inf, length(truepos))
  for(pairix in seq(1, length(truepos), 2)){
      chps = truepos[pairix]
      chpe = truepos[pairix+1]
      # pick the closest estimated chp for start and end columns
      d1 = min(abs(pos1-chps))
      d2 = min(abs(pos2-chpe))
      d[pairix] = d1
      d[pairix+1] = d2
  }
  # d<thr <=> there was a TP for this true chp within thr
  return(d)
}

allsegs = filter(allsegs, !is.na(V1))  # Drop anom/aPELT iters w/o detections

# get distances to the nearest true chp
segsperiter = group_by(allsegs, run, scenario, NPOINTS, i) %>%
    summarize(nsegs=max(nsegs), disterr=dist(V1, V2, scenario, NPOINTS)) %>%
    mutate(disterr=disterr/NPOINTS)

# fill in missing rows when an iteration returns 0 segs
temp = expand.grid(NPOINTS=ntotest, run=1:4, scenario=1:3, i=1:niter)
temp$nsegs = 0
temp$disterr = Inf
segsperiter = anti_join(temp, segsperiter, by=c("NPOINTS", "run", "scenario", "i")) %>% 
  bind_rows(segsperiter, .)

# calculate fraction of simulations reporting a chp within 0.05 of true seg ("TPR")
tpr_table = group_by(segsperiter, NPOINTS, scenario, run, i) %>%
  summarize(tpr = max(disterr)<0.05) %>%
  summarize(tpr = sum(tpr)/niter) %>%  # per penalty level
  ungroup

# Plot:
tpr_table %>%
    ggplot(aes(x=NPOINTS)) + geom_line(aes(y=tpr, lty=factor(run))) +
    facet_wrap(~scenario) + theme_bw()

# mean number of segs reported for each n x scenario x run
distrnsegs = group_by(allsegs, run, scenario, NPOINTS, i) %>%
    summarize(nsegs=max(nsegs)) %>%
    mutate(ntrue = ifelse(scenario==2, 3, 1)) %>%
    summarize(meannseg=mean(nsegs), ntrue=max(ntrue),
              ncorr=sum(nsegs==ntrue)/niter)  # add up iterations

distrnsegs %>%
  ggplot(aes(x=NPOINTS)) + geom_line(aes(y=meannseg, lty=factor(run))) +
  geom_hline(aes(yintercept=ntrue), col="green") +
  facet_wrap(~scenario) + theme_bw()
distrnsegs %>%
    filter(run!=1) %>%  # drop no-repeat
    mutate(alg=factor(run, labels=c("algorithm 1", "anomaly", "aPELT"))) %>%
    ggplot(aes(x=NPOINTS)) +
    geom_hline(aes(yintercept=0), col="grey30") +
    geom_line(aes(y=meannseg-ntrue, col=alg, lty=alg), lwd=1) +
    scale_linetype_manual(name="detections:", values=c(1,3,2)) +
    scale_color_manual(name="detections:", values=RColorBrewer::brewer.pal(3, "Set1")) +
    facet_wrap(~scenario, labeller=labeller(scenario=c("1"="one segment", "2"="multiple", "3"="heavy tail")),
               scales="free_y") +
    theme_bw() + xlab("n") + ylab(expression(E(~hat(k)-k))) +
    theme(legend.position = "bottom", text=element_text(size=14),
          legend.text=element_text(size=13), plot.margin = unit(c(0.1,0.3,0,0.1), "cm"))
ggsave("../drafts/changepoint-method/v21dec/fig1b.eps", width=17, height=8, units="cm")

# see raw segments for one case
filter(allsegs, scenario==2, run==2, i<=100) %>%
  ungroup() %>%
  ggplot(aes(y=i)) +
  geom_vline(xintercept=c(0.2, 0.3, 0.5, 0.6, 0.7, 0.8), col="green2", lwd=1) +
    geom_point(aes(x=V1/NPOINTS), alpha=0.6) + geom_point(aes(x=V2/NPOINTS), alpha=0.6) + 
  facet_wrap(~NPOINTS) +
  theme_minimal()
  
# Table 1:
# remember that run2 is the correct (full) run!
t1 = distrnsegs[,c("scenario", "NPOINTS", "run", "meannseg")] %>%
  spread(key="run", value="meannseg")
t2 = tpr_table[,c("scenario", "NPOINTS", "run", "tpr")] %>%
  spread(key="run", value="tpr")
colnames(t1)[3:6] = c("meann.norepeat", "meann.full", "meann.anom", "meann.apelt")
colnames(t2)[3:6] = c("tpr.norepeat", "tpr.full", "tpr.anom", "tpr.apelt")
full_join(t1, t2, by=c("scenario", "NPOINTS")) %>%
  filter(NPOINTS %in% c(30, 90, 180, 440, 750)) %>%
  select(-ends_with(c("norepeat"))) %>%
  print.data.frame
# Suppl. Table:
full_join(t1, t2, by=c("scenario", "NPOINTS")) %>%
  filter(NPOINTS %in% c(30, 90, 180, 440, 750)) %>%
  select(-ends_with(c("anom", "apelt")))
  
# Table S2:
# Cost (mean and SD) estimated in each setting, for the full Alg 1, or the online version,
# or using theoretical segment positions.
tableS2 = group_by(detsalln, n, scenario) %>%
  summarize(mean(likEst2), sd(likEst2), mean(likEst1), sd(likEst1), mean(likTheor), sd(likTheor))
tableS2 %>%
  filter(n %in% c(30, 90, 180, 440, 750)) %>%
  arrange(scenario, n) %>%
  print.data.frame(digits=3)


#### ------ SIMULATIONS 2 --------
# main processing oved to simulations2-parsing.R
# this is just leftover testing scripts

# see raw segments for one case
# filter(allsegs2, NPOINTS==60, scen==2) %>%
#   ggplot(aes(y=i+0.2*(segtype=="seg"), col=segtype)) +
#   geom_vline(xintercept=c(0.2, 0.3, 0.5, 0.7), col="green2", lwd=1) +
#   geom_segment(aes(x=(V1-1)/NPOINTS, xend=V2/NPOINTS, yend=i+0.2*(segtype=="seg")), alpha=0.8, lwd=1) +
#   facet_wrap(~alg) + scale_color_manual(values=c("black", "blue")) +
#   theme_minimal()
filter(allsegs2, NPOINTS==220, scen==1, alg=="sparse", i<100) %>%
  ggplot(aes(y=i+0.2*(segtype=="seg"), col=segtype)) +
  geom_vline(xintercept=c(0.2, 0.3, 0.5, 0.7), col="green2", lwd=1) +
  geom_segment(aes(x=(V1-1)/NPOINTS, xend=V2/NPOINTS, yend=i+0.2*(segtype=="seg")), alpha=0.8, lwd=1) +
  facet_wrap(~alg) + scale_color_manual(values=c("black", "blue")) +
  theme_minimal()
filter(allsegs2, NPOINTS==60, scen==1, alg=="full", i<100) %>%
    ggplot(aes(y=i+0.2*(segtype=="seg"), col=segtype)) +
    geom_vline(xintercept=c(0.2, 0.3, 0.5, 0.7), col="green2", lwd=1) +
    geom_segment(aes(x=(V1-1)/NPOINTS, xend=V2/NPOINTS, yend=i+0.2*(segtype=="seg")), alpha=0.8, lwd=1) +
    facet_wrap(~alg) + scale_color_manual(values=c("black", "blue")) +
    theme_minimal()

filter(allsegs2, NPOINTS==60, scen==3, alg %in% c("full", "anomaly"), i<60) %>%
    ggplot(aes(y=i+0.2*(segtype=="seg"), col=segtype)) +
    geom_vline(xintercept=c(0.2, 0.3, 0.5, 0.7), col="green2", lwd=1) +
    geom_segment(aes(x=(V1-1)/NPOINTS, xend=V2/NPOINTS, yend=i+0.2*(segtype=="seg")), alpha=0.8, lwd=1) +
    facet_wrap(~alg) + scale_color_manual(values=c("black", "blue")) +
    theme_minimal()

# NOT reports [start-1, end] positions it seems, so this fix may be needed
# allsegs2$V1[allsegs2$alg=="NOT"] = allsegs2$V1[allsegs2$alg=="NOT"]+1
# allsegs2$V2[allsegs2$alg=="NOT"] = allsegs2$V2[allsegs2$alg=="NOT"]+1

# extracts distances to true chps for TPR, ignoring the structure
get_dist_errors2 = function(pos1, pos2, scen, n, tpdf){
    truepos = tpdf$truepos[tpdf$scen==scen[1] & tpdf$n==n[1]]
    
    # for each true segment, calculate distance to chp, in points:
    d = rep(Inf, length(truepos))
    for(pairix in seq(length(truepos))){
        chp = truepos[pairix]
        # pick the closest estimated chp from either start or end columns
        d1 = min(abs(pos1-chp))
        d2 = min(abs(pos2-chp))
        d[pairix] = min(d1, d2)
    }
    # d<thr <=> there was a TP for this true chp within thr
    return(d)
}
# Same but checks start vs start and end vs end
get_dist_errors2_ordered = function(pos1, pos2, scen, n, tpdf){
    truepos = tpdf$truepos[tpdf$scen==scen[1] & tpdf$n==n[1]]
    
    # for each true segment, calculate distance to chp, in points:
    d = rep(Inf, length(truepos))
    for(pairix in seq(1, length(truepos), 2)){
        chps = truepos[pairix]
        chpe = truepos[pairix+1]
        # pick the closest estimated chp for start and end columns
        d1 = min(abs(pos1-chps))
        d2 = min(abs(pos2-chpe))
        d[pairix] = d1
        d[pairix+1] = d2
    }
    # d<thr <=> there was a TP for this true chp within thr
    return(d)
}
# Same but only reports one output per pair of chps (max error)
get_dist_errors2_paired = function(pos1, pos2, scen, n, tpdf){
    truepos = tpdf$truepos[tpdf$scen==scen[1] & tpdf$n==n[1]]
    
    # for each true segment, calculate distance to chp, in points:
    d = rep(Inf, length(truepos))
    for(pairix in seq(1, length(truepos), 2)){
        chps = truepos[pairix]
        chpe = truepos[pairix+1]
        # pick the closest estimated chp for start and end columns
        errs1 = abs(pos1-chps)
        errs2 = abs(pos2-chpe)
        ds = pmax(errs1, errs2) # distance(true seg, est seg) := max(err left, err right)
        ds = min(ds)  # minimum over est segments
        d[pairix] = ds
        d[pairix+1] = ds
    }
    # d<thr <=> there was a TP for this true chp within thr
    return(d)
}
# Same as 1st but reports the errors per each detected chp, not per true chp
get_dist_errors_est = function(pos1, pos2, scen, n, tpdf, segtype){
    d = rep(Inf, length(pos1)*2)
    if (segtype[1]=="nuis"){ return(d)}
    truepos = tpdf$truepos[tpdf$scen==scen[1] & tpdf$n==n[1]]
    
    # for each true segment, calculate distance to chp, in points:
    for(pairix in seq(length(pos1))){
        chp1 = pos1[pairix]
        chp2 = pos2[pairix]
        # pick the closest true chp for each estimated one
        d1 = min(abs(truepos-chp1))
        d2 = min(abs(truepos-chp2))
        d[2*pairix-1] = d1
        d[2*pairix] = d2
    }
    # d<thr <=> there was a TP for this est chp within thr
    return(d)
}

