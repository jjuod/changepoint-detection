options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
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
system("Rscript simulations1-separate.R")

ntotest = c(30, 60, 90, 130, 180, 240, 320, 440, 550, 750)
niter = 500

# read output:
detsalln = read.table("../drafts/changepoint-method/results-sim/sim1-estimates.tsv", h=T)
allsegs = read.table("../drafts/changepoint-method/results-sim/sim1-detections.tsv", h=T)

## PLOT RESULTS
## 1. test if the final estimate theta0 converges (FIGURE 1)

head(detsalln)
temp1 = group_by(detsalln, n, scenario) %>%
  summarize(estimator="alg1", q2.5 = quantile(theta0, 0.025), q25=quantile(theta0, 0.25),
            q75=quantile(theta0, 0.75), q97.5=quantile(theta0, 0.975))
temp2 = group_by(detsalln, n, scenario) %>%
  summarize(estimator="median", q2.5 = quantile(median, 0.025), q25=quantile(median, 0.25),
            q75=quantile(median, 0.75), q97.5=quantile(median, 0.975))
temp3 = detsalln[,c("n", "scenario", "fBG")] %>%
  unique() %>%
  mutate(se = ifelse(scenario==3, sqrt(3), 1) / sqrt(fBG*n) ) %>%
  mutate(estimator="mean", q2.5 = qnorm(0.025, sd=se), q25=qnorm(0.25, sd=se),
            q75=qnorm(0.75, sd=se), q97.5=qnorm(0.975, sd=se))

dets_sum = bind_rows(temp1, temp2, temp3)

ungroup(dets_sum) %>%
  filter(n<700) %>%
  ggplot() +
  geom_linerange(aes(x=n, col=estimator, ymin=q2.5, ymax=q97.5), alpha=0.4, size=2, position=position_dodge(15)) +
  geom_linerange(aes(x=n, col=estimator, ymin=q25, ymax=q75), alpha=0.8, size=2, position=position_dodge(15)) +
  facet_grid(scenario~., scales = "free_y", labeller="label_both") + 
  scale_x_continuous(breaks=ntotest, minor_breaks = NULL) +
  ylab("quantile value") + theme_bw(base_size = 15) + theme(panel.grid.major.x = element_blank(), legend.text=element_text(size=13))
# ggsave("../drafts/changepoint-method/results-sim/fig1.png", width=23, height=12, units="cm", dpi=150)
ggsave("../drafts/changepoint-method/results-sim/fig1.eps", width=23, height=12, units="cm", device=cairo_ps, fallback_resolution=600)


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
  maxd = 0
  # for each true chp:
  for(chp in truepos){
    # pick the closest estimated chp
    d1 = min(abs(pos1-chp))/n
    d2 = min(abs(pos2-chp))/n
    d = min(d1, d2)
    # get the worst-case (i.e. max over all truepos) d
    maxd = max(maxd, d)
  }
  # max d<thr <=> there was a TP for each true chp within thr
  return(maxd)
}

# summary of each iteration (dist for TPR and nsegs)
segsperiter = summarize(allsegs, nsegs=max(nsegs),
                       maxd=dist(V1, V2, scenario, NPOINTS))
# fill in missing rows when an iteration returns 0 segs
temp = expand.grid(NPOINTS=ntotest, run=1:2, scenario=1:3, i=1:500)
temp$nsegs = 0
temp$maxd = 1
segsperiter = anti_join(temp, segsperiter, by=c("NPOINTS", "run", "scenario", "i")) %>% 
  bind_rows(segsperiter, .)

# overall summary for each method
distrnsegs = group_by(segsperiter, NPOINTS, scenario, run) %>%
  mutate(ntrue = ifelse(scenario==2, 3, 1)) %>%
  summarize(meannseg = mean(nsegs), ncorr=sum(nsegs==ntrue)/max(i),
            tpr = mean(maxd<0.05), meand = mean(maxd), ntrue=max(ntrue)) %>%
  ungroup

# mean number of segs reported for each n x scenario x run
distrnsegs %>%
  ggplot(aes(x=NPOINTS)) + geom_line(aes(y=meannseg, lty=factor(run))) +
  geom_hline(aes(yintercept=ntrue), col="green") +
  facet_wrap(~scenario) + theme_bw()

# fraction of simulations reporting a chp within 0.05 of true seg ("TPR")
distrnsegs %>%
  ggplot(aes(x=NPOINTS)) + geom_line(aes(y=tpr, lty=factor(run))) +
  facet_wrap(~scenario) + theme_bw()

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
t2 = distrnsegs[,c("scenario", "NPOINTS", "run", "tpr")] %>%
  spread(key="run", value="tpr")
colnames(t1)[3:4] = c("meann.run1", "meann.run2")
colnames(t2)[3:4] = c("tpr.run1", "tpr.run2")
full_join(t1, t2, by=c("scenario", "NPOINTS")) %>%
  filter(NPOINTS %in% c(30, 90, 180, 440, 750)) %>%
  print.data.frame


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
source("MAIN-algorithm.R")

## RUN THE SIMULATIONS:
# Recommend running that outside Rstudio, as older versions of it sometimes
# don't manage to do garbage collection efficiently
system("Rscript simulations2-separate.R")

# read output:
allsegs2 = read.table("../drafts/changepoint-method/results-sim/sim2-detections.tsv", h=T)

ntotest2 = c(30, 60, 100, 150, 220)
niter2 = 500

## PLOT RESULTS

# NOT needs post-processing as it does not distinguish background and segments
# so we arbitrarily define any segs w/ mean within 1 SD of theta as background:
allsegs2 = filter(allsegs2, alg!="NOT" | abs(V3)>1)

# some summaries
nrow(allsegs2)
group_by(allsegs2, NPOINTS, alg, scen) %>% summarize(n())

## Does pruning make any difference?
dfF = filter(allsegs2, alg=="full")
dfP = filter(allsegs2, alg=="pruned")
diffis = bind_rows(anti_join(dfF, dfP, by=c("NPOINTS", "i", "segtype", "scen", "V1", "V2")),
                   anti_join(dfP, dfF, by=c("NPOINTS", "i", "segtype", "scen", "V1", "V2")))

# Table S2:
# all iterations that differed
semi_join(allsegs2, diffis, by=c("NPOINTS", "scen", "i")) %>%
  filter(alg %in% c("full", "pruned"), scen %in% 1:2) %>%   # skip scenario 3 as that one has more differences
  mutate(segtype=ifelse(segtype=="nuis", "N", "S")) %>%
  .[,c("alg", "scen", "NPOINTS", "i", "segtype", "V1", "V2")]

# how frequently something differed?
diffis_u = unique(diffis[,c("NPOINTS", "scen", "i")])
table(diffis_u$scen)
nrow(diffis_u) / niter2 / length(ntotest2) / 3  # 3 scenarios

# add the number of segments detected in each iteration
allsegs2 = group_by(allsegs2, NPOINTS, alg, segtype, i, scen) %>%
  mutate(nsegs=n())

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


# extracts distances to true chps for TPR
dist2 = function(pos1, pos2, segtype, scen, n){
  n = n[1]
  segtype = segtype[1]
  if(scen[1]==1){
    if(segtype=="seg"){
      truepos = floor(c(0.3*n+1, 0.5*n))
    } else {
      truepos = floor(c(0.2*n+1, 0.7*n))
    }  
  } else if(scen[1]==2){
    if(segtype=="seg"){
      truepos = floor(c(0.5*n+1, 0.6*n, 0.7*n+1, 0.8*n))
    } else {
      truepos = floor(c(0.2*n+1, 0.4*n))
    }
  } else {
    l1 = floor(n*0.1)
    l2 = floor(n*0.05)
    l3 = ceiling(n*0.05)
    # Note that the starts and ends are interleaved to maintain pairs
    truepos = c(rbind(l1+1 + 0:8*(l2+l3),  # starts
                      l1 + l2 + 0:8*(l2+l3)))  # ends
  }
  
  maxd = 0
  # for each true segment:
  for(pairix in seq(1, length(truepos))){
    chpst = truepos[pairix]
    # chpe = truepos[pairix+1]
    # pick the two estimated chps closest to true start and end
    d1 = min(abs(pos1-chpst))/n
    d2 = min(abs(pos2-chpst))/n
    d = min(d1, d2)
    # get the worst-case (i.e. max over all pairs) d
    maxd = max(maxd, d)
  }
  # max d<thr <=> there was a TP for each true chp within thr
  return(maxd)
}

# summary of each iteration (dist for TPR and nsegs)
segsperiter = summarize(allsegs2, nsegs=max(nsegs),
                        maxd=dist2(V1, V2, segtype, scen, NPOINTS))
# fill in missing rows when an iteration returns 0 segs
temp = expand.grid(NPOINTS=ntotest2, alg=c("pruned", "full", "anomaly", "NOT", "aPELT", "sparse"), segtype=c("seg", "nuis"), i=1:niter2, scen=1:3)
temp$nsegs = 0
temp$maxd = 1
segsperiter = anti_join(temp, segsperiter, by=c("NPOINTS", "alg", "segtype", "scen", "i")) %>% 
  bind_rows(segsperiter, .)

# overall summary for each method
distrnsegs = group_by(segsperiter, NPOINTS, alg, segtype, scen) %>%
  mutate(ntrue = ifelse(scen==3, ifelse(segtype=="seg", 9, 0),
                  ifelse(scen==2, ifelse(segtype=="seg", 2, 1),
                         1))) %>%
  summarize(meannseg = mean(nsegs), ncorr=sum(nsegs==ntrue)/max(i),
            tpr = mean(maxd<0.05), meand = mean(maxd), ntrue=max(ntrue)) %>%
  ungroup

# mean number of segs reported for each n x scenario x run
distrnsegs %>%
  filter(alg!="full", segtype=="seg" | alg=="pruned") %>%  # no nuisance segs in standard chp algs
  ggplot(aes(x=NPOINTS)) +
  geom_hline(aes(yintercept=ntrue, lty=segtype), col="black") +
  geom_line(aes(y=meannseg, col=factor(alg), lty=segtype), lwd=1) +
  scale_linetype_discrete(limits=c("seg", "nuis"), labels=c("signal", "nuisance"), name="ground truth") + 
  scale_color_discrete(labels=c("aPELT", "anomaly", "NOT", "sparse", "proposed"),
                       breaks=c("aPELT", "anomaly", "NOT", "sparse", "pruned"), name="detections") + 
  facet_wrap(~scen, labeller = labeller(scen=c("1"="scenario 1", "2"="scenario 2", "3"="scenario 3"))) + theme_bw()

# fraction of simulations reporting a chp within 0.05 of true seg ("TPR")
distrnsegs %>%
  ggplot(aes(x=NPOINTS)) + geom_line(aes(y=tpr, col=factor(alg))) +
  facet_grid(scen~segtype) + theme_bw()

# Table 2:
t1 = distrnsegs[,c("scen", "NPOINTS", "alg", "meannseg", "segtype")] %>%
  filter(alg!="full", segtype=="seg" | alg=="pruned") %>%  # no nuisance segs in standard chp algs
  mutate(algseg = paste(alg, segtype, sep="_")) %>%
  select(-one_of(c("alg", "segtype"))) %>%
  spread(key="algseg", value="meannseg") %>%
  select("scen", "NPOINTS", "pruned_seg", "pruned_nuis", "anomaly_seg", "aPELT_seg", "sparse_seg", "NOT_seg")
t2 = distrnsegs[,c("scen", "NPOINTS", "alg", "tpr", "segtype")] %>%
  filter(alg!="full", segtype=="seg" | alg=="pruned") %>%  # no nuisance segs in standard chp algs
  mutate(algseg = paste(alg, segtype, sep="_")) %>%
  select(-one_of(c("alg", "segtype"))) %>%
  spread(key="algseg", value="tpr") %>%
  select("scen", "NPOINTS", "pruned_seg", "pruned_nuis", "anomaly_seg", "aPELT_seg", "sparse_seg", "NOT_seg")
print.data.frame(t1)
print.data.frame(t2)

# Alternatively:
# Figure 2: bias in the number of segs reported for each n x scenario x run
distrnsegs %>%
  filter(alg!="full", segtype=="seg" | alg=="pruned") %>%  # no nuisance segs in standard chp algs
  mutate(bias = meannseg - ntrue, algseg=factor(ifelse(segtype=="seg", alg, "nuisance"),
                                          levels=c("anomaly", "aPELT", "sparse", "NOT", "pruned", "nuisance"),
                                          labels=c("anomaly", "aPELT", "sparse", "not", "proposed", "proposed (nuisance)"))) %>%
  ggplot(aes(x=NPOINTS)) +
  geom_hline(aes(yintercept=0), col="grey30") +
  geom_line(aes(y=bias, col=algseg, lty=algseg), lwd=1) +
  geom_point(aes(y=bias, shape=algseg, col=algseg, alpha=algseg), size=3) +
  scale_alpha_manual(name="detections:", values=c("anomaly"=0, "aPELT"=1, "sparse"=1, "not"=1, "proposed"=0, "proposed (nuisance)"=0)) +
  scale_shape_manual(name="detections:", values=c("anomaly"=17, "aPELT"=4, "sparse"=15, "not"=20, "proposed"=1, "proposed (nuisance)"=1)) +
  scale_linetype_manual(values=c("anomaly"=1, "aPELT"=1, "sparse"=1, "not"=1, "proposed"=1, "proposed (nuisance)"=3),
                        name="detections:") +
  scale_color_manual(name="detections:", values =c("dodgerblue2", "blue4", "mediumturquoise", "skyblue1", "coral1", "coral1")) +
  facet_wrap(~scen, labeller = labeller(scen=c("1"="scenario 1", "2"="scenario 2", "3"="scenario 3")), scales="free_y") +
  theme_bw() + xlab("n") + ylab(expression(E(~hat(k)-k))) +
  theme(legend.position = "bottom", text=element_text(size=14), legend.text=element_text(size=13), plot.margin = unit(c(0.1,0.3,0,0.1), "cm"))
# ggsave("../drafts/changepoint-method/results-sim/fig2.png", width=17, height=10, units="cm", dpi=150)
ggsave("../drafts/changepoint-method/results-sim/fig2.eps", width=17, height=10, units="cm")


# Compare theta_i for the first true segment:
# (using the signal segment w/ most overlap)
best1seg = filter(allsegs2, scen==1) %>%
  group_by(alg, segtype, NPOINTS, i) %>%
  mutate(firstend = pmin(V2, floor(0.5*NPOINTS)),
         laststart = pmax(V1, floor(0.3*NPOINTS+1))) %>%
  mutate(overl = pmax(firstend-laststart, 0)) %>%
  top_n(1, rank(overl, ties.method = "r"))
extract_theta = function(V1, V2, V3, segtype){
  if("nuis" %in% segtype & "seg" %in% segtype){
    ni = segtype=="nuis"
    si = segtype=="seg"
    # if segment overlaps a nuisance:
    if(V1[ni]<=V1[si] & V2[ni]>=V2[si]){
      thetaN = V3[segtype=="nuis"]
      thetaS = V3[segtype=="seg"] - thetaN  
    } else{
      thetaS = V3[si] 
    }
  } else if ("seg" %in% segtype){
    thetaS = V3[segtype=="seg"]
  } else {
    thetaS = NA
  }
  return(thetaS)
}
group_by(best1seg, alg, NPOINTS, i) %>%
  summarize(thetaS = extract_theta(V1, V2, V3, segtype)) %>%
  mutate(thetaS = ifelse(alg=="anomaly", sqrt(thetaS), thetaS)) %>%  # anomaly reports squared estimate
  summarize(m=mean(thetaS, na.rm=T), s=sd(thetaS, na.rm=T), n=sum(!is.na(thetaS))) %>%
  print.data.frame

