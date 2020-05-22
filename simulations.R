# Simulations 1

library(dplyr)
library(tidyr)
library(ggplot2)
library(microbenchmark)
source("CLEAN-algorithm.R")


## For checking consistency and convergence (FIGURE 1 and TABLE 1):
# generate time series from 3 different scenarios,
# store the final estimates of w_t from our alg,
# median of the entire dataset for comparison w/ CAPA,
# and segment positions estimated w/ or w/o the last step.


## DATA GENERATION SCENARIOS
# Scenario 1: Gaussian w/ one strong segment
run_scen1 = function(n){
  l1 = floor(n*0.3)
  l2 = floor(n*0.2)
  l3 = floor(n*0.5)
  
  lookback = floor(n*0.5)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 3, 1), rnorm(l3, 0, 1))
  pen = autoset_penalty(ts)
  res = shortsgd(ts, lookback, pen, 1)
  finalwt = res$wt[length(res$wt)]
  res2 = shortfixed(ts, finalwt, lookback, pen, 1)
  
  return(list(wt = finalwt, med = median(ts), segs1 = res$segs, segs2 = res2$segs))
}

# Scenario 2: Gaussian w/ multiple weak segments
run_scen2 = function(n){
  l1 = floor(n*0.2)
  l2 = floor(n*0.1) # seg 0.2:0.3
  l3 = floor(n*0.2)
  l4 = floor(n*0.1) # seg 0.5:0.6
  l5 = floor(n*0.1)
  l6 = floor(n*0.1) # seg 0.7:0.8
  l7 = floor(n*0.2)
  
  lookback = floor(n*0.5)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, -1, 1),
         rnorm(l3, 0, 1), rnorm(l4, 1, 1),
         rnorm(l5, 0, 1), rnorm(l6, -1, 1),
         rnorm(l7, 0, 1))
  pen = autoset_penalty(ts)
  res = shortsgd(ts, lookback, pen, 1)
  finalwt = res$wt[length(res$wt)]
  res2 = shortfixed(ts, finalwt, lookback, pen, 1)
  
  return(list(wt = finalwt, med = median(ts), segs1 = res$segs, segs2 = res2$segs))
}

# Scenario 3: t-distribution w/ one large segment
run_scen3 = function(n){
  l1 = floor(n*0.2)
  l2 = floor(n*0.4)
  l3 = floor(n*0.4)
  
  lookback = floor(n*0.5)
  ts = c(rt(l1, 3)+0, rt(l2, 3)+2, rt(l3, 3)+0)
  pen = autoset_penalty(ts)
  res = shortsgd(ts, lookback, pen, sqrt(3))
  finalwt = res$wt[length(res$wt)]
  res2 = shortfixed(ts, finalwt, lookback, pen, sqrt(3))
  
  return(list(wt = finalwt, med = median(ts), segs1 = res$segs, segs2 = res2$segs))
}

## RUN THE SIMULATIONS:
detsalln = data.frame()
allsegs = matrix(0, ncol=7)
ntotest = c(30, 60, 90, 130, 180, 240, 320, 440, 550, 750)
niter = 500

# at each n:
for(NPOINTS in ntotest){
  print(sprintf("Working on n=%d", NPOINTS))
  dets1 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, scenario=1, fBG=0.8)
  dets2 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, scenario=2, fBG=0.7)
  dets3 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, scenario=3, fBG=0.6)
  for(i in 1:niter){
    print(i)
    res = run_scen1(NPOINTS)
    dets1[i, "theta0"] = res$wt
    dets1[i, "median"] = res$med
    
    allsegs = rbind(allsegs, cbind(res$segs1, NPOINTS, i, 1, 1))
    allsegs = rbind(allsegs, cbind(res$segs2, NPOINTS, i, 2, 1))
    
    res = run_scen2(NPOINTS)
    dets2[i, "theta0"] = res$wt
    dets2[i, "median"] = res$med
    
    allsegs = rbind(allsegs, cbind(res$segs1, NPOINTS, i, 1, 2))
    allsegs = rbind(allsegs, cbind(res$segs2, NPOINTS, i, 2, 2))
    
    res = run_scen3(NPOINTS)
    dets3[i, "theta0"] = res$wt
    dets3[i, "median"] = res$med
    
    allsegs = rbind(allsegs, cbind(res$segs1, NPOINTS, i, 1, 3))
    allsegs = rbind(allsegs, cbind(res$segs2, NPOINTS, i, 2, 3))
  }
  detsalln = bind_rows(detsalln, dets1)
  detsalln = bind_rows(detsalln, dets2)
  detsalln = bind_rows(detsalln, dets3)
}

head(detsalln)

# allsegs had an empty row due to initialization
head(allsegs)
allsegs = allsegs[-1,]
allsegs = data.frame(allsegs)
colnames(allsegs)[6] = "run"
colnames(allsegs)[7] = "scenario"

# save output:
write.table(detsalln, "results-sim/sim1-estimates.tsv", quote=F, sep="\t", col.names=T, row.names=F)
write.table(allsegs, "results-sim/sim1-detections.tsv", quote=F, sep="\t", col.names=T, row.names=F)

detsalln = read.table("results-sim/sim1-estimates.tsv", h=T)
allsegs = read.table("results-sim/sim1-detections.tsv", h=T)

## PLOT RESULTS
## 1. test if the final estimate theta0 converges (FIGURE 1)

head(detsalln)
temp1 = group_by(detsalln, n, scenario) %>%
  summarize(estimator="alg1", q2.5 = quantile(theta0, 0.025), q25=quantile(theta0, 0.25),
            q75=quantile(theta0, 0.75), q97.5=quantile(theta0, 0.975))
temp2 = group_by(detsalln, n, scenario) %>%
  summarize(estimator="median", q2.5 = quantile(median, 0.025), q25=quantile(median, 0.25),
            q75=quantile(median, 0.75), q97.5=quantile(median, 0.975))
temp3 = expand.grid(n=ntotest, scenario=1:3) %>%
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
  ylab("quantile value") + theme_bw(base_size = 15) + theme(panel.grid.major.x = element_blank())
ggsave("results-sim/fig1.png", width=23, height=12, units="cm", dpi=150)


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
# fraction of simulations reporting the right nubmer of segs
# distrnsegs %>%
#   ggplot(aes(x=NPOINTS)) + geom_line(aes(y=ncorr, lty=factor(run))) +
#   facet_wrap(~scenario) + theme_bw()

# mean absolute distance from each chp to closest estimated one
distrnsegs %>%
  ggplot(aes(x=NPOINTS)) + geom_line(aes(y=meand*NPOINTS, lty=factor(run))) +
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



###

set.seed(54321)
## 3. just a basic test

run_scen1(100)

mb = microbenchmark(n30 = run_scen1(30),
               n60 = run_scen1(60), 
               n120 = run_scen1(120), 
               n240 = run_scen1(240), 
               times=30)
mb
autoplot(mb) + scale_y_continuous()

# time here is in ns
mbdf = as.data.frame(mb)
mbdf$n = as.numeric(sub("n([0-9]*)", "\\1", mbdf$expr))
ggplot(mbdf, aes(x=n, y=time)) + geom_point() + geom_smooth(method="lm") +
  geom_smooth(method="lm", formula="y ~ I(x^2) + x", col="purple")
summary(lm(time ~ n, data=mbdf))
# mean runtime = -0.75 ms + 65 us * n


n = 100
l1 = floor(n*0.2)
l2 = floor(n*0.1)
l3 = floor(n*0.2)
l4 = floor(n*0.1)
l5 = floor(n*0.4)
ts = c(rnorm(l1, 0, 1), rnorm(l2, 2, 1), rnorm(l3, 5, 1), rnorm(l4, 2, 1), rnorm(l5, 0, 1))
plot(ts)

lookback = floor(n*0.33)
pen = autoset_penalty(ts)

res.known = shortfixed(ts, 0, lookback, pen, 1)
res.online = shortsgd(ts, lookback, pen, 1)
finalwt = res$wt[length(res$wt)]
res.plugin = shortfixed(ts, finalwt, lookback, pen, 1)

res.full = fulldetector_noprune_reference(ts, theta0=0, lookback, PEN=pen, PEN2=pen, SD=1)
res.full.opt = fulldetector_noprune(ts, theta0=0, lookback, PEN=pen, PEN2=pen, SD=1)
all(res.full$segs == res.full.opt$segs)
microbenchmark(ref=fulldetector_noprune_reference(ts, theta0=0, lookback, PEN=pen, PEN2=pen, SD=1), 
               opt=fulldetector_noprune(ts, theta0=0, lookback, PEN=pen, PEN2=pen, SD=1), times=10)


run_scen_new = function(n){
  l1 = floor(n*0.3)
  l2 = floor(n*0.2)
  l3 = floor(n*0.5)
  
  lookback = floor(n*0.5)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 3, 1), rnorm(l3, 0, 1))
  pen = autoset_penalty(ts)
  res = shortsgd(ts, lookback, pen, 1)

  res = detector(ts[1], length(ts), lookback, pen, 1)
  for(t in 2:length(ts)){
    res = detector.step(res, ts[t])
  }
  
  return(list(wt = res$wt[length(res$wt)], med = median(ts), segs1 = res$segs))
}

microbenchmark(run_scen_new(240), times=30)
