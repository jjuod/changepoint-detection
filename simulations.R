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
  
  lookback = floor(n*0.333)
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
  
  lookback = floor(n*0.333)
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
  
  lookback = floor(n*0.333)
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
ntotest = c(30, 60, 90, 130, 180, 240, 320, 440)
niter = 300

# at each n:
for(NPOINTS in ntotest){
  print(sprintf("Working on n=%d", NPOINTS))
  dets1 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, scenario=1, fBG=0.8)
  dets2 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, scenario=2, fBG=0.7)
  dets3 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, scenario=3, fBG=0.6)
  for(i in 1:niter){
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
write.table(detsalln, "~/Documents/gitrep/changepoint-detection/results-sim/sim1-estimates.tsv", quote=F, sep="\t", col.names=T, row.names=F)
write.table(allsegs, "~/Documents/gitrep/changepoint-detection/results-sim/sim1-detections.tsv", quote=F, sep="\t", col.names=T, row.names=F)


## PLOT RESULTS
## 1. test if the final estimate theta0 converges (FIGURE 1)

head(detsalln)
# fraction of background points to get n_B at each NPOINTS
detsalln$fBG = 0.8
temp1 = group_by(detsalln, n) %>%
  summarize(estimator="alg1", q2.5 = quantile(theta0, 0.025), q25=quantile(theta0, 0.25),
            q75=quantile(theta0, 0.75), q97.5=quantile(theta0, 0.975))
temp2 = group_by(detsalln, n) %>%
  summarize(estimator="median", q2.5 = quantile(median, 0.025), q25=quantile(median, 0.25),
            q75=quantile(median, 0.75), q97.5=quantile(median, 0.975))
temp3 = data.frame(n=ntotest) %>%
  mutate(estimator="mean", q2.5 = qnorm(0.025)/sqrt(fBG*n), q25=qnorm(0.25)/sqrt(fBG*n),
            q75=qnorm(0.75)/sqrt(fBG*n), q97.5=qnorm(0.975)/sqrt(fBG*n))

dets_sum = bind_rows(temp1, temp2, temp3) %>%
  group_by(n, estimator)

ungroup(dets_sum) %>%
  mutate(n=ifelse(estimator=="median", n, ifelse(estimator=="mean", n+4, n+8))) %>%
  ggplot() +
  geom_segment(aes(x=n, xend=n, col=estimator, y=q2.5, yend=q97.5), alpha=0.3, size=3) +
  geom_segment(aes(x=n, xend=n, col=estimator, y=q25, yend=q75), alpha=0.7, size=3) +
  ylab("quantile value") + theme_minimal()


## 2. test if changepoints are estimated consistently (TABLE 1)

# add the number of segments detected in each iteration
allsegs = group_by(allsegs, NPOINTS, i, run) %>%
  mutate(nsegs=n())

# distribution of nsegs for each method x sample size
distrnsegs = summarize(allsegs, nsegs=max(nsegs)) %>%
  group_by(NPOINTS, run) %>%
  summarize(freq1=sum(nsegs==1)/niter, freqMore=sum(nsegs>1)/niter) %>%
  mutate(freqLess=1-freq1-freqMore)

# see all segments raw
ungroup(allsegs) %>%
  mutate(xpos=NPOINTS+(run==2)*5+i*0.1, run=factor(run)) %>%
  ggplot(aes(x=xpos, col=run)) +
  geom_segment(aes(y=V1/NPOINTS, yend=V2/NPOINTS, xend=xpos), alpha=0.5) +
  geom_hline(yintercept=c(0.3, 0.5), col="green") +
  scale_color_manual(values=c("deepskyblue", "blue2")) +
  theme_minimal() + xlab("length of time series") + ylab("pos along the series") 

# show only the positions when 1 segment is reported
filter(allsegs, nsegs==1) %>%
  ungroup() %>%
  mutate(xpos=NPOINTS+(run==2)*5+i*0.1, run=factor(run)) %>%
  ggplot(aes(x=xpos, col=run)) +
  geom_segment(aes(y=V1/NPOINTS, yend=V2/NPOINTS, xend=xpos), alpha=0.3) +
  geom_step(aes(y=floor(0.3*NPOINTS+1)/NPOINTS), col="orange") +
  geom_hline(yintercept=0.5, col="orange") +
  geom_text(aes(x=NPOINTS+(run==2)*10, y=0.20, label=round(100*freqLess, 1)),
            col="black", size=4, data=distrnsegs) +
  geom_text(aes(x=NPOINTS+(run==2)*10, y=0.25, label=round(100*freqMore, 1)),
            col="black", size=4, data=distrnsegs) +
  scale_color_manual(values=c("deepskyblue", "blue2")) +
  theme_minimal() + xlab("length of time series") + ylab("pos along the series") 



###


niter = 30
dets = data.frame(iter=1:niter, ndets = 0, s=0, e=0, theta=0)
for(i in 1:niter){
  res = run(30)
  if(is.null(res$segs)){
    next
  } else if(!is.matrix(res$segs)){
    dets[i,] = c(i, 1, res$segs[1], res$segs[2], res$segs[3])
  } else {
    dets[i, "ndets"] = nrow(res$segs)
  }
}
dets


ggplot(dets) + geom_point(aes(x=X1, y=iter), col="blue") +
  geom_point(aes(x=X2, y=iter), col="blue") +
  theme_minimal()


## 1. optimizing runtime

run = function(n){
  ts = c(rnorm(10, 0, 1), rnorm(5, 3, 1), rnorm(15, 0, 1))
  ts = rep(ts, length.out=n)
  pen = autoset_penalty(ts)
  res = shortsgd(ts, 10, pen, 1)
  return(res)
}

mb = microbenchmark(n30 = run(30),
               n45 = run(45),
               n60 = run(60), 
               n90 = run(90),
               n120 = run(120), 
               n180 = run(180), 
               n240 = run(240), times=60)
mb
autoplot(mb) + scale_y_continuous()

# time here is in ns
mbdf = as.data.frame(mb)
mbdf$n = as.numeric(sub("n([0-9]*)", "\\1", mbdf$expr))
ggplot(mbdf, aes(x=n, y=time)) + geom_point() + geom_smooth(method="lm") +
  geom_smooth(method="lm", formula="y ~ I(x^2) + x", col="purple")
summary(lm(time ~ n, data=mbdf))
# mean runtime = -0.75 ms + 65 us * n






ts = c(rnorm(40, 0, 1), rnorm(5, 3, 1), rnorm(20, 0, 1), rnorm(5, 3, 1), rnorm(10, 0, 1))
plot(ts)

plot(ts)

res = shortsgd(ts, 10, pen, 1)
res$segs
plot(res$wt, type="l", col="grey40")
tail(res$wt)

lastw = rep(0, 50)
plot(NULL, ylim=c(-1,1), xlim=c(1, length(ts)), xlab="t", ylab=expression(w[t]))
for(i in 1:50){
  ts = c(rnorm(100, 0, 1), rnorm(5, 3, 1), rnorm(3000, 0, 1))

  res = shortsgd(ts, 5, 20, 1)
  lines(res$wt, col="grey40")
  lastw[i] = res$wt[length(ts)]
}
summary(lastw)
sd(lastw)
1/sqrt(3000) # SD of mean over realizations of 3000-long ts