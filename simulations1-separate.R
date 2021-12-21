options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(anomaly)
source("MAIN-algorithm.R")


# Our implementation of aPELT-profile from Zhao & Yau (2019).
# Runs PELT with a range of background parameter values, and chooses the one with min cost.
# Simplified to avoid the version which starts with non-background segments
# as that is not true in our simulation scenarios, and we can reuse our PELT function then.
# Also allows limiting the max length of accepted segments.
aPELT_profile = function(ts, MAXLOOKBACK, PEN, SD){
    opt = optim(par=mean(ts), fn=function(theta0) shortfixed(ts, theta0, MAXLOOKBACK, PEN, SD)$cost,
                control=list(warn.1d.NelderMead=F))
    cat(paste("Took", opt$counts[1], "evaluations to optimize\n"))
    cat(paste("Final theta: ", opt$par, "\n"))
    res = shortfixed(ts, opt$par, MAXLOOKBACK, PEN, SD)
    return(list(segs=res$segs, theta=opt$par))
}


#### ------ SIMULATIONS 1 --------
## For checking consistency and convergence (FIGURE 1 and TABLE 1):
# generate time series from 3 different scenarios,
# store the final estimates of w_t from our alg,
# median of the entire dataset for comparison w/ CAPA,
# and segment positions estimated w/ or w/o the last step.

## Calculate log likelihood (i.e. final cost) for a segmentation
## assuming a normal model with given background mean & SD.
#  - df: a dataframe of only the segments detected (first and second col=start and end)
#  - ts: original time series
#  - pen: penalty for each segment
getLik = function(df, ts, sigma0, pen, mu0=FALSE){
  nparams = 0
  logl = 0
  
  isBg = rep(T, length(ts))
  
  if(nrow(df)>0){
    for(i in 1:nrow(df)){
      start = df[i,1]
      end = df[i,2]
      # mask any S segment points so that they would be excluded
      # from bg likelihood later
      isBg[start:end] = F
      
      nparams = nparams+1
      
      # get cost of a segment
      segx = ts[start:end]
      logl = logl - sum(dnorm(segx, mean(segx), sigma0, log=T))
    }
  }
  
  # deal with bg points
  bgpoints = ts[isBg]
  if(mu0){
    bglogl = -sum(dnorm(bgpoints, mu0, sigma0, log=T))  
  } else {
    bglogl = -sum(dnorm(bgpoints, mean(bgpoints), sigma0, log=T))
  }
  
  # num parameters = num segments x3 (for mean and two endpoints in each) + 2
  klogn = nrow(df)*pen
  cat(sprintf("Log-likelihood of background points: %.1f\n", bglogl))
  cat(sprintf("Log-likelihood of non-bg points: %.1f\n", logl))
  cat(sprintf("penalty k log(n): %.1f\n", klogn))
  return(logl+bglogl+klogn)
}

apply_all_detectors = function(ts, lookback, pen, SD=1){
    res = shortsgd(ts, lookback, pen, SD)
    finalwt = res$wt[length(res$wt)]
    res2 = shortfixed(ts, finalwt, lookback, pen, SD)
    
    # classical epidemic detector implemented in R:anomaly.
    med = median(ts)
    ts.norm = (ts-med)/SD
    res.anom = capa.uv(ts.norm, beta=pen, beta_tilde=pen, type="mean", min_seg_len=2,
                       max_seg_len=lookback, transform=identity)
    res.anom.c = collective_anomalies(res.anom)[,1:3]
    res.anom.p = point_anomalies(res.anom)
    # attach length 1 anomalies, if any
    if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
        res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location,
                                mean.change=res.anom.p$strength) 
        res.anom.c = rbind(res.anom.c, res.anom.p)
    }
    if(nrow(res.anom.c)==0) res.anom.c = matrix(c(NA, NA, NA), nrow=1)
    res.anom.c = as.matrix(res.anom.c)
    
    # aPELT_profile detector, our implementation
    res.apelt = aPELT_profile(ts, MAXLOOKBACK=lookback, PEN=pen, SD=SD)
    theta.apelt = res.apelt$theta
    res.apelt = res.apelt$segs
    if(nrow(res.apelt)==0) res.apelt = matrix(c(NA, NA, NA), nrow=1)

    # costs:
    # using the estimated bg parameter and segments, no step 13
    likEst1 = getLik(res$segs, ts, 1, pen, finalwt)
    # using the estimated bg parameter and segments, full alg
    likEst2 = getLik(res2$segs, ts, 1, pen, finalwt)
    
    return(list(wt = finalwt, med = med, theta.apelt = theta.apelt, # all the bg estimates
                segs1 = res$segs, segs2 = res2$segs,
                segsan = res.anom.c, segsap = res.apelt,
                likEst1=likEst1, likEst2=likEst2))
}

## DATA GENERATION SCENARIOS
# Scenario 1: Gaussian w/ one strong segment
run_scen1 = function(n){
  l1 = floor(n*0.3)
  l2 = floor(n*0.2)
  l3 = floor(n*0.5)
  
  # bg estimates and segmentations
  lookback = floor(n*0.5)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 3, 1), rnorm(l3, 0, 1))
  pen = autoset_penalty(ts)
  out = apply_all_detectors(ts, lookback, pen)
  
  # attach cost using the set bg parameter and segments
  print("calculating likelihood")
  out$likTheor = getLik(data.frame("s"=l1+1, "e"=l1+l2), ts, 1, pen)
  return(out)
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
  out = apply_all_detectors(ts, lookback, pen)
  
  # using the set bg parameter and segments
  print("calculating likelihood 2")
  out$likTheor = getLik(data.frame("s"=c(l1, l1+l2+l3, l1+l2+l3+l4+l5)+1,
                               "e"=c(l1+l2, l1+l2+l3+l4, l1+l2+l3+l4+l5+l6)), ts, 1, pen)
  return(out)
}

# Scenario 3: t-distribution w/ one large segment
run_scen3 = function(n){
  l1 = floor(n*0.2)
  l2 = floor(n*0.4)
  l3 = floor(n*0.4)
  
  lookback = floor(n*0.5)
  ts = c(rt(l1, 3)+0, rt(l2, 3)+2, rt(l3, 3)+0)
  pen = autoset_penalty(ts)
  out = apply_all_detectors(ts, lookback, pen, sqrt(3))
  
  # using the set bg parameter and segments
  print("calculating likelihood 3")
  out$likTheor = getLik(data.frame("s"=l1+1, "e"=l1+l2), ts, sqrt(3), pen)
  return(out)
}

## RUN THE SIMULATIONS:
detsalln = data.frame()
allsegs = matrix(0, ncol=7)
ntotest = c(30, 60, 90, 130, 180, 240, 320, 440, 550, 750)
niter = 500

# at each n:
for(NPOINTS in ntotest){
  print(sprintf("Working on n=%d", NPOINTS))
  dets1 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, median=NA, theta.apelt=NA, 
                     scenario=1, fBG=0.8, likTheor=NA, likEst1=NA, likEst2=NA)
  dets2 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, median=NA, theta.apelt=NA,
                     scenario=2, fBG=0.7, likTheor=NA, likEst1=NA, likEst2=NA)
  dets3 = data.frame(n=NPOINTS, iter=1:niter, theta0=NA, median=NA, theta.apelt=NA,
                     scenario=3, fBG=0.6, likTheor=NA, likEst1=NA, likEst2=NA)
  for(i in 1:niter){
    print(i)
    res = run_scen1(NPOINTS)
    dets1[i, "theta0"] = res$wt
    dets1[i, "median"] = res$med
    dets1[i, "theta.apelt"] = res$theta.apelt
    dets1[i, "likTheor"] = res$likTheor
    dets1[i, "likEst1"] = res$likEst1
    dets1[i, "likEst2"] = res$likEst2
    
    allsegs = rbind(allsegs, cbind(res$segs1, NPOINTS, i, 1, 1))
    allsegs = rbind(allsegs, cbind(res$segs2, NPOINTS, i, 2, 1))
    allsegs = rbind(allsegs, cbind(res$segsan, NPOINTS, i, 3, 1))
    allsegs = rbind(allsegs, cbind(res$segsap, NPOINTS, i, 4, 1))
    print("--- done 1 ---")
    
    res = run_scen2(NPOINTS)
    dets2[i, "theta0"] = res$wt
    dets2[i, "median"] = res$med
    dets2[i, "theta.apelt"] = res$theta.apelt
    dets2[i, "likTheor"] = res$likTheor
    dets2[i, "likEst1"] = res$likEst1
    dets2[i, "likEst2"] = res$likEst2

    allsegs = rbind(allsegs, cbind(res$segs1, NPOINTS, i, 1, 2))
    allsegs = rbind(allsegs, cbind(res$segs2, NPOINTS, i, 2, 2))
    allsegs = rbind(allsegs, cbind(res$segsan, NPOINTS, i, 3, 2))
    allsegs = rbind(allsegs, cbind(res$segsap, NPOINTS, i, 4, 2))
    print("--- done 2 ---")
    
    res = run_scen3(NPOINTS)
    dets3[i, "theta0"] = res$wt
    dets3[i, "median"] = res$med
    dets3[i, "theta.apelt"] = res$theta.apelt
    dets3[i, "likTheor"] = res$likTheor
    dets3[i, "likEst1"] = res$likEst1
    dets3[i, "likEst2"] = res$likEst2
    
    allsegs = rbind(allsegs, cbind(res$segs1, NPOINTS, i, 1, 3))
    allsegs = rbind(allsegs, cbind(res$segs2, NPOINTS, i, 2, 3))
    allsegs = rbind(allsegs, cbind(res$segsan, NPOINTS, i, 3, 3))
    allsegs = rbind(allsegs, cbind(res$segsap, NPOINTS, i, 4, 3))
    print("--- done 3 ---")
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
write.table(detsalln, "../drafts/changepoint-method/results-sim-new/sim1-estimates.tsv", quote=F, sep="\t", col.names=T, row.names=F)
write.table(allsegs, "../drafts/changepoint-method/results-sim-new/sim1-detections.tsv", quote=F, sep="\t", col.names=T, row.names=F)
