options(stringsAsFactors = F)
source("MAIN-algorithm.R")
library(anomaly)
library(not)

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
  return(res)
}

# Our implementation of Jeng, Cai & Li (2010):
# Optimal Sparse Segment Identification With Application in Copy Number Variation Analysis.
segmenter_jeng = function(ts, theta0, MAXLOOKBACK, PEN, SD){
  # 0. original publication assumes standardized data (i.e. known parameters).
  # We assume these parameters are known and passed in.
  ts = (ts-theta0)/SD
  
  # 1. enumerate all segments w/ length <=MAXLOOKBACK, i.e. n x L elements
  segs = expand.grid(list("start"=1:length(ts), "length"=1:MAXLOOKBACK))
  segs$end = segs$start + segs$length - 1 
  # remove edge segments that would stick outside data length
  segs = segs[segs$end<=length(ts),]

  # 2. calculate the contrast for each segment (assumes standardized data)
  segs$contr = apply(segs, 1, function(x) abs(sum(ts[x[1]:x[3]])/sqrt(x[2])))
  
  # 3. filter the segments for contrast > "thr"=PEN
  segs = segs[segs$contr > PEN,]
  
  outsegs = data.frame("start"=0, "end"=0)
  
  cat("Starting loop...\n")
  while(nrow(segs)>0){
    # 4. find and store I = argmax_segment contrast
    bestseg = segs[which.max(segs$contr),c(1,3)]
    outsegs = rbind(outsegs, bestseg)
    # 5. filter segments to not overlap I
    segs = segs[segs$start>bestseg[1,2] | segs$end<bestseg[1,1], ]
    # 6. repeat 4-5 until no more segments
  }
  cat("loop finished\n")
  return(list(segs=as.matrix(outsegs)[-1,,drop=F]))
}

apply_all_detectors = function(ts, lookback, pen, capatype="mean"){
    if(lookback==0.33){
        # our method, no pruning
        res.full = fulldetector_prune(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=1, prune=0)
        if(nrow(res.full$segs)==0) res.full$segs = matrix(c(NA, NA, NA, NA), nrow=1)
    }
    # our method, full pruning
    res.full.pr = fulldetector_prune(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=1, prune=2)
    if(nrow(res.full.pr$segs)==0) res.full.pr$segs = matrix(c(NA, NA, NA, NA), nrow=1)
    
    # classical epidemic detector implemented in R:anomaly
    res.anom = capa.uv(ts, beta=pen, beta_tilde=pen, type=capatype, min_seg_len=2, max_seg_len=lookback, transform=identity)
    res.anom.c = collective_anomalies(res.anom)[,1:3]
    res.anom.p = point_anomalies(res.anom)
    # attach length 1 anomalies, if any
    if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
        res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location, mean.change=res.anom.p$strength) 
        res.anom.c = rbind(res.anom.c, res.anom.p)
    }
    # segment "type"
    if(nrow(res.anom.c)==0) res.anom.c = data.frame(start=NA, end=NA, mean.change=NA)
    res.anom.c$segtype = 1
    
    if(lookback==0.33){
        # Narrowest-over-threshold detector
        res.not = not(ts, method="not", contrast="pcwsConstMean")
        res.not = features(res.not)$cpt
        # adjust b/c it reports the last point of each segment as changepoint
        res.not = cbind(res.not+1, c(res.not[2:length(res.not)], length(ts)))
        res.not = cbind(res.not, apply(res.not, 1, function(x){
            if(any(is.na(x))) { NA }
            else { mean(ts[x[1]:x[2]]) }
        }))
        # segment "type"
        if(nrow(res.not)==0) res.not = data.frame(start=NA, end=NA, mean.change=NA)
        res.not = cbind(res.not, 1)
    }
    
    # aPELT_profile detector, our implementation
    res.apelt = aPELT_profile(ts, MAXLOOKBACK=lookback, PEN=pen, SD=1)$segs
    # segment "type"
    if(nrow(res.apelt)==0) res.apelt = matrix(c(NA, NA, NA), nrow=1)
    res.apelt = cbind(res.apelt, 1)
    
    # Sparse detector by Jeng, Cai & Li (2010)
    # We use their penalty of sqrt(2log(nL))? or our penalty of 3log(n)^1.1?
    pen_jeng = sqrt(2*log(length(ts)*lookback))
    res.sparse = segmenter_jeng(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen_jeng, SD=1)$segs
    # segment "type" and "effect"
    if(nrow(res.sparse)==0) res.sparse = matrix(c(NA, NA), nrow=1)
    res.sparse = cbind(res.sparse, NA, 1)
    
    if(lookback==0.33){
        return(list(segsf = res.full$segs, segsp = res.full.pr$segs,
                    segsan = as.matrix(res.anom.c), segsnot = res.not,
                    segsapelt = res.apelt, segssparse = res.sparse))
    } else {
        return(list(segsp = res.full.pr$segs, segsan = as.matrix(res.anom.c),
                    segsapelt = res.apelt, segssparse = res.sparse))
    }
}


## DATA GENERATION SCENARIOS
run_scen1_l2 = function(n, lfactor=0.33){
  l1 = floor(n*0.2)
  l2 = floor(n*0.1)
  l3 = floor(n*0.2)
  l4 = floor(n*0.2)
  l5 = floor(n*0.3)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 2, 1), rnorm(l3, 4, 1), rnorm(l4, 2, 1), rnorm(l5, 0, 1))
  
  lookback = floor(n*lfactor)
  pen = autoset_penalty(ts)
  
  res = apply_all_detectors(ts, lookback, pen)
  return(res)
}

run_scen1b_l2 = function(n){
    l1 = floor(n*0.2)
    l2 = floor(n*0.1)
    l3 = floor(n*0.2)
    l4 = floor(n*0.2)
    l5 = floor(n*0.3)
    ts = c(rnorm(l1, 0, 1), rnorm(l2, 0, 3), rnorm(l3, 2, 3), rnorm(l4, 0, 3), rnorm(l5, 0, 1))
    
    lookback = floor(n*0.33)
    pen = autoset_penalty(ts)
    
    res = apply_all_detectors(ts, lookback, pen, capatype="meanvar")
    return(res)
}

run_scen2_l2 = function(n){
  l1 = floor(n*0.2)
  l2 = floor(n*0.2) # N only
  l3 = floor(n*0.1)
  l4 = floor(n*0.1) # S only
  l5 = floor(n*0.1)
  l6 = floor(n*0.1) # S only
  l7 = floor(n*0.2)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 1.5, 1), rnorm(l3, 0, 1), rnorm(l4, 3, 1),
         rnorm(l5, 0, 1), rnorm(l6, -3, 1), rnorm(l7, 0, 1))
  
  lookback = floor(n*0.15)
  pen = autoset_penalty(ts)
  
  res = apply_all_detectors(ts, lookback, pen)
  return(res)
}

run_scen3_l2 = function(n){
  l1 = floor(n*0.1)
  l2 = floor(n*0.05) # S (to be repeated)
  l3 = ceiling(n*0.05) # background again (to be repeated)
  
  ts = c(rnorm(l1, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1),
         rnorm(l2, runif(1, -4, 4), 1), rnorm(l3, 0, 1))
  
  lookback = floor(n*0.2)
  pen = autoset_penalty(ts)
  
  res = apply_all_detectors(ts, lookback, pen)
  return(res)
}


## RUN THE SIMULATIONS:
allsegs2 = matrix(0, ncol=8)
ntotest2 = c(30, 60, 100, 150, 220)
niter2 = 500

# Args: N at which to run, and vector of iteration numbers
cycle_l2 = function(NPOINTS, iters, scen, lfactor=0.33){
  temp = matrix(0, ncol=8)
  print(sprintf("Working on n=%d", NPOINTS))
  for(i in iters){
    print(i)
    res = NULL
    if(scen==1){
      res = run_scen1_l2(NPOINTS, lfactor)
    } else if(scen==2){
      res = run_scen2_l2(NPOINTS)
    } else if(scen==3){
      res = run_scen3_l2(NPOINTS)
    } else if(scen==4){
      res = run_scen1b_l2(NPOINTS)
    }
    
    # These do not need to be done in the length tests
    if(lfactor==0.33){
        temp = rbind(temp, cbind(res$segsf, NPOINTS, i, 1, scen))
        temp = rbind(temp, cbind(res$segsnot, NPOINTS, i, 4, scen))
    }
    temp = rbind(temp, cbind(res$segsp, NPOINTS, i, 2, scen))
    temp = rbind(temp, cbind(res$segsan, NPOINTS, i, 3, scen))
    temp = rbind(temp, cbind(res$segsapelt, NPOINTS, i, 5, scen))
    temp = rbind(temp, cbind(res$segssparse, NPOINTS, i, 6, scen))
  }
  return(temp)
}


# at each n:
for(npoints in ntotest2){
  allsegs2 = rbind(allsegs2, cycle_l2(npoints, 1:niter2, 4))
  allsegs2 = rbind(allsegs2, cycle_l2(npoints, 1:niter2, 1))
  allsegs2 = rbind(allsegs2, cycle_l2(npoints, 1:niter2, 2))
  allsegs2 = rbind(allsegs2, cycle_l2(npoints, 1:niter2, 3))
}

nrow(allsegs2)
# (Recommend running this in a separate R session, RStudio doesn't manage to do garbage collection fast enough)

# allsegs had an empty row due to initialization
# also drop NAs which occur when an algorithm did not detect anything
head(allsegs2)
allsegs2 = allsegs2[which(allsegs2[,1]!=0),]
allsegs2bkp = allsegs2
allsegs2 = data.frame(allsegs2)
rownames(allsegs2) = NULL
colnames(allsegs2)[4] = "segtype"
colnames(allsegs2)[7] = "alg"
colnames(allsegs2)[8] = "scen"
allsegs2$segtype = ifelse(allsegs2$segtype==1, "seg", "nuis")
methodnames = c("1"="full", "2"="pruned", "3"="anomaly", "4"="NOT", "5"="aPELT", "6"="sparse")
allsegs2$alg = methodnames[as.character(allsegs2$alg)]


# save output:
write.table(allsegs2, "../drafts/changepoint-method/results-sim-new/sim2-detections.tsv", quote=F, sep="\t", col.names=T, row.names=F)
save("allsegs2bkp", file="../drafts/changepoint-method/results-sim-new/sim2-raw.RData")
# allsegs2 = read.table("../drafts/changepoint-method/results-sim/sim2-detections.tsv", h=T)


# ----------- LENGTH ROBUSTNESS TESTS --------------

# For producing ROC:
NPOINTS = 60
allsegs2_test = matrix(0, ncol=8)
for(lenf in c(0.25, 0.4)){
    tmp = cycle_l2(NPOINTS, 1:niter2, 1, lenf)
    tmp[,8] = lenf   # store the penalty factor instead of scen
    allsegs2_test = rbind(allsegs2_test, tmp)
}

nrow(allsegs2_test)
# (Recommend running this in a separate R session, RStudio doesn't manage to do garbage collection fast enough)

# allsegs had an empty row due to initialization
# also drop NAs which occur when an algorithm did not detect anything
head(allsegs2_test)
allsegs2_test = allsegs2_test[which(allsegs2_test[,1]!=0),]
allsegs2_testbkp = allsegs2_test
allsegs2_test = data.frame(allsegs2_test)
rownames(allsegs2_test) = NULL
colnames(allsegs2_test)[4] = "segtype"
colnames(allsegs2_test)[7] = "alg"
colnames(allsegs2_test)[8] = "lenfactor"
allsegs2_test$segtype = ifelse(allsegs2_test$segtype==1, "seg", "nuis")
methodnames = c("1"="full", "2"="pruned", "3"="anomaly", "4"="NOT", "5"="aPELT", "6"="sparse")
allsegs2_test$alg = methodnames[as.character(allsegs2_test$alg)]

# save output:
write.table(allsegs2_test, "../drafts/changepoint-method/results-sim-new/sim2-forlen.tsv", quote=F, sep="\t", col.names=T, row.names=F)
save("allsegs2_testbkp", file="../drafts/changepoint-method/results-sim-new/sim2-forlen.RData")