options(stringsAsFactors = F)
source("CLEAN-algorithm.R")
library(anomaly)
library(not)

## DATA GENERATION SCENARIOS
run_scen1_l2 = function(n){
  l1 = floor(n*0.2)
  l2 = floor(n*0.1)
  l3 = floor(n*0.2)
  l4 = floor(n*0.2)
  l5 = floor(n*0.3)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 2, 1), rnorm(l3, 4, 1), rnorm(l4, 2, 1), rnorm(l5, 0, 1))
  
  lookback = floor(n*0.33)
  pen = autoset_penalty(ts)
  
  # our method, no pruning
  res.full = fulldetector_prune(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=1, prune=0)
  # our method, full pruning
  res.full.pr = fulldetector_prune(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=1, prune=2)
  # classical epidemic detector implemented in R:anomaly
  res.anom = capa.uv(ts, beta=pen, beta_tilde=pen, type="mean", min_seg_len=2, max_seg_len=lookback, transform=identity)
  res.anom.c = collective_anomalies(res.anom)[,1:3]
  res.anom.p = point_anomalies(res.anom)
  # attach length 1 anomalies, if any
  if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
    res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location, mean.change=res.anom.p$strength) 
    res.anom.c = rbind(res.anom.c, res.anom.p)
  }
  # segment "type"
  res.anom.c$segtype = 1
  
  # Narrowest-over-threshold detector
  res.not = not(ts, method="not", contrast="pcwsConstMean")
  res.not = features(res.not)$cpt
  res.not = cbind(res.not, c(res.not[2:length(res.not)]-1, length(ts)))
  res.not = cbind(res.not, apply(res.not, 1, function(x){
    if(any(is.na(x))) { NA }
    else { mean(ts[x[1]:x[2]]) }
  }))
  # segment "type"
  res.not = cbind(res.not, 1)
  
  return(list(segsf = res.full$segs, segsp = res.full.pr$segs, segsan = as.matrix(res.anom.c), segsnot = res.not))
}

run_scen2_l2 = function(n){
  l1 = floor(n*0.2)
  l2 = floor(n*0.2) # N only
  l3 = floor(n*0.1)
  l4 = floor(n*0.1) # S only
  l5 = floor(n*0.1)
  l6 = floor(n*0.1) # S only
  l7 = floor(n*0.2)
  ts = c(rnorm(l1, 0, 1), rnorm(l2, 1, 1), rnorm(l3, 0, 1), rnorm(l4, 3, 1),
         rnorm(l5, 0, 1), rnorm(l6, -3, 1), rnorm(l7, 0, 1))
  
  lookback = floor(n*0.15)
  pen = autoset_penalty(ts)
  
  # our method, no pruning
  res.full = fulldetector_prune(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=1, prune=0)
  # our method, full pruning
  res.full.pr = fulldetector_prune(ts, theta0=0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=1, prune=2)
  # classical epidemic detector implemented in R:anomaly
  res.anom = capa.uv(ts, beta=pen, beta_tilde=pen, type="mean", min_seg_len=2, max_seg_len=lookback, transform=identity)
  res.anom.c = collective_anomalies(res.anom)[,1:3]
  res.anom.p = point_anomalies(res.anom)
  # attach length 1 anomalies, if any
  if(nrow(res.anom.p)>1 | !is.na(res.anom.p$location[1])){
    res.anom.p = data.frame(start=res.anom.p$location, end=res.anom.p$location, mean.change=res.anom.p$strength) 
    res.anom.c = rbind(res.anom.c, res.anom.p)
  }
  # segment "type"
  res.anom.c$segtype = 1

  # Narrowest-over-threshold detector
  res.not = not(ts, method="not", contrast="pcwsConstMean")
  res.not = features(res.not)$cpt
  res.not = cbind(res.not, c(res.not[2:length(res.not)]-1, length(ts)))
  res.not = cbind(res.not, apply(res.not, 1, function(x){
    if(any(is.na(x))) { NA }
    else { mean(ts[x[1]:x[2]]) }
  }))
  # segment "type"
  res.not = cbind(res.not, 1)
  
  return(list(segsf = res.full$segs, segsp = res.full.pr$segs, segsan = as.matrix(res.anom.c), segsnot = res.not))
}


## RUN THE SIMULATIONS:
allsegs2 = matrix(0, ncol=8)
ntotest2 = c(30, 60, 100, 160, 240)
niter2 = 500

# Args: N at which to run, and vector of iteration numbers
cycle_l2 = function(NPOINTS, iters){
  temp = matrix(0, ncol=8)
  print(sprintf("Working on n=%d", NPOINTS))
  for(i in iters){
    print(i)
    res = run_scen1_l2(NPOINTS)
    
    temp = rbind(temp, cbind(res$segsf, NPOINTS, i, 1, 1))
    temp = rbind(temp, cbind(res$segsp, NPOINTS, i, 2, 1))
    temp = rbind(temp, cbind(res$segsan, NPOINTS, i, 3, 1))
    temp = rbind(temp, cbind(res$segsnot, NPOINTS, i, 4, 1))
  }
  return(temp)
}

cycle2_l2 = function(NPOINTS, iters){
  temp = matrix(0, ncol=8)
  print(sprintf("Working on n=%d", NPOINTS))
  for(i in iters){
    print(i)
    res = run_scen2_l2(NPOINTS)
    
    temp = rbind(temp, cbind(res$segsf, NPOINTS, i, 1, 2))
    temp = rbind(temp, cbind(res$segsp, NPOINTS, i, 2, 2))
    temp = rbind(temp, cbind(res$segsan, NPOINTS, i, 3, 2))
    temp = rbind(temp, cbind(res$segsnot, NPOINTS, i, 4, 2))
  }
  return(temp)
}

# at each n:
allsegs2 = rbind(allsegs2, cycle_l2(30, 1:niter2))
allsegs2 = rbind(allsegs2, cycle_l2(60, 1:niter2))
allsegs2 = rbind(allsegs2, cycle_l2(100, 1:niter2))
allsegs2 = rbind(allsegs2, cycle_l2(160, 1:niter2))
allsegs2 = rbind(allsegs2, cycle_l2(240, 1:niter2))

allsegs2 = rbind(allsegs2, cycle2_l2(30, 1:niter2))
allsegs2 = rbind(allsegs2, cycle2_l2(60, 1:niter2))
allsegs2 = rbind(allsegs2, cycle2_l2(100, 1:niter2))
allsegs2 = rbind(allsegs2, cycle2_l2(160, 1:niter2))
allsegs2 = rbind(allsegs2, cycle2_l2(240, 1:niter2))

nrow(allsegs2)
# (Recommend running this in a separate R session, RStudio doesn't manage to do garbage collection fast enough)

# allsegs had an empty row due to initialization
head(allsegs2)
allsegs2 = allsegs2[allsegs2[,1]!=0,]
allsegs2bkp = allsegs2
allsegs2 = data.frame(allsegs2)
rownames(allsegs2) = NULL
colnames(allsegs2)[4] = "segtype"
colnames(allsegs2)[7] = "alg"
colnames(allsegs2)[8] = "scen"
allsegs2$segtype = ifelse(allsegs2$segtype==1, "seg", "nuis")
methodnames = c("1"="full", "2"="pruned", "3"="anomaly", "4"="NOT")
allsegs2$alg = methodnames[as.character(allsegs2$alg)]

# save output:
write.table(allsegs2, "../drafts/changepoint-method/results-sim/sim2-detections.tsv", quote=F, sep="\t", col.names=T, row.names=F)
save("allsegs2bkp", file="../drafts/changepoint-method/results-sim/sim2-raw.RData")
# allsegs2 = read.table("../drafts/changepoint-method/results-sim/sim2-detections.tsv", h=T)
