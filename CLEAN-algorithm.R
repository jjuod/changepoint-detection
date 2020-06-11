## Changepoint detection, OUR ALGORITHM (unknown theta0)
## for change in mean of normal

# empty environment for caching
cache.env = new.env(parent=emptyenv())

# cost of points from Gaussian w/ given mean, sd
cost0sq <- function(xs, theta0, sd, log2pi){
  ssq = sum((xs-theta0)^2)
  length(xs)*log2pi + ssq / sd^2
}
# cost of points from Gaussian w/ free mean, sd
cost1sq <- function(xs, sd, log2pi){
  ssq = sum((xs-mean(xs))^2)
  length(xs)*log2pi + ssq / sd^2
}

# sx2 and sx: sums of x^2 and x corresponding to this seg
# n: seg length
cost1sq.rec <- function(sx2, sx, n, sd, log2pi){
  c1 = n*log2pi
  c2 = sx2 - sx^2/n
  c1+c2/sd^2
}

autoset_penalty <- function(data, alpha=3, delta=1.1){
  return(alpha*log(length(data)^delta))
}

# estimate theta0 (and putative segments)
# for a time series, using our SGD-style alg.
shortsgd <- function(ts, MAXLOOKBACK, PEN, SD){
  # Initialize:
  # sequence of estimates {w_t}
  wt = rep(ts[1], length(ts))
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # K: possible starts of segments to consider
  possibleKs = c(1)
  # sum(delta) for keeping track of background size easily
  bgsizes = rep(1, length(ts))
  
  # precompute:
  LOG2PI = log(2*pi*SD^2)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  
  # Main loop:
  for(t in 2:length(ts)){
    # cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0sq(ts[t], wt[t-1], SD, LOG2PI)
    
    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    # (maxlookback limited if t short)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      t2 = possibleKs[i]
      segcost[i] = bestcost[t2] + cost1sq(ts[(t2+1):t], SD, LOG2PI)
    }
    
    # proposed segment cost:    
    bestsegstart = which.min(segcost)
    bestsegcost = segcost[bestsegstart] + PEN
    # First x of the proposed segment:
    bestsegstart = possibleKs[bestsegstart]+1
    # cat(sprintf("Best segment was %d-%d: Fsplit = %.1f + %.1f + P
    #             vs F0 = %.1f + %.1f | theta0: %.2f\n",
    #             bestsegstart, t, bestcost[bestsegstart-1], cost1(ts[bestsegstart:t], SD),
    #             bestcost[t-1], cost0(ts[t], wt[t-1], SD), wt[t-1]))
    
    # Is background better than seg?
    # Fill out bestcost, segs, wt for this t
    if(bgcost < bestsegcost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
      # update number of est. bg points
      bgsizes[t] = bgsizes[t-1] + 1
      # update theta0
      wt[t] = wt[t-1] - (wt[t-1] - ts[t])/bgsizes[t]
    } else {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, mean(ts[bestsegstart:t]))
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
      # select bg size and parameter estimates as they were at k = argmin F(s)
      wt[t] = wt[bestsegstart-1]
      bgsizes[t] = bgsizes[bestsegstart-1]
    }
    
    # update and prune possible segment starts
    toprune = rep(F, length(possibleKs))
    # cat("\nsegcosts:\n")
    # print(segcost)
    for(i in seq_along(possibleKs)){
      if(bestcost[t] <= segcost[i]){
        toprune[i] = T
      }
    }
    # remove one possible segment start to limit lookback:
    if(t+1-possibleKs[1]>MAXLOOKBACK){
      toprune[1] = T
    }
    # print(toprune)
    possibleKs = possibleKs[!toprune]
    possibleKs = c(possibleKs, t)
    # cat(sprintf("After cycle %d, %d possible starts remain, from %d to %d",
                # t, length(possibleKs), min(possibleKs), max(possibleKs)))
    # print(possibleKs)
  }
  #cat(sprintf("Estimated theta0: %.2f\n", wt[t]))

  # form output object
  output = list(segs = segs[[length(segs)]][-1,,drop=F], wt = wt, cost = bestcost[length(bestcost)])
  return(output)
}

# estimate chps with known theta0 (to get certain segments)
shortfixed <- function(ts, theta0, MAXLOOKBACK, PEN, SD){
  # Initialize:
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # K: possible starts of segments to consider
  possibleKs = c(1)
  
  # precompute:
  LOG2PI = log(2*pi*SD^2)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  
  # Main loop:
  for(t in 2:length(ts)){
    # cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0sq(ts[t], theta0, SD, LOG2PI)
    
    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    # (maxlookback limited if t short)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      # coord relative to t (1=at t):
      t2 = possibleKs[i]
      segcost[i] = bestcost[t2] + cost1sq(ts[(t2+1):t], SD, LOG2PI)
    }
    
    # proposed segment cost:    
    bestsegstart = which.min(segcost)
    bestsegcost = segcost[bestsegstart] + PEN
    # First x of the proposed segment:
    bestsegstart = possibleKs[bestsegstart]+1
    # cat(sprintf("Best segment was %d-%d: Fsplit = %.1f + %.1f + P
    #             vs F0 = %.1f + %.1f \n",
    #             bestsegstart, t, bestcost[bestsegstart-1], cost1(ts[bestsegstart:t], SD),
    #             bestcost[t-1], cost0(ts[t], theta0, SD)))
    
    # Is background better than seg?
    # Fill out bestcost, segs, wt for this t
    if(bgcost < bestsegcost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
    } else {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, mean(ts[bestsegstart:t]))
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
    }
    
    # update and prune possible segment starts
    toprune = rep(F, length(possibleKs))
    # cat("\nsegcosts:\n")
    # print(segcost)
    for(i in seq_along(possibleKs)){
      if(bestcost[t] <= segcost[i]){
        toprune[i] = T
      }
    }
    # remove one possible segment start to limit lookback:
    if(t+1-possibleKs[1]>MAXLOOKBACK){
      toprune[1] = T
    }
    # print(toprune)
    possibleKs = possibleKs[!toprune]
    possibleKs = c(possibleKs, t)
    # cat(sprintf("After cycle %d, %d possible starts remain, from %d to %d",
    # t, length(possibleKs), min(possibleKs), max(possibleKs)))
    # print(possibleKs)
  }

  # form output object
  output = list(segs = segs[[length(segs)]][-1,,drop=F], cost = bestcost[length(bestcost)])
  return(output)
}

# Constructor for Alg 1 (unknown-background short segment detector)
# firstpoint: will be used as the initial guess of theta0 (w_1)
# datalen: max length of input data (for allocating memory)
# MAXLOOKBACK: max length of signal seg
# PEN: signal segment penalty
# SD: sigma of fN, fNS
detector <- function(firstpoint, datalen, MAXLOOKBACK, PEN, SD){
  # Initialize:
  # sequence of estimates {w_t}
  wt = rep(firstpoint, datalen)
  # bestcost[t] := F(all x[1:t])
  bestcost = rep(0, datalen)
  # K: possible starts of segments to consider
  possibleKs = c(1)
  # for easier tracking of number of est. bg points
  bgsizes = rep(1, datalen)
  
  # precompute:
  LOG2PI = log(2*pi*SD^2)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=datalen)
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  
  # Detector object:
  structure(list(ts=firstpoint, bgsizes=bgsizes, wt=wt,
                 SD=SD, LOG2PI=LOG2PI, PEN=PEN, MAXLOOKBACK=MAXLOOKBACK,
                 bestcost=bestcost, segs=segs, possibleKs=possibleKs),
            class="detector")
}

# Add a new point to the time series being analyzed
# Args:
# d: a detector object
# newpoint: x_t
# cache: using cached cost values (REQUIRES: cache.env)
# startadj: pointer such that startadj+1 value would be t=1 for this detector
detector.step <- function(d, newpoint, cache=T, startadj=NULL){
  d$ts = c(d$ts, newpoint)
  SD = d$SD
  LOG2PI = d$LOG2PI
  PEN = d$PEN
  MAXLOOKBACK = d$MAXLOOKBACK
  t = length(d$ts)
  # C^S_{i:t} will be tested for each i in this:
  possibleKs = d$possibleKs
  
  # cat(sprintf("\nCycle %d\n", t))
  # Cost if t came from background:
  bgcost = d$bestcost[t-1] + cost0sq(newpoint, d$wt[t-1], SD, LOG2PI)

  # Cost if t came from segment:
  # over all possible segments from t:t to t-MLB+1:t
  # (seg length >= 1)
  # (maxlookback limited if t short)
  segcost = rep(Inf, length(possibleKs))
  for(i in seq_along(possibleKs)){
    t2 = possibleKs[i]
    segcost[i] = i
    if(!cache){
      # recalculate cost every time
      segcost[i] = d$bestcost[t2] + cost1sq(d$ts[(t2+1):t], SD, LOG2PI)
    } else {
      cachedvalue = cache.env$cost1cache[t2+1+startadj, t+startadj]
      if(is.infinite(cachedvalue)){
        # first calculation (happens when this segment start was pruned from the main loop). cache
        # cat(sprintf("Cache miss for relative t=%d at matrix %d , %d\n", t2+1, t2+1+startadj, t+startadj))
        sx = cache.env$sxcache[t+startadj] - cache.env$sxcache[t2+startadj]
        sx2 = cache.env$sx2cache[t+startadj] - cache.env$sx2cache[t2+startadj]
        cache.env$cost1cache[t2+1+startadj, t+startadj] = cost1sq.rec(sx2, sx, t-t2, SD, LOG2PI)
        segcost[i] = d$bestcost[t2] + cache.env$cost1cache[t2+1+startadj, t+startadj]
      } else {
        # cat(sprintf("Looking up cached value for relative t=%d at matrix %d , %d\n", t2+1, t2+1+startadj, t+startadj))
        segcost[i] = d$bestcost[t2] + cachedvalue
      }
    }
  }
  
  # proposed segment cost:    
  bestsegstart = which.min(segcost)
  bestsegcost = segcost[bestsegstart] + PEN
  # First x of the proposed segment:
  bestsegstart = possibleKs[bestsegstart]+1
  # cat(sprintf("Best segment was %d-%d: Fsplit = %.1f + %.1f + P
  #             vs F0 = %.1f + %.1f | theta0: %.2f\n",
  #             bestsegstart, t, d$bestcost[bestsegstart-1], cost1sq(d$ts[bestsegstart:t], SD, LOG2PI),
  #             d$bestcost[t-1], cost0sq(d$ts[t], d$wt[t-1], SD, LOG2PI), d$wt[t-1]))
  
  # Is background better than seg?
  # Fill out bestcost, segs, wt for this t
  if(bgcost < bestsegcost){
    d$bestcost[t] = bgcost
    # no new changepoints
    d$segs[[t]] = d$segs[[t-1]]
    # update number of est. bg points
    d$bgsizes[t] = d$bgsizes[t-1] + 1
    # update theta0
    d$wt[t] = d$wt[t-1] - (d$wt[t-1] - newpoint)/d$bgsizes[t]
  } else {
    d$bestcost[t] = bestsegcost
    # add a segment
    newseg = c(bestsegstart, t, mean(d$ts[bestsegstart:t]))
    d$segs[[t]] = rbind(d$segs[[bestsegstart-1]], newseg)
    # select bg size and parameter estimates as they were at k = argmin F(s)
    d$wt[t] = d$wt[bestsegstart-1]
    d$bgsizes[t] = d$bgsizes[bestsegstart-1]
  }
  
  # update and prune possible segment starts
  toprune = rep(F, length(possibleKs))
  # cat("\nsegcosts:\n")
  # print(segcost)
  for(i in seq_along(possibleKs)){
    if(d$bestcost[t] <= segcost[i]){
      toprune[i] = T
    }
  }
  # remove one possible segment start to limit lookback:
  if(t+1-possibleKs[1]>MAXLOOKBACK){
    toprune[1] = T
  }
  # print(toprune)
  possibleKs = possibleKs[!toprune]
  d$possibleKs = c(possibleKs, t)
  # cat(sprintf("After cycle %d, %d possible starts remain, from %d to %d",
  # t, length(possibleKs), min(possibleKs), max(possibleKs)))
  # print(possibleKs)
  return(d)
}

# Add a set of points to the time series w/o segmenting
# Args:
# d: a detector object
# newpoint: {x_t}
# cache: using cached cost values (REQUIRES: cache.env)
# startadj: pointer such that startadj+1 value would be t=1 for this detector
detector.burnin <- function(d, newpoint, cache=T, startadj=NULL){
  prevt = length(d$ts)
  d$ts = c(d$ts, newpoint)
  SD = d$SD
  LOG2PI = d$LOG2PI
  PEN = d$PEN
  MAXLOOKBACK = d$MAXLOOKBACK
  t = length(d$ts)
  
  # cat(sprintf("\nCycle %d\n", t))
  # update theta0, assuming all points are from B
  d$wt[t] = mean(d$ts)
  
  # Cost if t came from background:
  bgcost = d$bestcost[prevt] + cost0sq(newpoint, d$wt[prevt], SD, LOG2PI)
  
  # Fill out bestcost, segs, wt for this t
  d$bestcost[t] = bgcost
  # no new changepoints
  d$segs[[t]] = d$segs[[prevt]]
  # update number of est. bg points
  d$bgsizes[t] = d$bgsizes[prevt] + length(newpoint)
  
  # Override any previous segment starts as the newpoints are assumed to be B:
  d$possibleKs = c(t)
  # cat(sprintf("After cycle %d, %d possible starts remain, from %d to %d",
  # t, length(possibleKs), min(possibleKs), max(possibleKs)))
  # print(possibleKs)
  return(d)
}

print.detector <- function(d){
  cat("Detector setup:\n")
  cat(sprintf("SD = %f / Max segment lenth = %i / Penalty = %f\n", d$SD, d$MAXLOOKBACK, d$PEN))
  cat(sprintf("Current data length: %i (of %i max allowed)\n", length(d$ts), length(d$wt)))
  cat("Current segments estimated:\n")
  print(d$segs[[length(d$segs)]])
  cat(sprintf("Last theta0 estimate: %f based on %i points\n", d$wt[length(d$ts)], d$bgsizes[length(d$ts)]))
  cat(sprintf("Last pruning set: %i elements between %i-%i\n", length(d$possibleKs), min(d$possibleKs), max(d$possibleKs)))
  cat(sprintf("Most recent estimates of F(n):\n"))
  print(d$bestcost[max(1, length(d$ts)-5):length(d$ts)])
}

# Full method for detecting in presence of nuisance segments
# Reference version, with some optimization (n^2.5)
# theta0: mean of fB
# SD: sigma of fB, fN, fS, fNS
# PEN: signal segment penalty
# PEN2: nuisance penalty
fulldetector_noprune_reference <- function(ts, theta0, MAXLOOKBACK, PEN, PEN2, SD){
  # Initialize:
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # possible starts of segments to consider
  possibleKs = c(1)
  # possible starts of nuisances to consider
  possibleKNs = c()
  
  # For storing the inner loop detectors starting at each t
  detectors = vector(mode="list", length=length(ts))
  
  # precompute:
  LOG2PI = log(2*pi*SD^2)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=4)
  
  # Main loop:
  for(t in 2:length(ts)){
    cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0sq(ts[t], theta0, SD, LOG2PI)
    
    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      t2 = possibleKs[i]
      segcost[i] = bestcost[t2] + cost1sq(ts[(t2+1):t], SD, LOG2PI)
    }
    
    # Cost if t came from nuisance (w/ or w/o segments, C'):
    # over all possible nuisances from 1:t to t-MLB:t
    # (length > maxlookback)
    
    nuiscost = rep(Inf, length(possibleKNs))
    for(i in seq_along(possibleKNs)){
      t2 = possibleKNs[i]

      if(is.null(detectors[[t2+1]])){
        # init the detector for new time point
        # cat(sprintf("\nCreating a detector at %d\n", t2+1))
        detectors[[t2+1]] = detector(ts[t2+1], length(ts)-t2, MAXLOOKBACK, PEN, SD)
        # step it up until current t
        # (could instead enforce that all this initial part is not checked for S later)
        for(t3 in (t2+2):t){
          detectors[[t2+1]] = detector.step(detectors[[t2+1]], ts[t3], cache=F)
        }
      } else {
        # or add one new point
        detectors[[t2+1]] = detector.step(detectors[[t2+1]], ts[t], cache=F)
      }
      
      # cat(sprintf("Estimated theta0: %.2f\n", wt[t]))

      # save other parameters from the inner loop
      lastpos = length(detectors[[t2+1]]$ts)
      # cat(sprintf("Inner loop cost: %.1f , nseg: %d \n", ressgd$cost, nrow(ressgd$segs)))
      nuiscost[i] = bestcost[t2] + detectors[[t2+1]]$bestcost[lastpos]
    }
    
    # proposed segment cost:    
    bestsegstart = which.min(segcost)
    bestsegcost = segcost[bestsegstart] + PEN
    # First x of the proposed segment:
    bestsegstart = possibleKs[bestsegstart]+1
    # cat(sprintf("Best segment was %d-%d: Fseg = %.1f + %.1f + P
    #             vs F0 = %.1f + %.1f \n",
    #             bestsegstart, t, bestcost[bestsegstart-1], cost1sq(ts[bestsegstart:t], SD, LOG2PI),
    #             bestcost[t-1], cost0sq(ts[t], theta0, SD, LOG2PI)))
    
    # proposed nuisance cost:
    if(length(possibleKNs)>0){
      bestnuisstart = which.min(nuiscost)   # KN index, not actual positions
      bestnuiscost = nuiscost[bestnuisstart] + PEN2
      # First x of the proposed nuisance:
      bestnuisstart = possibleKNs[bestnuisstart]+1
      
      # Other info of the proposed nuisance:
      bestdet = detectors[[bestnuisstart]]
      lastpos = length(bestdet$ts)
      bestnuistheta = bestdet$wt[lastpos]
      bestnuissegs = bestdet$segs[[lastpos]][-1,,drop=F]
      
      # (to avoid recalculating or storing it:)
      # actualFNS = bestnuiscost-PEN2-bestcost[bestnuisstart-1]
      # cat(sprintf("Best nuisance was %d-%d: Fnuis = %.1f + %.1f + P2
      #           vs F0 = %.1f + %.1f \n",
      #             bestnuisstart, t, bestcost[bestnuisstart-1], actualFNS,
      #             bestcost[t-1], cost0sq(ts[t], theta0, SD, LOG2PI)))
      # cat(sprintf("FB %.1f / FS %.1f / FN %.1f \n", bgcost, bestsegcost, bestnuiscost))
    } else {
      bestnuisstart = 0
      bestnuiscost = Inf
    }
    
    # Is background better than signal or nuisance?
    # Fill out bestcost, segs, wt for this t
    if(bgcost < bestsegcost & bgcost < bestnuiscost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
    } else if (bestsegcost < bestnuiscost) {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, mean(ts[bestsegstart:t]), 1)
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
    } else {
      bestcost[t] = bestnuiscost
      # add a nuisance
      newseg = c(bestnuisstart, t, bestnuistheta, 2)
      # add signals overlapping this nuisance
      if(nrow(bestnuissegs)>0){
        bestnuissegs = cbind(bestnuissegs, 1)
        # adjust start-end pos, b/c inner loop reports relative to its start:
        bestnuissegs[,1] = bestnuisstart + bestnuissegs[,1] - 1
        bestnuissegs[,2] = bestnuisstart + bestnuissegs[,2] - 1
        newseg = rbind(newseg, bestnuissegs)
      }
      segs[[t]] = rbind(segs[[bestnuisstart-1]], newseg)
    }
    
    # update and prune possible segment starts
    toprune = rep(F, length(possibleKs))
    # cat("\nsegcosts:\n")
    # print(segcost)
    for(i in seq_along(possibleKs)){
      if(bestcost[t] <= segcost[i]){
        toprune[i] = T
      }
    }
    # remove one possible segment start to limit lookback:
    if(t+1-possibleKs[1]>MAXLOOKBACK){
      toprune[1] = T
    }
    possibleKs = possibleKs[!toprune]
    possibleKs = c(possibleKs, t)
    
    cat(sprintf("After cycle %d, %d possible segment starts remain, from %d to %d\n",
                t, length(possibleKs), min(possibleKs), max(possibleKs)))
    # print(possibleKs)
    
    # for nuisance segment starts, just add one here - no pruning
    newNstart = t-MAXLOOKBACK
    if(newNstart > 0){
      possibleKNs = c(possibleKNs, newNstart)
    }
    cat(sprintf("After cycle %d, %d possible nuisance starts remain\n",
                t, length(possibleKNs)))
  }
  
  # form output object
  output = list(segs = segs[[length(segs)]][-1,,drop=F])
  return(output)
}

# Full method for detecting in presence of nuisance segments
# Optimized version
# theta0: mean of fB
# SD: sigma of fB, fN, fS, fNS
# PEN: signal segment penalty
# PEN2: nuisance penalty
fulldetector_noprune <- function(ts, theta0, MAXLOOKBACK, PEN, PEN2, SD){
  # Initialize:
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # possible starts of segments to consider
  possibleKs = c(1)
  # possible starts of nuisances to consider
  possibleKNs = c()
  
  # For storing the inner loop detectors starting at each t
  detectors = vector(mode="list", length=length(ts))
  # For caching seg costs and cumsums of x, x^2
  # (b/c the same seg might be tested in C0 and multiple C' cycles)
  # cost0cache?
  assign("cost1cache", matrix(Inf, nrow=length(ts), ncol=length(ts)), envir = cache.env)
  assign("sx2cache", cumsum(ts^2), envir = cache.env)
  assign("sxcache", cumsum(ts), envir = cache.env)
  
  # precompute:
  LOG2PI = log(2*pi*SD^2)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=4)
  
  PRUNETR_K = 0
  PRUNETR_KN = 0
  
  # Main loop:
  for(t in 2:length(ts)){
    cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0sq(ts[t], theta0, SD, LOG2PI)

    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      t2 = possibleKs[i]
      # This is the first computation with x_t, so cache it
      # \sum_{i=t2+1}^{t} x_i
      sx = cache.env$sxcache[t] - cache.env$sxcache[t2]
      # \sum_{i=t2+1}^{t} x^2_i
      sx2 = cache.env$sx2cache[t] - cache.env$sx2cache[t2]
      cache.env$cost1cache[t2+1, t] = cost1sq.rec(sx2, sx, t-t2, SD, LOG2PI)
      segcost[i] = bestcost[t2] + cache.env$cost1cache[t2+1, t]
    }

    # Cost if t came from nuisance (w/ or w/o segments, C'):
    # over all possible nuisances from 1:t to t-MLB:t
    # (length > maxlookback)

    nuiscost = rep(Inf, length(possibleKNs))
    for(i in seq_along(possibleKNs)){
      t2 = possibleKNs[i]

      if(is.null(detectors[[t2+1]])){
        # init the detector for new time point
        # cat(sprintf("\nCreating a detector at %d\n", t2+1))
        detectors[[t2+1]] = detector(ts[t2+1], length(ts)-t2, MAXLOOKBACK, PEN, SD)
        # step it up until current t
        # (could instead enforce that all this initial part is not checked for S later)
        for(t3 in (t2+2):t){
          detectors[[t2+1]] = detector.step(detectors[[t2+1]], ts[t3], cache=T, startadj=t2)
        }
      } else {
        # or add one new point
        detectors[[t2+1]] = detector.step(detectors[[t2+1]], ts[t], cache=T, startadj=t2)
      }

      # cat(sprintf("Estimated theta0: %.2f\n", wt[t]))

      # save other parameters from the inner loop
      lastpos = length(detectors[[t2+1]]$ts)
      # cat(sprintf("Inner loop cost: %.1f , nseg: %d \n", ressgd$cost, nrow(ressgd$segs)))
      nuiscost[i] = bestcost[t2] + detectors[[t2+1]]$bestcost[lastpos]
    }

    # proposed segment cost:
    bestsegstart = which.min(segcost)
    bestsegcost = segcost[bestsegstart] + PEN
    # First x of the proposed segment:
    bestsegstart = possibleKs[bestsegstart]+1
    # cat(sprintf("Best segment was %d-%d: Fseg = %.1f + %.1f + P
    #             vs F0 = %.1f + %.1f \n",
    #             bestsegstart, t, bestcost[bestsegstart-1], cost1sq(ts[bestsegstart:t], SD, LOG2PI),
    #             bestcost[t-1], cost0sq(ts[t], theta0, SD, LOG2PI)))

    # proposed nuisance cost:
    if(length(possibleKNs)>0){
      bestnuisstart = which.min(nuiscost)   # KN index, not actual positions
      bestnuiscost = nuiscost[bestnuisstart] + PEN2
      # First x of the proposed nuisance:
      bestnuisstart = possibleKNs[bestnuisstart]+1

      # Other info of the proposed nuisance:
      bestdet = detectors[[bestnuisstart]]
      lastpos = length(bestdet$ts)
      bestnuistheta = bestdet$wt[lastpos]
      bestnuissegs = bestdet$segs[[lastpos]][-1,,drop=F]

      # (to avoid recalculating or storing it:)
      # actualFNS = bestnuiscost-PEN2-bestcost[bestnuisstart-1]
      # cat(sprintf("Best nuisance was %d-%d: Fnuis = %.1f + %.1f + P2
      #           vs F0 = %.1f + %.1f \n",
      #             bestnuisstart, t, bestcost[bestnuisstart-1], actualFNS,
      #             bestcost[t-1], cost0sq(ts[t], theta0, SD, LOG2PI)))
      # cat(sprintf("FB %.1f / FS %.1f / FN %.1f \n", bgcost, bestsegcost, bestnuiscost))
    } else {
      bestnuisstart = 0
      bestnuiscost = Inf
    }

    # Is background better than signal or nuisance?
    # Fill out bestcost, segs, wt for this t
    if(bgcost < bestsegcost & bgcost < bestnuiscost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
    } else if (bestsegcost < bestnuiscost) {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, mean(ts[bestsegstart:t]), 1)
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
    } else {
      bestcost[t] = bestnuiscost
      # add a nuisance
      newseg = c(bestnuisstart, t, bestnuistheta, 2)
      # add signals overlapping this nuisance
      if(nrow(bestnuissegs)>0){
        bestnuissegs = cbind(bestnuissegs, 1)
        # adjust start-end pos, b/c inner loop reports relative to its start:
        bestnuissegs[,1] = bestnuisstart + bestnuissegs[,1] - 1
        bestnuissegs[,2] = bestnuisstart + bestnuissegs[,2] - 1
        newseg = rbind(newseg, bestnuissegs)
      }
      segs[[t]] = rbind(segs[[bestnuisstart-1]], newseg)
    }

    # update and prune possible segment starts
    toprune = rep(F, length(possibleKs))
    # cat("\nsegcosts:\n")
    # print(segcost)
    for(i in seq_along(possibleKs)){
      if(bestcost[t] <= segcost[i]){
        toprune[i] = T
      }
    }
    # remove one possible segment start to limit lookback:
    if(t+1-possibleKs[1]>MAXLOOKBACK){
      toprune[1] = T
    }
    possibleKs = possibleKs[!toprune]
    possibleKs = c(possibleKs, t)
    PRUNETR_K = PRUNETR_K + length(possibleKs)

    # cat(sprintf("After cycle %d, %d possible segment starts remain, from %d to %d\n",
    #             t, length(possibleKs), min(possibleKs), max(possibleKs)))
    # # print(possibleKs)

    # for nuisance segment starts, just add one here - no pruning
    newNstart = t-MAXLOOKBACK
    if(newNstart > 0){
      possibleKNs = c(possibleKNs, newNstart)
    }
    PRUNETR_KN = PRUNETR_KN + length(possibleKNs)
    
    # cat(sprintf("After cycle %d, %d possible nuisance starts remain\n",
    #             t, length(possibleKNs)))
  }
  
  cat(sprintf("total segments checked: %d (%.2f on average)\n", PRUNETR_K, PRUNETR_K/length(ts)))
  cat(sprintf("total nuisances checked: %d (%.2f on average)\n", PRUNETR_KN, PRUNETR_KN/length(ts)))
  # form output object
  rm("cost1cache", envir=cache.env)
  rm("sx2cache", envir=cache.env)
  rm("sxcache", envir=cache.env)
  output = list(segs = segs[[length(segs)]][-1,,drop=F])
  return(output)
}

# Full method for detecting in presence of nuisance segments
# Optimzed with pruning
# theta0: mean of fB
# SD: sigma of fB, fN, fS, fNS
# PEN: signal segment penalty
# PEN2: nuisance penalty
# BURNIN: number of burnin iterations (points assumed to be from fN at the start of each nuisance).
# prune: 0/1/2 = none/partial/full pruning of nuisance starts
fulldetector_prune <- function(ts, theta0, MAXLOOKBACK, PEN, PEN2, SD, BURNIN, prune){
  # Initialize:
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # possible starts of segments to consider
  possibleKs = c(1)
  # possible starts of nuisances to consider
  possibleKNs = c()
  # could tweak the inner loop to avoid this requirement, but for now:
  BURNIN = max(2, min(BURNIN, MAXLOOKBACK))
  
  # For storing the inner loop detectors starting at each t
  detectors = vector(mode="list", length=length(ts))
  # For caching seg costs and cumsums of x, x^2
  # (b/c the same seg might be tested in C0 and multiple C' cycles)
  # cost0cache?
  assign("cost1cache", matrix(Inf, nrow=length(ts), ncol=length(ts)), envir = cache.env)
  assign("sx2cache", cumsum(ts^2), envir = cache.env)
  assign("sxcache", cumsum(ts), envir = cache.env)
  
  # precompute:
  LOG2PI = log(2*pi*SD^2)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=4)
  
  PRUNETR_K = 0
  PRUNETR_KN = 0
  
  # Main loop:
  for(t in 2:length(ts)){
    cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0sq(ts[t], theta0, SD, LOG2PI)
    
    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      t2 = possibleKs[i]
      # This is the first computation with x_t, so cache it
      # \sum_{i=t2+1}^{t} x_i
      sx = cache.env$sxcache[t] - cache.env$sxcache[t2]
      # \sum_{i=t2+1}^{t} x^2_i
      sx2 = cache.env$sx2cache[t] - cache.env$sx2cache[t2]
      cache.env$cost1cache[t2+1, t] = cost1sq.rec(sx2, sx, t-t2, SD, LOG2PI)
      segcost[i] = bestcost[t2] + cache.env$cost1cache[t2+1, t]
    }
    
    # Cost if t came from nuisance (w/ or w/o segments, C'):
    # over all possible nuisances from 1:t to t-MLB:t
    # (length > maxlookback)
    # BURNIN=2
    nuiscost = rep(Inf, length(possibleKNs))
    for(i in seq_along(possibleKNs)){
      t2 = possibleKNs[i]
      
      if(is.null(detectors[[t2+1]])){
        # init the detector for new time point
        # cat(sprintf("\nCreating a detector at %d\n", t2+1))
        det = detector(ts[t2+1], length(ts)-t2, MAXLOOKBACK, PEN, SD)
        # step it up until current t
        # det = detector.burnin(det, ts[(t2+2):(t2+BURNIN)], cache=T, startadj=t2)
        # for(t3 in (t2+BURNIN+1):t){
        for(t3 in (t2+2):t){
          det = detector.step(det, ts[t3], cache=T, startadj=t2)
        }
        detectors[[t2+1]] = det
      } else {
        # or add one new point
        detectors[[t2+1]] = detector.step(detectors[[t2+1]], ts[t], cache=T, startadj=t2)
      }
      
      # cat(sprintf("Estimated theta0: %.2f\n", wt[t]))
      
      # save other parameters from the inner loop
      lastpos = length(detectors[[t2+1]]$ts)
      # cat(sprintf("Inner loop cost: %.1f , nseg: %d \n", ressgd$cost, nrow(ressgd$segs)))
      nuiscost[i] = bestcost[t2] + detectors[[t2+1]]$bestcost[lastpos]
    }
    
    # proposed segment cost:
    bestsegstart = which.min(segcost)
    bestsegcost = segcost[bestsegstart] + PEN
    # First x of the proposed segment:
    bestsegstart = possibleKs[bestsegstart]+1
    # cat(sprintf("Best segment was %d-%d: Fseg = %.1f + %.1f + P
    #             vs F0 = %.1f + %.1f \n",
    #             bestsegstart, t, bestcost[bestsegstart-1], cost1sq(ts[bestsegstart:t], SD, LOG2PI),
    #             bestcost[t-1], cost0sq(ts[t], theta0, SD, LOG2PI)))
    
    # proposed nuisance cost:
    if(length(possibleKNs)>0){
      bestnuisstart = which.min(nuiscost)   # KN index, not actual positions
      bestnuiscost = nuiscost[bestnuisstart] + PEN2
      # First x of the proposed nuisance:
      bestnuisstart = possibleKNs[bestnuisstart]+1
      
      # Other info of the proposed nuisance:
      bestdet = detectors[[bestnuisstart]]
      lastpos = length(bestdet$ts)
      bestnuistheta = bestdet$wt[lastpos]
      bestnuissegs = bestdet$segs[[lastpos]][-1,,drop=F]
      
      # (to avoid recalculating or storing it:)
      # actualFNS = bestnuiscost-PEN2-bestcost[bestnuisstart-1]
      # cat(sprintf("Best nuisance was %d-%d: Fnuis = %.1f + %.1f + P2
      #           vs F0 = %.1f + %.1f \n",
      #             bestnuisstart, t, bestcost[bestnuisstart-1], actualFNS,
      #             bestcost[t-1], cost0sq(ts[t], theta0, SD, LOG2PI)))
      # cat(sprintf("FB %.1f / FS %.1f / FN %.1f \n", bgcost, bestsegcost, bestnuiscost))
    } else {
      bestnuisstart = 0
      bestnuiscost = Inf
    }
    
    # Is background better than signal or nuisance?
    # Fill out bestcost, segs, wt for this t
    if(bgcost < bestsegcost & bgcost < bestnuiscost){
      bestcost[t] = bgcost
      # no new changepoints
      segs[[t]] = segs[[t-1]]
    } else if (bestsegcost < bestnuiscost) {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, mean(ts[bestsegstart:t]), 1)
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
    } else {
      bestcost[t] = bestnuiscost
      # add a nuisance
      newseg = c(bestnuisstart, t, bestnuistheta, 2)
      # add signals overlapping this nuisance
      if(nrow(bestnuissegs)>0){
        bestnuissegs = cbind(bestnuissegs, 1)
        # adjust start-end pos, b/c inner loop reports relative to its start:
        bestnuissegs[,1] = bestnuisstart + bestnuissegs[,1] - 1
        bestnuissegs[,2] = bestnuisstart + bestnuissegs[,2] - 1
        newseg = rbind(newseg, bestnuissegs)
      }
      segs[[t]] = rbind(segs[[bestnuisstart-1]], newseg)
    }
    
    # update and prune possible segment starts
    toprune = rep(F, length(possibleKs))
    # cat("\nsegcosts:\n")
    # print(segcost)
    for(i in seq_along(possibleKs)){
      if(bestcost[t] <= segcost[i]){
        toprune[i] = T
      }
    }
    # remove one possible segment start to limit lookback:
    if(t+1-possibleKs[1]>MAXLOOKBACK){
      toprune[1] = T
    }
    possibleKs = possibleKs[!toprune]
    possibleKs = c(possibleKs, t)
    PRUNETR_K = PRUNETR_K + length(possibleKs)
    
    # cat(sprintf("After cycle %d, %d possible segment starts remain, from %d to %d\n",
    #             t, length(possibleKs), min(possibleKs), max(possibleKs)))
    # # print(possibleKs)
    
    # update and prune possible nuisance starts
    if(prune>0){
      toprune = rep(F, length(possibleKNs))
      # cat("\n nuiscosts:\n")
      # print(nuiscost)
      if(prune==2){
        for(i in seq_along(possibleKNs)){
          if(bestcost[t] <= nuiscost[i]){
            toprune[i] = T
          }
        }
      }
      # TO DO: prune=1 option here
      
      possibleKNs = possibleKNs[!toprune]
    }
    # once we reach sufficient length, add one each cycle
    newNstart = t-MAXLOOKBACK
    if(newNstart > 0){
      possibleKNs = c(possibleKNs, newNstart)
    }
    PRUNETR_KN = PRUNETR_KN + length(possibleKNs)
    
    cat(sprintf("After cycle %d, %d possible nuisance starts remain\n",
                t, length(possibleKNs)))
  }
  
  cat(sprintf("total segments checked: %d (%.2f on average)\n", PRUNETR_K, PRUNETR_K/length(ts)))
  cat(sprintf("total nuisances checked: %d (%.2f on average)\n", PRUNETR_KN, PRUNETR_KN/length(ts)))
  # form output object
  rm("cost1cache", envir=cache.env)
  rm("sx2cache", envir=cache.env)
  rm("sxcache", envir=cache.env)
  output = list(segs = segs[[length(segs)]][-1,,drop=F], bestcost)
  return(output)
}
