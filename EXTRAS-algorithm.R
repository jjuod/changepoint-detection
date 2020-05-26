# Here we keep various parts of our algorithm that are not part
# of the default method, but may be useful in customizing it.


# Cost function templates for easily replacing the likelihood with any d* available in R:
# cost of points from Gaussian w/ given mean, sd
cost0 <- function(xs, theta0, sd){
  liks = dnorm(xs, theta0, sd)
  -2 * sum(log(liks))
}
# cost of points from Gaussian w/ free mean, sd
cost1 <- function(xs, sd){
  liks = dnorm(xs, mean(xs), sd)
  -2 * sum(log(liks))
}

# Full method for detecting in presence of nuisance segments
# Original version, no optimization at all (n^4 basically)
# theta0: mean of fB
# MAXLOOKBACK: max length of signal seg
# SD: sigma of fB, fN, fS, fNS
# PEN: signal segment penalty
# PEN2: nuisance penalty
fulldetector_noprune_noopt <- function(ts, theta0, MAXLOOKBACK, PEN, PEN2, SD){
  # Initialize:
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # possible starts of segments to consider
  possibleKs = c(1)
  # possible starts of nuisances to consider
  possibleKNs = c()
  
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
    nuissegs = vector(mode="list", length=length(possibleKNs))
    nuistheta = rep(-99, length(possibleKNs))
    for(i in seq_along(possibleKNs)){
      t2 = possibleKNs[i]
      # save other parameters from the inner loop
      ressgd = shortsgd(ts[(t2+1):t], MAXLOOKBACK, PEN, SD)
      
      nuissegs[[i]] = ressgd$segs
      nuistheta[i] = ressgd$wt[length(ressgd$wt)]
      # cat(sprintf("Inner loop cost: %.1f , nseg: %d \n", ressgd$cost, nrow(ressgd$segs)))
      nuiscost[i] = bestcost[t2] + ressgd$cost
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
      bestnuissegs = nuissegs[[bestnuisstart]]
      bestnuistheta = nuistheta[bestnuisstart]
      # First x of the proposed segment:
      bestnuisstart = possibleKNs[bestnuisstart]+1
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
    if(t > MAXLOOKBACK){
      possibleKNs = c(possibleKNs, t-MAXLOOKBACK) 
    }
    cat(sprintf("After cycle %d, %d possible nuisance starts remain\n",
                t, length(possibleKNs)))
  }
  
  # form output object
  output = list(segs = segs[[length(segs)]][-1,,drop=F])
  return(output)
}
