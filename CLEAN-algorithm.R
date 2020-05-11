## Changepoint detection, OUR ALGORITHM (unknown theta0)
## for change in mean of normal

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
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  
  # Main loop:
  for(t in 2:length(ts)){
    # cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0(ts[t], wt[t-1], SD)
    
    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    # (maxlookback limited if t short)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      # coord relative to t (1=at t):
      t2 = possibleKs[i]
      segcost[i] = bestcost[t2] + cost1(ts[(t2+1):t], SD)
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
      #wt[t] = mean(ts[bgpoints])
      wt[t] = wt[t-1] - (wt[t-1] - ts[t])*1/length(bgpoints)
    } else {
      bestcost[t] = bestsegcost
      # add a segment
      newseg = c(bestsegstart, t, mean(ts[bestsegstart:t]))
      segs[[t]] = rbind(segs[[bestsegstart-1]], newseg)
      wt[t] = wt[bestsegstart-1]
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
  output = list(segs = segs[[length(segs)]][-1,,drop=F], wt = wt)
  return(output)
}

# estimate chps with known theta0 (to get certain segments)
shortfixed <- function(ts, theta0, MAXLOOKBACK, PEN, SD){
  # Initialize:
  # bestcost[t] := F(all x[1:t])
  bestcost = c(0)
  # K: possible starts of segments to consider
  possibleKs = c(1)
  
  # Output:
  # for each t, a matrix of segment starts-ends up to t
  segs = vector(mode="list", length=length(ts))
  segs[[1]] = matrix(0, nrow=1, ncol=3)
  
  # Main loop:
  for(t in 2:length(ts)){
    # cat(sprintf("\nCycle %d\n", t))
    # Cost if t came from background:
    bgcost = bestcost[t-1] + cost0(ts[t], theta0, SD)
    
    # Cost if t came from segment:
    # over all possible segments from t:t to t-MLB+1:t
    # (seg length >= 1)
    # (maxlookback limited if t short)
    segcost = rep(Inf, length(possibleKs))
    for(i in seq_along(possibleKs)){
      # coord relative to t (1=at t):
      t2 = possibleKs[i]
      segcost[i] = bestcost[t2] + cost1(ts[(t2+1):t], SD)
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
  output = list(segs = segs[[length(segs)]][-1,,drop=F])
  return(output)
}
