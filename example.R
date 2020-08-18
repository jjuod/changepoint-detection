# Example usage of epidemic changepoint detectors (Algorithms 1 and 2)
# as described in Juodakis J. and Marsland S.:
# Epidemic changepoint detection in the presence of nuisance changes

# Load the detector functions
source("MAIN-algorithm.R")

# Generate some data
ts = c(rnorm(30, 0, 1),
       rnorm(25, 2, 1),
       rnorm(15, 0, 1),
       rnorm(10, -2, 1),
       rnorm(20, 0, 1))
plot(ts)

# Set penalty to default (3 * log(n) ^ 1.1)
pen = autoset_penalty(ts)

# Set maximum segment length
lookback = 30

# Provide an estimate of the background sigma
SD = 1

# Provide an estimate of the global background mean (Algorithm 2 only)
mu0 = median(ts) 


# Apply Algorithm 1 (online version):
res = shortsgd(ts, lookback, pen, SD)
# Final detections
res$segs
# Final background estimate
res$wt[length(res$wt)]


# Apply Algorithm 2 (with global pruning):
res2 = fulldetector_prune(ts, theta0=mu0, MAXLOOKBACK=lookback, PEN=pen, PEN2=pen, SD=SD, prune=2)
# Final detections
res2$segs

# With nicer formatting
res2df = data.frame(res2$segs)
colnames(res2df) = c("start", "end", "theta", "segtype")
rownames(res2df) = NULL
res2df$segtype = ifelse(res2df$segtype==1, "S", "N")
res2df
