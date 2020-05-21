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
