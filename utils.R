rgmrf <- function(n, mean = 0, Q) {
  if(length(mean) == 1) mean <- rep(mean, nrow(Q))
  samples <- RTMB:::rgmrf0(n, Q)
  samples <- samples + mean
  t(samples)
}
