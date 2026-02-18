library(LaMa)       # for HMM functions
library(RTMBdist)   # for ExGaussian distribution
library(fmesher)    # for mesh and FEM matrices
library(Matrix)     # for sparse matrices
# source("./functions/forward_banded.R") # sourcing banded forward algorithm (not in LaMa yet)
# source("./functions/hmm_corr.R")       # sourcing functions for HMM state decoding

### colors for plotting
color = c("#00000050", "red", "orange")


### loading data
data <- read.csv("./data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[,c("TIME","PDCSAP_FLUX")]
colnames(data) = c("time", "y")

data <- data[-1, ] # exclude first observation because NA

# interpolate time
time_na <- is.na(data$time)
data$time[time_na] <- approx(which(!time_na), data$time[!time_na], xout = which(time_na))$y

end_before_gap <- 9524
start_after_gap <- 10341

data <- data[-((end_before_gap + 1):(start_after_gap-1)), ]
data$trackID <- 2
data$trackID[1:end_before_gap] <- 1


# subset data
# data <- data[2:9524,] # first observation is NA, after 9524, long series of missing
# na_rle <- rle(is.na(data$y))
# cbind(na_rle$values, na_rle$lengths)[na_rle$values == 1, ]

# linearly interpolating missing observations
# data$y[is.na(data$y)] <- approx(data$time, data$y, data$time[is.na(data$y)])$y

# centering data
data$y <- scale(data$y, scale = FALSE)


# data <- data[1:4000,]

### creating 2meshs and finite element matrices
mesh1 <- fm_mesh_1d(data$time[data$trackID == 1])
spde1 <- fm_fem(mesh1)

mesh2 <- fm_mesh_1d(data$time[data$trackID == 2])
spde2 <- fm_fem(mesh2)


# ### Simple analysis without HMM
# nll <- function(par) {
#   getAll(par, dat)
#
#   ## observation model
#   # parameter transformations
#   sigma = exp(log_sigma); REPORT(sigma)
#   f <- f0 + mu; REPORT(f0); REPORT(f) # f is the smooth function
#   nll <- -sum(dnorm(y, f, sigma, log = TRUE)) # observation liklihood
#
#   ## GP model
#   # parameter transformations
#   tau <- exp(log_tau); REPORT(tau)
#   kappa <- exp(log_kappa); REPORT(kappa)
#   rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
#   omega <- plogis(logit_omega); REPORT(omega)
#   Q <- tau^2 * (kappa^4 * c0 + 2 * cos(pi*omega) * kappa^2 * g1 + g2); REPORT(Q)
#   nll <- nll - dgmrf(f - mu, 0, Q, log = TRUE) # GP likelihood
#
#   nll
# }
#
# # initial parameter list
# par <- list(
#   log_sigma = log(7),
#   log_tau = log(0.005),
#   log_kappa = log(30),
#   f0 = numeric(nrow(spde$c0)),
#   logit_omega = qlogis(0.9),
#   mu = -0.5
# )
#
# # data list
# dat <- list(
#   y = data$y,
#   c0 = spde$c0, g1 = spde$g1, g2 = spde$g2
# )
#
# obj_simple <- MakeADFun(nll, par, random = "f0")
# opt_simple <- nlminb(obj_simple$par, obj_simple$fn, obj_simple$gr)
# mod_simple <- obj_simple$report()
# mod_simple$tau
# mod_simple$kappa
# mod_simple$omega
#
# ### visualising results
# idx = 1:nrow(data)
# plot(data$time[idx], data$y[idx], pch = 16,
#      xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Stellar flare detection")
# lines(data$time[idx], mod_simple$f[idx], lwd = 2, col = "plum")
# # smooth is not very periodic
# # erratic behaviour whenever there is a flare


### HMM analysis
# HMM likelihood function
jnll <- function(par) {
  getAll(par, dat)

  #### state process ####
  # restricted tpm
  Gamma <- diag(3)
  Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] <- exp(eta)
  Gamma <- Gamma / rowSums(Gamma)
  # estimated initial distributions
  delta1 <- c(1, exp(logit_delta1))
  delta1 <- delta1 / sum(delta1)
  delta2 <- c(1, exp(logit_delta2))
  delta2 <- delta2 / sum(delta2)
  Delta <- rbind(delta1, delta2)

  #### state-dependent process ####
  # parameter transformations
  sigma <- exp(log_sigma); REPORT(sigma)
  r <- plogis(logit_r); REPORT(r)
  lambda <- exp(log_lambda); REPORT(lambda)
  # state-dependent densities
  f1 <- w1 + mu; REPORT(f1)
  f2 <- w2 + mu; REPORT(f2)
  f <- c(f1, f2); REPORT(f) # total smooth
  z <- y - f; REPORT(mu)

  # latent zs
  z[is.na(y)] <- z.star; REPORT(z)
  y.star <- z.star + f[is.na(y)]; REPORT(y.star)

  n <- length(z)
  # indices of the two separate time series
  idx <- list(2:(trackInd[2]-1), (trackInd[2]+1):n)

  lallprobs <- matrix(0, n, 3)
  for(ind in idx) {
    for(t in ind) {
      if(!is.na(z[t]) & !is.na(z[t-1])) {
        # regular measurement error
        lallprobs[t,1] <- dnorm(z[t], 0, sigma, log = TRUE)
        # firing
        lallprobs[t,2] <- dexgauss(z[t], z[t-1], sigma, lambda, log = TRUE)
        # decaying
        lallprobs[t,3] <- dnorm(z[t], r * z[t-1], sigma, log = TRUE)
      }
    }
  }

  # banded forward algorithm - HMM likelihood
  nll <- -forward(Delta, Gamma, lallprobs,
                  trackID = trackID, bw = bw, logspace = TRUE)

  ### GP part ###
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  rho <- sqrt(8) / kappa; REPORT(rho) # distance where corr has dropped to 0.1
  # omega <- plogis(logit_omega); REPORT(omega)
  cos_pi_omega <- 2 * plogis(u) - 1
  omega <- acos(cos_pi_omega) / pi; REPORT(omega)

  Q1 <- tau_sq * (kappa_sq*kappa_sq * spde1$c0 + 2 * cos_pi_omega * kappa_sq * spde1$g1 + spde1$g2)
  Q2 <- tau_sq * (kappa_sq*kappa_sq * spde2$c0 + 2 * cos_pi_omega * kappa_sq * spde2$g1 + spde2$g2)

  nll <- nll - dgmrf(w1, 0, Q1, log = TRUE) # GP likelihood 1
  nll <- nll - dgmrf(w2, 0, Q2, log = TRUE) # GP likelihood 2

  nll
}


# initial parameter list
par <- list(
  eta = rep(-2, 4),
  logit_delta1 = rep(0, 2),
  logit_delta2 = rep(0, 2),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  w1 = rep(0, nrow(spde1$c0)),
  w2 = rep(0, nrow(spde2$c0)),
  u = -5,
  mu = -0.5,
  z.star = rnorm(sum(is.na(data$y)), 0, 5)
)

# data list
dat <- list(
  y = data$y,
  spde1 = spde1,
  spde2 = spde2,
  bw = 15,
  trackID = data$trackID,
  trackInd = calc_trackInd(data$trackID)
)

t1 <- Sys.time()
obj <- MakeADFun(jnll, par, random = c("w1", "w2", "z.star"))

system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)
Sys.time()-t1

mod <- obj$report()
mod$kappa
mod$tau
mod$omega

states <- viterbi(mod = mod)
stateprobs <- stateprobs(mod = mod)
flare <- states != 1

# sdr <- sdreport(obj, getJointPrecision = TRUE)
# par <- as.list(sdr, "Est")

# cholP <- Cholesky(sdr$jointPrecision, LDL = FALSE, Imult = 1e-10)  # returns a "Cholesky" object
# pars <- lapply(1:1000, function(i){
#   z <- rnorm(nrow(sdr$jointPrecision))
#   p <- solve(cholP, z, system = "Lt")  # system="Lt" solves L^T x = z
#   p + c(sdr$par.fixed, sdr$par.random)
# })
# pars <- lapply(pars, obj3$env$parList)
#
# allstates <- matrix(NA, nrow(data), length(pars))
# allstateprobs <- array(dim = c(nrow(data), 3, length(pars)))
# for(i in 1:length(pars)){
#   getAll(pars[[i]], dat, warn = FALSE)
#   Gamma = diag(3)
#   Gamma[cbind(c(1:3, 3), c(2, 3, 1, 2))] = exp(eta)
#   Gamma = Gamma / rowSums(Gamma)
#   delta = c(1, exp(logit_delta))
#   delta = delta / sum(delta)
#   sigma = exp(log_sigma)
#   r = plogis(logit_r)
#   lambda = exp(log_lambda)
#   f <- f0 + mu
#   z <- y - f
#   n <- length(z); idx = 2:n
#   allprobs = matrix(0, n, 3)
#   allprobs[idx,1] = dnorm(z[idx], 0, sigma, log = TRUE)
#   allprobs[idx,2] = dexgauss(z[idx], z[idx-1], sigma, lambda, log = TRUE)
#   allprobs[idx,3] = dnorm(z[idx], r * z[idx-1], sigma, log = TRUE)
#   allstates[,i] <- viterbi(delta, Gamma, exp(allprobs))
#   allstateprobs[,,i] <- stateprobs(delta, Gamma, exp(allprobs))
# }
# state_proportions <- t(apply(allstates, 1, function(x) {
#   tab <- table(factor(x, levels = 1:3))
#   prop <- tab / length(x)
#   return(prop)
# }))
# stateprobs_unc <- apply(allstateprobs, c(1,2), mean, na.rm = TRUE)
prop_mat <- t(stateprobs)

### Plot result
# choose what to plot here
idx <- 7000:7350
# idx <- 7020:7150
# idx <- 11000:nrow(data)
# idx <- 4000:7500

# idx <- 1:nrow(data)
# prop_mat <- t(state_proportions[idx, ])  # rows = time, cols = states
# prop_mat <- t(stateprobs_unc[idx,])

# Stacked barplot of state probabilities
# par(mfrow = c(2,1), mar = c(5,4,2,2)+0.1)
# barplot(prop_mat[,idx], col = color, border = "white", space = 0,
#         ylab = "State proportion",
#         main = "Local state probabilities")

# Decoded time series
par(mfrow = c(1,1), mar = c(5,4,2,2)+0.1)
plot(data$time[idx], data$y[idx], col = color[states[idx]], pch = 16,
     xlab = "Time (sec)", ylab = "Flux", bty = "n", main = "Light Curve")
# lines(data$time[idx], mod_simple$f[idx], lwd = 2, lty = 3, col = "blue")
# y_aug <- data$y
# y_aug[is.na(data$y)] <- mod$y.star
# points(data$time[idx][is.na(data$y[idx])],
#       y_aug[idx][is.na(data$y[idx])], col = "blue")
lines(data$time[idx], mod$f[idx], lwd = 2, col = "plum")
legend("topleft", legend = c("Quiet", "Firing", "Decaying"), pch = 16, col = color, bty = "n")


# plot for paper
# Decoded time series
pdf("./figs/decoded_ts.pdf", width = 6, height = 2.7)
par(mfrow = c(1,1), mar = c(4,4,0.2,2)+0.1)
plot(data$time[idx], data$y[idx], col = color[states[idx]], pch = 16,
     xlab = "Time", ylab = "Flux", bty = "n")
# lines(data$time[idx], mod_simple$f[idx], lwd = 2, lty = 3, col = "blue")
lines(data$time[idx], mod$f[idx], lwd = 3, col = "plum")
legend("topright",
       legend = c("Quiet", "Firing", "Decaying", "Trend"),
       pch = c(rep(16, 3), NA),
       lwd = c(rep(NA, 3), 3),
       col = c(color, "plum"), bty = "n")
dev.off()



# simulate from fitted model to check acf
set.seed(123)
Gamma <- mod$Gamma
delta <- stationary(Gamma)
sigma <- mod$sigma
r <- mod$r
lambda <- mod$lambda

nObs <- 1e5
s <- z <- rep(NA, nObs)

s[1] <- sample(1:3, 1, prob = delta)
z[1] <- 0
for(t in 2:nObs) {
  s[t] <- sample(1:3, 1, prob = Gamma[s[t-1], ])
  if(s[t] == 1){
    z[t] <- rnorm(1, 0, sigma)
  }
  if(s[t] == 2){
    z[t] <- rexgauss(1, z[t-1], sigma, lambda)
  }
  if(s[t] == 3){
    z[t] <- r * z[t-1] + rnorm(1, 0, sigma)
  }
}

acf(z)

# simulate from fitted GP
sim <- rgmrf(10, mod$mu, mod$Q)

par(mfrow = c(1,1))
plot(sim[1,], type = "l", col = "#00000050", ylim = c(-15,15))
for(i in 2:nrow(sim)) lines(sim[i,], col = "#00000050")
