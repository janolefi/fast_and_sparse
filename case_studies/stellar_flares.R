# installing dev version of packages
# devtools::install_github("janolefi/LaMa")
# devtools::install_github("kaskr/RTMB", subdir = "RTMB")
# devtools::install_github("janolefi/RTMBdist")

library(LaMa)       # for HMM functions
library(RTMBdist)   # for ExGaussian distribution
library(fmesher)    # for mesh and FEM matrices
library(Matrix)     # for sparse matrices

source("utils.R")

### colors for plotting
color = c("#00000070", "red", "orange")


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

# centering data
data$y <- scale(data$y, scale = FALSE)


### creating 2meshs and finite element matrices
mesh1 <- fm_mesh_1d(data$time[data$trackID == 1])
spde1 <- fm_fem(mesh1)

mesh2 <- fm_mesh_1d(data$time[data$trackID == 2])
spde2 <- fm_fem(mesh2)

### Simple analysis without state switching
nll0 <- function(par) {
  getAll(par, dat)

  sigma <- exp(log_sigma); REPORT(sigma)
  f1 <- w1 + mu
  f2 <- w2 + mu
  f <- c(f1, f2); REPORT(f) # total smooth
  nll <- -sum(dnorm(y, f, sigma, log = TRUE), na.rm = TRUE) # observation liklihood

  ### GP part ###
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  cos_pi_omega <- 2 * plogis(u) - 1; omega <- acos(cos_pi_omega) / pi; REPORT(omega)

  Q1 <- tau_sq * (kappa_sq*kappa_sq * spde1$c0 + 2 * cos_pi_omega * kappa_sq * spde1$g1 + spde1$g2)
  Q2 <- tau_sq * (kappa_sq*kappa_sq * spde2$c0 + 2 * cos_pi_omega * kappa_sq * spde2$g1 + spde2$g2)

  nll <- nll - dgmrf(w1, 0, Q1, log = TRUE) # GP likelihood 1
  nll <- nll - dgmrf(w2, 0, Q2, log = TRUE) # GP likelihood 2

  nll
}


# initial parameter list
par <- list(
  log_sigma = log(7),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  u = -5,
  mu = -0.5,
  w1 = numeric(nrow(spde1$c0)),
  w2 = numeric(nrow(spde2$c0))
)

# data list
dat <- list(
  y = data$y,
  spde1 = spde1,
  spde2 = spde2
)

obj0 <- MakeADFun(nll0, par, random = c("w1", "w2"))
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
mod0 <- obj0$report()
mod0$tau
mod0$kappa
mod0$omega

### visualising results
idx = 6000:8000
plot(data$time[idx], data$y[idx], pch = 16, col = color[1],
     xlab = "Time (days)", ylab = "Flux", bty = "n")
lines(data$time[idx], mod0$f[idx], lwd = 3, col = "blue")
# smooth is not very periodic
# erratic behaviour whenever there is a flare


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
  Delta <- cbind(1, exp(logit_Delta))
  Delta <- Delta / rowSums(Delta)

  #### state-dependent process ####
  # parameter transformations
  sigma <- exp(log_sigma); REPORT(sigma)
  r <- plogis(logit_r); REPORT(r)
  lambda <- exp(log_lambda); REPORT(lambda)
  # state-dependent densities
  f1 <- w1 + mu
  f2 <- w2 + mu
  f <- c(f1, f2); REPORT(f) # quasi-periodic smooth
  z <- y - f; REPORT(mu)

  # latent z's where observations are mising -> integrated out
  z[is.na(y)] <- z.star; REPORT(z)

  n <- length(z)
  idx <- list(2:(trackInd[2]-1), (trackInd[2]+1):n) # idx of the 2 time series
  lallprobs <- matrix(0, n, 3)
  for(i in 1:2) { # loop over time series
    ind <- idx[[i]]
    # regular measurement error:
    lallprobs[ind,1] <- dnorm(z[ind], 0, sigma, log = TRUE)
    # firing:
    lallprobs[ind,2] <- dexgauss(z[ind], z[ind-1], sigma, lambda, log = TRUE)
    # decaying:
    lallprobs[ind,3] <- dnorm(z[ind], r * z[ind-1], sigma, log = TRUE)
  }

  ### HMM likelihood
  nll <- -forward(Delta, Gamma, lallprobs,
                  trackID = trackID, bw = bw, logspace = TRUE)

  ### GP likelihood
  # parameter transformations
  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  cos_pi_omega <- 2 * plogis(u) - 1; omega <- acos(cos_pi_omega)/pi; REPORT(omega)

  Q1 <- tau_sq * (kappa_sq*kappa_sq * spde1$c0 + 2 * cos_pi_omega * kappa_sq * spde1$g1 + spde1$g2)
  Q2 <- tau_sq * (kappa_sq*kappa_sq * spde2$c0 + 2 * cos_pi_omega * kappa_sq * spde2$g1 + spde2$g2)

  nll <- nll - dgmrf(w1, 0, Q1, log = TRUE) # GP likelihood 1
  nll <- nll - dgmrf(w2, 0, Q2, log = TRUE) # GP likelihood 2

  nll
}


# initial parameter list
par <- list(
  eta = rep(-2, 4),
  logit_Delta = matrix(0, 2, 2),
  log_sigma = log(7),
  logit_r = qlogis(0.8),
  log_lambda = log(0.025),
  log_tau_sq = log(0.005^2),
  log_kappa_sq = log(30^2),
  w1 = rep(0, nrow(spde1$c0)),
  w2 = rep(0, nrow(spde2$c0)),
  u = -5,
  mu = -0.5,
  z.star = rnorm(sum(is.na(data$y)), 0, 3)
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
Sys.time() - t1

mod <- obj$report()
mod$kappa
mod$tau
mod$omega

states <- viterbi(mod = mod)
stateprobs <- stateprobs(mod = mod)
flare <- states != 1


### Plot result
# choose what to plot here
idx <- 7030:7350

pdf("./figs/flare_result.pdf", width = 7, height = 4.5)

# Stacked barplot of state probabilities
layout(matrix(1:2, ncol = 1), heights = c(0.6, 1))
# par(mfrow = c(2,1), mar = c(5,4,2,2)+0.1)

par(mar = c(1, 4, 1, 2))
barplot(t(stateprobs[idx,]), col = c("black", "red", "orange"), border = "white", space = 0,
        ylab = "State probability",
        main = "", yaxt = "n")
axis(2, at = c(0, 0.5, 1), labels = c(0, 0.5, 1))

# Decoded time series
par(mar = c(5, 4, 0.5, 2), xpd = NA)
plot(data$time[idx], data$y[idx], col = color[states[idx]], pch = 16,
     xlab = "Time (days)", ylab = "Flux", bty = "n")
lines(data$time[idx], mod$f[idx], lwd = 3, col = "plum")

# legend
legend(x = 1335.415, y = 190,
       legend = c("Quiet", "Firing", "Decaying", "Trend"),
       pch = c(rep(16, 3), NA),
       lwd = c(rep(NA, 3), 3),
       col = c(color, "plum"), bty = "n")

dev.off()



# simulate from fitted model to check acf
set.seed(123)
Gamma <- mod$Gamma
delta <- mod$delta[1,]
sigma <- mod$sigma
r <- mod$r
lambda <- mod$lambda

nObs <- 1e6
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

# Look at PACF: empirical measure for conditional dependence
pacf(z)
# significant lags up to 5 -> bw = 15 seems enough

# simulate from fitted GP
Q <- mod$tau^2 * (mod$kappa^4 * spde1$c0 + 2 * cos(pi * mod$omega) * mod$kappa^2 * spde1$g1 + spde1$g2)
sim <- rgmrf(10, mod$mu, Q)

par(mfrow = c(1,1))
plot(sim[1,], type = "l", col = "#00000050", ylim = c(-15,15))
for(i in 2:nrow(sim)) lines(sim[i,], col = "#00000050")
