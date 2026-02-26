library(LaMa)
library(RTMBdist)
library(fmesher)
library(viridis)
library(scales)
library(sf)

set.seed(38475)

# True spatial function for mean step length
true_field <- function(x,y){
  2 * (sin(2*pi*x/40) + cos(2*pi*y/40))
}

true_par <- list(
  delta = c(0.5, 0.5),
  mu = c(0.2, 5),
  sigma = c(0.5, 3),
  beta0 = qlogis(rep(0.2, 2))
)

curve(dgamma2(x, true_par$mu[1], true_par$sigma[1]), xlim = c(0, 8), n = 500,
      bty = "n", ylab = "Density", xlab = "Step length")
curve(dgamma2(x, true_par$mu[2], true_par$sigma[2]), add = TRUE, n = 500)


sim_data <- function(n, kappa_pull = 0.25, p = true_par) {
  s <- rep(NA, n)
  s[1] <- sample(1:2, size = 1, prob = p$delta) # first state
  loc <- matrix(0, n, 2) # initialise location matrix
  step <- rep(NA, n) # initialise step length vector
  phi <- 0 # initialise absolute heading
  for(t in 1:(n-1)) {
    # simulate step
    step[t] <- rgamma2(1, p$mu[s[t]], p$sigma[s[t]])
    # angle: add a pull towards c(0,0)
    pull_angle <- atan2(-loc[t,2], -loc[t,1]) - phi
    # draw absolute angle (independent of state)
    angle <- rvm(1, pull_angle, kappa_pull)
    # update heading
    phi <- phi + angle
    # next location
    loc[t+1,1] <- loc[t,1] + step[t] * cos(phi)
    loc[t+1,2] <- loc[t,2] + step[t] * sin(phi)
    # tpm evalutated at next location
    eta <- p$beta0
    eta[1] <- eta[1] + true_field(loc[t+1,1], loc[t+1,2])
    Gamma <- tpm(eta)
    # print(Gamma[2,1])
    # sample next state based on that tpm
    s[t+1] <- sample(1:2, 1, prob = Gamma[s[t], ])
  }
  data.frame(step = step, x = loc[,1], y = loc[,2], state = s)
}

data <- sim_data(10000)


# Make a prediction grid
x_seq <- seq(min(data$x), max(data$x), length.out = 512)
y_seq <- seq(min(data$y), max(data$y), length.out = 512)
grid_data <- expand.grid(x = x_seq, y = y_seq)
grid <- as.matrix(grid_data)
z <- outer(x_seq, y_seq, true_field)

par(mfrow = c(1,1))
image(x_seq, y_seq, z,
      xlab = "x", ylab = "y",
      col = viridis(25),
      main = "True spatial field", asp = 1, bty = "n")
lines(data$x, data$y, col = "white")
points(0, 0, col = "black", cex = 1, pch = 16)


jnll <- function(par) {
  getAll(par, dat, warn = FALSE)

  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)

  beta <- cbind(beta0, w); REPORT(w)
  Gamma <- tpm_g(cbind(1, X), beta)
  delta <- c(0.5, 0.5)

  lallprobs <- matrix(0, length(step), 2)
  ind <- which(!is.na(step))
  for(j in 1:2) {
    lallprobs[ind,j] <- dgamma2(step[ind], mu[j], sigma[j], log = TRUE)
  }

  nll <- -forward_g(delta, Gamma, lallprobs, bw = bw, logspace = TRUE)

  tau_sq <- exp(log_tau_sq); tau <- sqrt(tau_sq); REPORT(tau)
  kappa_sq <- exp(log_kappa_sq); kappa <- sqrt(kappa_sq); REPORT(kappa)
  Q <- tau_sq * (kappa_sq*kappa_sq * c0 + 2 * kappa_sq * g1 + g2)
  nll <- nll - sum(dgmrf(w, 0, Q, log = TRUE))

  nll
}

## Spatial part part of model
loc <- cbind(data$x, data$y)  # Spatial coordinates
bnd1 <- fm_nonconvex_hull(loc, convex = -0.025)
bnd2 <- fm_nonconvex_hull(loc, convex = -0.2)
mesh <- fm_mesh_2d(
  loc=loc,
  boundary=list(bnd1, bnd2),
  max.edge=c(10, 50),
  cutoff=2,
  plot.delay=0.5
)
# plot(mesh)
points(data$x, data$y)
spde <- fm_fem(mesh)
nrow(spde$c0)
X <- fm_basis(mesh, loc)

dat <- list(
  step = data$step,
  X = X,
  c0 = spde$c0, g1 = spde$g1, g2 = spde$g2,
  bw = 5
)

range <- max(dist(loc)) / 5
kappa0 <- sqrt(8)/range
sigma0 <- 3
tau0 <- 1/(sigma0 * sqrt(4*pi) * kappa0)

par <- list(
  log_mu = log(true_par$mu),
  log_sigma = log(true_par$sigma),
  beta0 = true_par$beta0,
  w = matrix(0, 2, nrow(spde$c0)),
  log_tau_sq = log(tau0^2),
  log_kappa_sq = log(kappa0^2)
)

environment(jnll) <- environment()

obj <- MakeADFun(jnll, par, random = "w")
opt <- nlminb(obj$par, obj$fn, obj$gr)

mod <- report(obj)

X_p <- fm_basis(mesh, grid)
fields <- X_p %*% t(mod$w)


# construct indicator which data points are inside bnd1
locs_sf <- st_as_sf(grid_data, coords = c("x", "y"), crs = st_crs(bnd1))
inside <- st_within(locs_sf, bnd1, sparse = FALSE)[,1]

z1 <- matrix(fields[,1], nrow = length(x_seq), ncol = length(y_seq))
z2 <- matrix(fields[,2], nrow = length(x_seq), ncol = length(y_seq))

g <- plogis(true_par$beta0[1] + z)
g1 <- plogis(true_par$beta0[1] + z1)
g2 <- plogis(true_par$beta0[1] + z2)

grange <- range(c(g, g1, g2), na.rm = TRUE)

par(mfrow = c(1,3))
image(x_seq, y_seq, g,
      xlab = "x", ylab = "y",
      col = viridis(100),
      zlim = grange,
      main = "True spatial field", asp = 1, bty = "n")
plot(st_geometry(bnd1), add = TRUE, border = "white", lty = 3)

image(x_seq, y_seq, g1,
      xlab = "x", ylab = "y",
      col = viridis(100),
      zlim = grange,
      main = "True spatial field", asp = 1, bty = "n")
plot(st_geometry(bnd1), add = TRUE, border = "white", lty = 3)

image(x_seq, y_seq, g2,
      xlab = "x", ylab = "y",
      col = viridis(100),
      zlim = grange,
      main = "True spatial field", asp = 1, bty = "n")
plot(st_geometry(bnd1), add = TRUE, border = "white", lty = 3)


g_inside <- as.numeric(g)[inside]
g1_inside <- as.numeric(g1)[inside]

cor(g_inside, g1_inside)

