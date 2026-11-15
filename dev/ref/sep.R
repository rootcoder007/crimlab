# Reference values of tests/testthat/test-sep.R: the spatial epidemiology
# of rmorie (the package the functions come from) on the inputs the tests
# build, and those of surveillance (the Bayes detector), fields
# (great-circle distances), MASS (generalised inverse, negative binomial),
# stats (Poisson), pscl (zero-inflated Poisson), deSolve (Runge-Kutta),
# survival (the Weibull model) and spdep (Moran's I). From the package
# root, with these packages in the library:
#   Rscript dev/ref/sep.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("rmorie", "surveillance", "fields", "MASS", "pscl", "deSolve", "survival", "spdep")
source(file.path("dev", "ref", "put.R"))

# the inputs of test-sep.R
set.seed(51)
obs <- rpois(70, 6) + c(rep(0, 66), 12, 0, 9, 0)
rec <- cbind(x = runif(6), y = runif(6))
src <- cbind(x = runif(5), y = runif(5))
sw <- c(2, 1, 3, 1, 0.5)
ll <- cbind(lon = -79.4 + runif(4, -0.1, 0.1), lat = 43.7 + runif(4, -0.1, 0.1))
ls <- cbind(lon = -79.4 + runif(3, -0.1, 0.1), lat = 43.7 + runif(3, -0.1, 0.1))
A <- rbind(c(0, 1, 1, 0, 0), c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1), c(0, 1, 1, 0, 1), c(0, 0, 1, 1, 0))
ne <- 40
ex <- cbind(x1 = rnorm(ne), x2 = runif(ne))
ee <- runif(ne, 3, 8)
ey <- rpois(ne, ee * exp(0.2 + 0.4 * ex[, 1] - 0.5 * ex[, 2])) * rbinom(ne, 1, 0.8)
hop <- rbind(c(0, 1, 2, -1), c(1, 0, 1, -1), c(2, 1, 0, -1), c(-1, -1, -1, 0))
W3 <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
i <- 0:29
fx <- matrix(sin(i))
ft <- 0.5 + ((i * 7) %% 11) / 5 * exp(-0.5 * fx[, 1] - 0.3 * (i %% 3))
fe <- as.numeric(i %% 5 != 0)
j <- 0:35
cx <- matrix(sin(j))
ct <- ifelse(j %% 4 == 0, 5, 0.2 + ((j * 7) %% 11) / 8)
ce <- as.numeric(!(j %% 4 == 0 | j %% 9 == 1))
P <- 1 * (abs(outer(1:12, 1:12, "-")) == 1)
P[1, 12] <- P[12, 1] <- 1
P[3, 9] <- P[9, 3] <- 1
fv <- sin(1:12 / 2) + 0.3 * cos(1:12)
fz <- cos(1:12 / 3)

# "the functions agree with rmorie"
fr <- rmorie::spatial_frailty(ft, fe, fx, i %% 3, W3)
cr <- rmorie::spatial_cure_rate(ct, ce, cx, j %% 3, W3)
sf1 <- rmorie::spatial_filter(fv, P, method = "moran")
sf2 <- rmorie::spatial_filter(fv, P, X = matrix(fz), method = "aic")
sf3 <- rmorie::spatial_filter(fv, P, method = "r2", tol = 0.01)
sir <- rmorie::spatial_sir(c(5, 0, 0), W3, c(1000, 800, 1200), beta = 0.5, gamma = 0.2, t_max = 20, dt = 0.5)
by <- rmorie::BayesOutbreak(obs, 52, 0, 6)
b2 <- rmorie::BayesOutbreak(obs, 10, 2, 2, act_y = FALSE, alpha = 0.1, time_points = 60:70)
rm <- list(
  inp = c(sum(obs), sum(rec), sum(src), sum(ll), sum(ls), sum(ex), sum(ee), sum(ey)),
  bay = by[c("upperbound", "alarm")],
  bay2 = b2[c("upperbound", "alarm")],
  buf = rmorie::BufferExposure(rec, src, 0.4, sw),
  bufll = rmorie::BufferExposure(ll[, 2:1], ls[, 2:1], 8, metric = "haversine"),
  kg = rmorie::KernelExposure(rec, src, 0.3, sw),
  kq = rmorie::KernelExposure(rec, src, 0.5, kernel = "quartic"),
  ler = rmorie::LerouxPrecision(A, 0.7, 2),
  bym = rmorie::Bym2Structure(A, 2, 0.4),
  pois = rmorie::EcologicalRegression(ey, ee, ex)[c("coefficients", "fitted", "loglik")],
  nb = rmorie::EcologicalRegression(ey, ee, ex, "negbin")[c("coefficients", "theta", "loglik")],
  zip = rmorie::EcologicalRegression(ey, ee, ex, "zip", Z = ex[, 1, drop = FALSE])[
    c("coefficients", "zero_coefficients", "loglik")],
  bias = rmorie::spatial_bias(hop, c(1, 0.5, 0.25), -5)$bias,
  sir = list(attack = sir$attack_rate, S = do.call(rbind, sir$S), I = do.call(rbind, sir$I),
             R = do.call(rbind, sir$R)),
  fr = fr[c("statistic", "shape", "coefficients", "frailties", "se", "iterations")],
  cr = cr[c("statistic", "shape", "log_scale", "coefficients", "regional_effects", "se", "cure_fraction")],
  sf1 = sf1[c("selected", "coefficients", "residual_moran", "aic", "r2", "filter_r2", "residuals")],
  sf2 = sf2[c("selected", "coefficients", "residual_moran", "aic", "r2", "filter_r2")],
  sf3 = sf3[c("selected", "coefficients", "r2")],
  mi = rmorie::MoranEigenvectors(P)$moran_i
)

# "the functions agree with the reference packages"
dp <- surveillance::create.disProg(week = seq_along(obs), observed = obs, state = rep(0, length(obs)),
                                   freq = 52)
ab <- suppressWarnings(surveillance::algo.bayes(dp, control = list(range = 7:70, b = 0, w = 6, actY = TRUE,
                                                                    alpha = 0.05)))
dll <- fields::rdist.earth(ll, ls, miles = FALSE, R = 6371.0088)
Q <- diag(rowSums(A)) - A
G <- MASS::ginv(Q)
df <- data.frame(y = ey, ex, e = ee)
gp <- stats::glm(y ~ x1 + x2 + offset(log(e)), family = stats::poisson(), data = df,
                 control = stats::glm.control(epsilon = 1e-14, maxit = 100))
gn <- MASS::glm.nb(y ~ x1 + x2 + offset(log(e)), data = df,
                   control = stats::glm.control(epsilon = 1e-14, maxit = 200))
gz <- pscl::zeroinfl(y ~ x1 + x2 + offset(log(e)) | x1, data = df, dist = "poisson",
                     control = pscl::zeroinfl.control(reltol = 1e-14, maxit = 5000))
N3 <- c(1000, 800, 1200)
C3 <- 0.9 * diag(3) + 0.1 * W3 / rowSums(W3)
rhs <- function(t, y, parms) {
  S <- y[1:3]
  I <- y[4:6]
  f <- as.vector(C3 %*% (I / N3))
  list(c(-0.5 * S * f, 0.5 * S * f - 0.2 * I, 0.2 * I))
}
od <- deSolve::ode(c(N3 - c(5, 0, 0), c(5, 0, 0), rep(0, 3)), seq(0, 20, by = 0.5), rhs, NULL,
                   method = "rk4")
sv <- survival::survreg(survival::Surv(ft, fe) ~ fx[, 1], dist = "weibull")
fr1 <- rmorie::spatial_frailty(ft, fe, fx, rep(0, 30), matrix(0, 1, 1), tau = 0)
lw <- spdep::mat2listw(P, style = "B")
sres <- stats::lm.fit(cbind(1, rmorie::MoranEigenvectors(P)$vectors[, sf1$selected + 1, drop = FALSE]), fv)
oth <- list(
  bay = as.vector(ab$upperbound), bayal = as.logical(ab$alarm),
  near = apply(dll, 1, min),
  ginv = G, scale = exp(mean(log(diag(G)))),
  pois = unname(stats::coef(gp)), poll = as.numeric(stats::logLik(gp)),
  nb = unname(stats::coef(gn)), theta = gn$theta,
  zip = unname(c(gz$coefficients$count, gz$coefficients$zero)), zll = as.numeric(stats::logLik(gz)),
  rk4 = unname(od[, -1]),
  weib = c(1 / sv$scale, -unname(stats::coef(sv)) / sv$scale),
  fr1 = c(fr1$shape, fr1$coefficients, fr1$frailties),
  moran = spdep::moran(sres$residuals, lw, 12, spdep::Szero(lw))$I,
  lmcoef = unname(sres$coefficients)
)

put(list(rm = rm, oth = oth), "sep", need)
