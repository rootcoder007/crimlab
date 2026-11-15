# Reference values of tests/testthat/test-krg.R: the geostatistics of
# rmorie (the package the functions come from) on the inputs the tests
# build, and those of gstat (semivariogram, point and block kriging,
# cross-validation), stats (trend surface), geoR (maximum likelihood and
# REML fits) and sf (the inverse Web Mercator projection). From the package
# root, with these packages in the library:
#   Rscript dev/ref/krg.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- c("rmorie", "gstat", "sp", "geoR", "sf")
source(file.path("dev", "ref", "put.R"))

# the inputs of test-krg.R
set.seed(41)
n <- 30
p <- cbind(x = runif(n), y = runif(n))
v <- sin(3 * p[, 1]) + cos(2 * p[, 2]) + rnorm(n, sd = 0.1)
at <- cbind(x = c(0.25, 0.5, 0.8), y = c(0.3, 0.55, 0.15))
vgm <- list(model = "exponential", nugget = 0.05, psill = 1, range = 0.3)
blk <- list(rbind(c(0, 0), c(0.5, 0.5)), rbind(c(0.5, 0.2), c(0.9, 0.6)))
z <- cbind(1, p[, 1])
# simulated Gaussian processes with an exponential covariance of range 0.2
set.seed(43)
q <- cbind(x = runif(40), y = runif(40))
D <- as.matrix(stats::dist(q))
e <- as.numeric(t(chol(exp(-D / 0.2) + diag(0.1, 40))) %*% stats::rnorm(40))
u0 <- 3 + e
u1 <- 1 + 2 * q[, 1] + as.numeric(t(chol(exp(-D / 0.2) + diag(1e-6, 40))) %*% stats::rnorm(40))
u2 <- 1 + 2 * q[, 1] + e

# "the semivariogram, the fits and kriging agree with rmorie"
vg <- rmorie::morie_spatial_variogram(p, v, n_bins = 10)
fit <- lapply(c("exponential", "spherical", "gaussian"), function(m) {
  f <- rmorie::morie_spatial_variogram_fit(q, u0, m)
  f[c("nugget", "psill", "range", "loglik")]
})
kr <- rmorie::morie_spatial_krige(p, v, at, vgm = vgm)
kw <- rmorie::morie_spatial_krige(p, v, at)
bk <- rmorie::morie_spatial_block_kriging(v, p, blk, n_quad = 16, nugget = 0,
                                          sill = 1, range_ = 0.3)
cv <- rmorie::morie_spatial_cross_validation(v, p, nugget = 0.05, sill = 1.05,
                                             range_ = 0.3)
tr <- rmorie::morie_spatial_trend_surface(v, p, order = 2)
sg <- rmorie::morie_spatial_glm(cbind(1, q[, 1]), u1, q)
sm <- rmorie::morie_spatial_mixed_model(cbind(1, q[, 1]), u2, q)
rm <- list(
  inp = c(sum(p), sum(v)),
  vg = as.list(vg),
  fit = fit,
  kr = as.list(kr),
  kw = as.list(kw),
  bk = bk[c("estimate", "se")],
  cv = cv$estimate,
  tr = tr[c("estimate", "se", "r2")],
  sg = sg[c("estimate", "se", "sigma2", "phi")],
  sm = sm[c("estimate", "se", "sigma2", "tau2", "phi")]
)

# "the spatial object agrees with rmorie"
d <- data.frame(lon = c(-79.4, -79.3, -79.5, -79.35), lat = c(43.6, 43.7, 43.65, 43.8),
                val = c(1, 2, 3, 4))
so <- rmorie::morie_spatial(d, c("lon", "lat"))
sw <- rmorie::morie_spatial_transform(so, 3857)
dm <- rmorie:::.sp_haversine(so$coords)
bw <- stats::median(dm[dm > 0])
dens <- rowSums(exp(-(dm / bw)^2))
pe <- as.matrix(stats::dist(sw$coords))
dp <- rowSums(exp(-(pe / stats::median(pe[pe > 0]))^2))
obj <- list(
  coords = unname(rmorie::morie_spatial_coords(so)),
  crs = rmorie::morie_spatial_crs(so),
  crsw = rmorie::morie_spatial_crs(sw),
  merc = unname(sw$coords),
  lf = lapply(rmorie::morie_spatial_leaflet_spec(sw)[c("center", "points")], unname),
  prob = unname(dens / sum(dens)),
  probm = unname(dp / sum(dp))
)

# "kriging, block kriging and cross-validation agree with gstat"
df <- data.frame(p, z = v)
sp::coordinates(df) <- ~ x + y
nd <- data.frame(at)
sp::coordinates(nd) <- ~ x + y
mod <- gstat::vgm(psill = 1, model = "Exp", range = 0.3, nugget = 0.05)
gk <- gstat::krige(z ~ 1, df, nd, model = mod, debug.level = 0)
br <- seq(0, max(stats::dist(p)) / 3, length.out = 11)
gv <- gstat::variogram(z ~ 1, df, boundaries = br)
m0 <- gstat::vgm(psill = 1, model = "Exp", range = 0.3)
# the 4 x 4 grids of the boxes, relative to their centres
g4 <- function(lo, hi) {
  e <- expand.grid(seq(lo[1], hi[1], length.out = 4), seq(lo[2], hi[2], length.out = 4))
  data.frame(x = e[, 1] - (lo[1] + hi[1]) / 2, y = e[, 2] - (lo[2] + hi[2]) / 2)
}
gb <- lapply(blk, function(b) {
  ctr <- data.frame(x = mean(b[, 1]), y = mean(b[, 2]))
  sp::coordinates(ctr) <- ~ x + y
  r <- gstat::krige(z ~ 1, df, ctr, model = m0, block = g4(b[1, ], b[2, ]),
                    debug.level = 0)
  c(r$var1.pred, r$var1.var)
})
gc <- gstat::krige.cv(z ~ 1, df, model = gstat::vgm(psill = 1, model = "Exp", range = 0.3,
                                                     nugget = 0.05), verbose = FALSE)
gst <- list(
  pred = gk$var1.pred, var = gk$var1.var,
  gamma = gv$gamma, np = gv$np,
  bk = gb,
  cv = gc$residual
)

# "the trend surface agrees with lm(), the fits with geoR and nlminb()"
l <- stats::lm(v ~ x + y + I(x^2) + I(y^2) + I(x * y), data = data.frame(p, v = v))
# vgfit() keeps the mean at the sample mean, as rmorie does, where geoR
# estimates it by GLS: its maximum is checked by another optimiser,
# nlminb(), on the same likelihood, from rmorie's estimate
nll <- function(lp) {
  C <- exp(lp[2]) * exp(-D / exp(lp[3])) + diag(exp(lp[1]) + 1e-8 * stats::var(u0), 40)
  L <- chol(C)
  a <- backsolve(L, u0 - mean(u0), transpose = TRUE)
  sum(log(diag(L))) + 0.5 * sum(a^2)
}
r0 <- rmorie::morie_spatial_variogram_fit(q, u0)
nm <- stats::nlminb(log(c(r0$nugget, r0$psill, r0$range)), nll)
gd1 <- geoR::as.geodata(cbind(q, u1))
gd2 <- geoR::as.geodata(cbind(q, u2))
lg <- NULL
utils::capture.output(lg <- suppressWarnings(geoR::likfit(gd1, trend = ~ q[, 1], ini.cov.pars = c(1, 0.2),
                                                         fix.nugget = TRUE, nugget = 0,
                                                         cov.model = "exponential", messages = FALSE)))
lr <- NULL
utils::capture.output(lr <- suppressWarnings(geoR::likfit(gd2, trend = ~ q[, 1], ini.cov.pars = c(1, 0.2),
                                                         nugget = 0.1, lik.method = "REML",
                                                         cov.model = "exponential", messages = FALSE)))
w3857 <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(data.frame(sw$coords), coords = c("lon", "lat"),
                                                        crs = 3857), 4326))
oth <- list(
  lm = unname(stats::coef(l)), lmse = unname(sqrt(diag(stats::vcov(l)))),
  ml = c(exp(nm$par), -nm$objective),
  glm = unname(c(lg$beta, lg$phi, lg$sigmasq)),
  reml = unname(c(lr$beta, lr$phi, lr$nugget / lr$sigmasq)),
  sf = unname(w3857)
)

put(list(rm = rm, obj = obj, gst = gst, oth = oth), "krg", need)
