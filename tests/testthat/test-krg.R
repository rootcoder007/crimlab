# Geostatistics against the values of rmorie and of gstat, stats, geoR and
# sf, computed once by dev/ref/krg.R and stored in ref/krg.txt.

kin <- function() {
  set.seed(41)
  n <- 30
  p <- cbind(x = runif(n), y = runif(n))
  v <- sin(3 * p[, 1]) + cos(2 * p[, 2]) + rnorm(n, sd = 0.1)
  list(
    p = p,
    v = v,
    at = cbind(x = c(0.25, 0.5, 0.8), y = c(0.3, 0.55, 0.15)),
    vgm = list(model = "exponential", nugget = 0.05, psill = 1, range = 0.3),
    blk = list(rbind(c(0, 0), c(0.5, 0.5)), rbind(c(0.5, 0.2), c(0.9, 0.6)))
  )
}

# simulated Gaussian processes with an exponential covariance of range 0.2
gpin <- function() {
  set.seed(43)
  q <- cbind(x = runif(40), y = runif(40))
  D <- as.matrix(dist(q))
  e <- as.numeric(t(chol(exp(-D / 0.2) + diag(0.1, 40))) %*% rnorm(40))
  list(
    q = q,
    u0 = 3 + e,
    u1 = 1 +
      2 * q[, 1] +
      as.numeric(t(chol(exp(-D / 0.2) + diag(1e-6, 40))) %*% rnorm(40)),
    u2 = 1 + 2 * q[, 1] + e
  )
}

# the Gaussian log-likelihood (less n / 2 log(2 pi)) of a variogram model
# about the sample mean, with vgfit()'s jitter of 1e-8 of the variance
gll <- function(q, u, m, nug, ps, rg) {
  h <- as.matrix(dist(q)) / rg
  s <- switch(
    m,
    exponential = exp(-h),
    spherical = ifelse(h >= 1, 0, 1 - 1.5 * h + 0.5 * h^3),
    gaussian = exp(-h^2)
  )
  C <- ps * s + diag(nug + 1e-8 * var(u), length(u))
  L <- chol(C)
  -(sum(log(diag(L))) +
    0.5 * sum(backsolve(L, u - mean(u), transpose = TRUE)^2))
}

test_that("the semivariogram, the fits and kriging agree with rmorie", {
  k <- kin()
  ref <- refget("krg")$rm
  expect_equal(
    c(sum(k$p), sum(k$v)),
    ref$inp,
    tolerance = 1e-12,
    info = "inputs changed"
  )
  g <- vgemp(k$p, v = k$v, nbins = 10)
  expect_equal(g$dist, ref$vg$dist, tolerance = 1e-12)
  expect_equal(g$gamma, ref$vg$gamma, tolerance = 1e-12)
  expect_identical(g$np, ref$vg$np)
  # rmorie's likelihood takes ||C^-1 r||^2 (backsolve() after
  # forwardsolve()) where the Gaussian likelihood has r' C^-1 r, so its fit
  # is not the maximum: the Gaussian likelihood is higher at vgfit()'s fit
  g <- gpin()
  for (m in 1:3) {
    mo <- c("exponential", "spherical", "gaussian")[m]
    f <- vgfit(g$q, v = g$u0, model = mo)
    r <- ref$fit[[m]]
    expect_equal(
      f$loglik,
      gll(g$q, g$u0, mo, f$nugget, f$psill, f$range),
      tolerance = 1e-10
    )
    expect_gt(f$loglik, gll(g$q, g$u0, mo, r$nugget, r$psill, r$range))
  }
  r <- okrge(k$p, v = k$v, at = k$at, vgm = k$vgm)
  expect_equal(r$pred, ref$kr$pred, tolerance = 1e-10)
  expect_equal(r$var, ref$kr$var, tolerance = 1e-10)
  # the default: an exponential model by weighted least squares, then kriging
  w <- okrge(k$p, v = k$v, at = k$at)
  expect_equal(w$pred, ref$kw$pred, tolerance = 1e-6)
  expect_equal(w$var, ref$kw$var, tolerance = 1e-6)
  # rmorie takes the distances between quadrature points as sqrt(|p|^2 +
  # |q|^2 - 2 p.q), which leaves about 1e-9 where they coincide: its
  # standard errors differ by about 2e-9 relative (gstat's agree to 1e-10)
  b <- bkrig(k$p, v = k$v, blocks = k$blk, nq = 16, range = 0.3)
  expect_equal(b$est, ref$bk$estimate, tolerance = 1e-10)
  expect_equal(b$se, ref$bk$se, tolerance = 1e-8)
  cv <- krigcv(k$p, v = k$v, nugget = 0.05, sill = 1.05, range = 0.3)
  expect_equal(cv$resid, ref$cv$residuals, tolerance = 1e-10)
  expect_equal(cv$mspe, ref$cv$MSPE, tolerance = 1e-10)
  expect_equal(cv$rmspe, ref$cv$RMSPE, tolerance = 1e-10)
  expect_equal(cv$mae, ref$cv$MAE, tolerance = 1e-10)
  t <- trsurf(k$p, v = k$v)
  expect_equal(unname(t$coef), ref$tr$estimate, tolerance = 1e-10)
  expect_equal(unname(t$se), ref$tr$se, tolerance = 1e-10)
  expect_equal(t$r2, ref$tr$r2, tolerance = 1e-10)
  z <- cbind(1, g$q[, 1])
  # Brent's method as optimize() runs it: phi to its tolerance
  s <- sglmx(g$u1, z, g$q)
  expect_equal(s$phi, ref$sg$phi, tolerance = 1e-6)
  expect_equal(unname(s$coef), ref$sg$estimate, tolerance = 1e-6)
  expect_equal(unname(s$se), ref$sg$se, tolerance = 1e-6)
  expect_equal(s$sigma2, ref$sg$sigma2, tolerance = 1e-6)
  m <- smixm(g$u2, z, g$q)
  expect_equal(unname(m$coef), ref$sm$estimate, tolerance = 1e-6)
  expect_equal(unname(m$se), ref$sm$se, tolerance = 1e-6)
  expect_equal(
    c(m$sigma2, m$tau2, m$phi),
    c(ref$sm$sigma2, ref$sm$tau2, ref$sm$phi),
    tolerance = 1e-6
  )
})

test_that("kriging, block kriging and cross-validation agree with gstat", {
  k <- kin()
  ref <- refget("krg")$gst
  r <- okrge(k$p, v = k$v, at = k$at, vgm = k$vgm)
  expect_equal(r$pred, ref$pred, tolerance = 1e-10)
  expect_equal(r$var, ref$var, tolerance = 1e-10)
  g <- vgemp(k$p, v = k$v, nbins = 10)
  expect_equal(g$gamma, ref$gamma, tolerance = 1e-12)
  expect_equal(g$np, ref$np)
  b <- bkrig(k$p, v = k$v, blocks = k$blk, nq = 16, range = 0.3)
  expect_equal(b$est, vapply(ref$bk, `[`, 0, 1), tolerance = 1e-10)
  expect_equal(b$se^2, vapply(ref$bk, `[`, 0, 2), tolerance = 1e-10)
  cv <- krigcv(k$p, v = k$v, nugget = 0.05, sill = 1.05, range = 0.3)
  expect_equal(cv$resid, ref$cv, tolerance = 1e-10)
})

test_that("the trend surface agrees with lm(), the fits with geoR, nlminb()", {
  k <- kin()
  ref <- refget("krg")$oth
  t <- trsurf(k$p, v = k$v)
  expect_equal(unname(t$coef), ref$lm, tolerance = 1e-10)
  expect_equal(unname(t$se), ref$lmse, tolerance = 1e-10)
  g <- gpin()
  q <- g$q
  u0 <- g$u0
  u1 <- g$u1
  u2 <- g$u2
  # the maximum found by nlminb() on the same likelihood
  f <- vgfit(q, v = u0)
  expect_equal(f$loglik, ref$ml[4], tolerance = 1e-6)
  expect_equal(c(f$nugget, f$psill, f$range), ref$ml[1:3], tolerance = 1e-2)
  # likfit() maximises the same profile likelihood with optim()'s L-BFGS-B;
  # its variance divides by n, sglmx()'s by n - p
  s <- sglmx(u1, cbind(1, q[, 1]), q)
  expect_equal(c(unname(s$coef), s$phi), ref$glm[1:3], tolerance = 1e-3)
  expect_equal(s$sigma2 * 38 / 40, ref$glm[4], tolerance = 1e-3)
  # REML: Nelder-Mead stops at a relative change of 1e-6 in the likelihood
  m <- smixm(u2, cbind(1, q[, 1]), q)
  expect_equal(
    c(unname(m$coef), m$phi, m$tau2 / m$sigma2),
    ref$reml,
    tolerance = 1e-2
  )
})

test_that("the spatial object agrees with rmorie and sf", {
  ref <- refget("krg")$obj
  d <- data.frame(
    lon = c(-79.4, -79.3, -79.5, -79.35),
    lat = c(43.6, 43.7, 43.65, 43.8),
    val = c(1, 2, 3, 4)
  )
  s <- spobj(d)
  expect_s3_class(s, "spobj")
  expect_equal(unname(spxy(s)), ref$coords, tolerance = 0)
  expect_identical(colnames(spxy(s)), c("lon", "lat"))
  expect_identical(crsof(s), ref$crs)
  expect_identical(names(s$data), "val")
  m <- spobj(
    data.frame(e = ref$merc[, 1], n = ref$merc[, 2], val = 1:4),
    c("e", "n"),
    3857
  )
  expect_identical(crsof(m), ref$crsw)
  expect_identical(m$units, "metres")
  expect_identical(crsof(spobj(d[1:2], crs = 2958))$wkt, "EPSG:2958")
  l <- lfspec(m)
  expect_equal(unname(l$points), ref$lf$points, tolerance = 1e-12)
  expect_equal(unname(l$center), unname(ref$lf$center), tolerance = 1e-12)
  expect_equal(unname(l$points), refget("krg")$oth$sf, tolerance = 1e-10)
  expect_equal(unname(lfspec(s)$points), ref$coords, tolerance = 0)
  expect_identical(l$zoom, 8L)
  # the sampling probabilities of rmorie's density weights
  expect_equal(
    attr(spsamp(s, 2, seed = 1), "prob"),
    ref$prob,
    tolerance = 1e-10
  )
  expect_equal(
    attr(spsamp(m, 2, seed = 1), "prob"),
    ref$probm,
    tolerance = 1e-10
  )
  expect_output(print(s), "spobj: 4 points, EPSG:4326 \\(degrees\\)")
  expect_identical(format(m), "spobj: 4 points, EPSG:3857 (metres)")
  expect_identical(summary(s)$min, c(-79.5, 43.6))
})

test_that("sampling draws without replacement, by density or uniformly", {
  s <- spobj(expand.grid(lon = 1:5, lat = 1:5))
  a <- spsamp(s, 25, seed = 3)
  expect_setequal(as.vector(a), 1:25)
  expect_identical(spsamp(s, 5, seed = 3), spsamp(s, 5, seed = 3))
  u <- spsamp(s, 4, dens = FALSE, seed = 9)
  expect_equal(attr(u, "prob"), rep(1 / 25, 25))
  set.seed(1)
  r1 <- spsamp(s, 3)
  set.seed(1)
  expect_identical(spsamp(s, 3), r1)
  # the corners of a grid are sparser than its centre
  pr <- attr(a, "prob")
  expect_lt(pr[1], pr[13])
  expect_length(spsamp(s, 0, seed = 1), 0)
})

test_that("rmorie's structural checks hold", {
  # a simulated Gaussian process with an exponential covariance
  set.seed(70)
  n <- 80
  xy <- cbind(x = runif(n), y = runif(n))
  C <- exp(-as.matrix(dist(xy)) / 0.3) + diag(0.01, n)
  y <- as.numeric(t(chol(C)) %*% rnorm(n)) + 5
  vg <- vgemp(xy, v = y)
  expect_true(all(vg$np > 0))
  expect_lt(mean(vg$gamma[1:3]), mean(vg$gamma[(nrow(vg) - 2):nrow(vg)]))
  fit <- vgfit(xy, v = y)
  expect_true(is.finite(fit$loglik))
  expect_gt(fit$range, 0.05)
  expect_lt(fit$range, 1.5)
  ll0 <- sum(dnorm(y, mean(y), sd(y), log = TRUE)) + n / 2 * log(2 * pi)
  expect_gt(fit$loglik, ll0)
  kr <- okrge(xy, v = y, at = xy[1:10, ], vgm = fit)
  expect_equal(kr$pred, y[1:10], tolerance = 0.15)
  k2 <- okrge(xy, v = y, at = cbind(0.5, 0.5), vgm = fit)
  expect_true(k2$pred > min(y) && k2$pred < max(y))
  k3 <- okrge(xy, v = y, at = cbind(5, 5), vgm = fit)
  expect_gt(k3$var, k2$var)
  # a line is fitted exactly (rmorie's canonical case, on two dimensions)
  p <- cbind(x = 0:4, y = c(0, 1, 0, 1, 0))
  expect_equal(
    unname(trsurf(p, v = 1 + (0:4), order = 1)$coef),
    c(1, 1, 0),
    tolerance = 1e-12
  )
  s <- sglmx(1 + 2 * (0:4), cbind(1, 0:4), p)
  expect_equal(unname(s$coef), c(1, 2), tolerance = 1e-6)
  m <- smixm(1 + 2 * (0:4) + c(0.1, -0.1, 0.05, 0, -0.05), cbind(1, 0:4), p)
  expect_equal(unname(m$coef), c(1, 2), tolerance = 0.1)
})

test_that("results print and summarise", {
  k <- kin()
  expect_output(print(vgfit(k$p, v = k$v)), "vgfit: exponential model")
  expect_output(
    print(krigcv(k$p, v = k$v, range = 0.3)),
    "krigcv: leave-one-out kriging of 30 values"
  )
  t <- trsurf(k$p, v = k$v, order = 3)
  expect_output(print(t), "trsurf: trend surface of order 3 on 30 points")
  st <- summary(t)
  expect_identical(
    st$term,
    c("1", "x", "y", "x2", "y2", "xy", "x3", "y3", "x2y", "xy2")
  )
  expect_equal(st$z, st$est / st$se)
  s <- sglmx(k$v, NULL, k$p)
  expect_identical(names(s$coef), "1")
  expect_output(print(s), "sglmx: spatial linear model")
  expect_named(summary(s), c("term", "est", "se", "z"))
  zz <- cbind(a = 1, b = k$p[, 1])
  m <- smixm(k$v, zz, k$p)
  expect_identical(names(m$coef), c("a", "b"))
  expect_output(
    print(m),
    "smixm: spatial linear mixed model by REML, 30 points"
  )
  expect_identical(summary(m)$term, c("a", "b"))
  expect_identical(names(sglmx(k$v, k$p[, 1], k$p)$coef), "z1")
  expect_identical(names(sglmx(k$v, as.data.frame(zz), k$p)$coef), c("a", "b"))
  expect_identical(nrow(bkrig(k$p, v = k$v, blocks = list(k$p[1:3, ]))), 1L)
})

test_that("geostatistics refuse bad input in words", {
  k <- kin()
  p <- k$p
  v <- k$v
  errs(
    "vgemp",
    list(
      list("at least 2 points", p[1, , drop = FALSE], NULL, 1),
      list("`nbins` must be one whole number from 1", p, NULL, v, 0),
      list(
        "`cutoff` must be one number strictly between 0 and Inf",
        p,
        NULL,
        v,
        15,
        -1
      ),
      list(
        "the points all coincide",
        cbind(x = c(1, 1), y = c(2, 2)),
        NULL,
        c(1, 2)
      )
    )
  )
  errs(
    "vgfit",
    list(
      list("at least 3 points", p[1:2, ], NULL, v[1:2]),
      list("`model` must be one of", p, NULL, v, "linear"),
      list("`v` must vary", p, NULL, rep(1, 30))
    )
  )
  errs(
    "okrge",
    list(
      list("at least 2 points to fit", p[1, , drop = FALSE], NULL, 1, p),
      list("`vgm` must be a list", p, NULL, v, p, 1),
      list(
        "`vgm\\$range` must be greater than 0",
        p,
        NULL,
        v,
        p,
        list(model = "exponential", nugget = 0, psill = 1, range = 0)
      ),
      list(
        "singular",
        rbind(c(0, 0), c(0, 0)),
        NULL,
        c(1, 2),
        p,
        list(model = "exponential", nugget = 0, psill = 1, range = 1)
      ),
      list("close enough", cbind(x = c(1, 1), y = c(2, 2)), NULL, c(1, 2), p)
    )
  )
  errs(
    "bkrig",
    list(
      list("`blocks` must be a list of blocks", p, NULL, v, 1),
      list(
        "`blocks\\[\\[1\\]\\]` must be a numeric matrix",
        p,
        NULL,
        v,
        list(1:3)
      ),
      list(
        "`blocks\\[\\[2\\]\\]` must hold finite numbers",
        p,
        NULL,
        v,
        list(p[1:3, ], rbind(c(0, NA), c(1, 1)))
      ),
      list(
        "`sill` must be at least `nugget`",
        p,
        NULL,
        v,
        list(p[1:3, ]),
        25,
        2,
        1
      ),
      list("singular", rbind(c(0, 0), c(0, 0)), NULL, c(1, 2), list(p[1:3, ]))
    )
  )
  errs(
    "krigcv",
    list(
      list("at least 3 points", p[1:2, ], NULL, v[1:2]),
      list(
        "`range` must be one number strictly between 0 and Inf",
        p,
        NULL,
        v,
        0,
        1,
        0
      ),
      list("singular", rbind(c(0, 0), c(0, 0), c(1, 1)), NULL, c(1, 2, 3))
    )
  )
  errs(
    "trsurf",
    list(
      list("`order` must be 0, 1, 2 or 3", p, NULL, v, 4),
      list("order 3 needs at least 10 points", p[1:5, ], NULL, v[1:5], 3),
      list("collinear", cbind(x = 1:4, y = 1:4), NULL, 1:4, 1)
    )
  )
  errs(
    "sglmx",
    list(
      list("`z` must have 30 rows", v, cbind(1, 1:3), p),
      list("`z` must be a numeric matrix or vector", v, "a", p),
      list("`z` must hold finite numbers", v, c(NA, 2:30), p),
      list("more values than `z` has columns", v[1:2], diag(2), p[1:2, ]),
      list("collinear", v, cbind(1, rep(1, 30)), p)
    )
  )
  errs("smixm", list(list("collinear", v, cbind(1, rep(1, 30)), p)))
  errs(
    "spobj",
    list(
      list("`data` must be a data frame", 1),
      list("`coords` must name two columns", data.frame(a = 1, b = 2)),
      list(
        "`crs` must be one EPSG code",
        data.frame(lon = 1, lat = 1),
        c("lon", "lat"),
        1.5
      ),
      list("finite numbers", data.frame(lon = NA, lat = 1))
    )
  )
  s <- spobj(data.frame(lon = 1:3, lat = 1:3))
  errs("spxy", list(list("the object spobj\\(\\) returns", 1)))
  errs("crsof", list(list("the object spobj\\(\\) returns", 1)))
  errs(
    "lfspec",
    list(
      list("the object spobj\\(\\) returns", 1),
      list(
        "EPSG:4326 and EPSG:3857",
        spobj(data.frame(e = 1, n = 1), c("e", "n"), 2958)
      ),
      list(
        "`x\\$coords` must be a matrix",
        structure(list(coords = 1:3, crs = 4326), class = "spobj")
      )
    )
  )
  errs(
    "spsamp",
    list(
      list("`size` must be at most the number of points, 3", s, 4),
      list("`dens` must be TRUE or FALSE", s, 1, NA),
      list(
        "at least 2 points to sample by density",
        spobj(data.frame(lon = 1, lat = 1)),
        1,
        TRUE,
        1
      ),
      list(
        "the points all coincide",
        spobj(data.frame(lon = c(1, 1), lat = c(2, 2))),
        1,
        TRUE,
        1
      )
    )
  )
  expect_error(
    lfspec(structure(
      list(coords = cbind(c(NA, 1), c(1, 1)), crs = 3857),
      class = "spobj"
    )),
    "must be finite"
  )
})
