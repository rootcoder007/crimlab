# idwint(), gtwreg() and asrule() against the methods written out here in
# base R: inverse-distance weights over the nearest monitors, zonal means;
# lm() with kernel weights at each observation, its hat values, the fit
# without the observation, the AICc and the local variances; association
# measures from table(); and the refusals.

rhav <- function(x0, y0, x1, y1) {
  r <- pi / 180
  h <- sin((y1 - y0) * r / 2)^2 +
    cos(y0 * r) * cos(y1 * r) * sin((x1 - x0) * r / 2)^2
  2 * 6371.0088 * asin(sqrt(pmin(1, h)))
}

# IDW at one target, written out
ridw <- function(mx, my, mv, x, y, p, k, rad, ll) {
  ok <- !is.na(mv)
  d <- if (ll) rhav(mx, my, x, y) else sqrt((mx - x)^2 + (my - y)^2)
  if (any(ok & d == 0)) return(mean(mv[ok & d == 0]))
  j <- which(ok & d <= rad)
  j <- j[order(d[j], j)]
  if (k > 0) j <- utils::head(j, k)
  if (!length(j)) return(NA_real_)
  w <- d[j]^-p
  sum(w * mv[j]) / sum(w)
}

test_that("idwint() is the weighted mean of the nearest monitors", {
  set.seed(51)
  m <- data.frame(x = runif(15, 0, 10), y = runif(15, 0, 10),
                  v = rnorm(15, 10))
  m$v[4] <- NA
  tg <- data.frame(x = c(runif(20, 0, 10), m$x[2], 50),
                   y = c(runif(20, 0, 10), m$y[2], 50))
  for (a in list(list(p = 2, k = 4, rad = 5), list(p = 1, k = 0, rad = Inf),
                 list(p = 0, k = 1, rad = 3))) {
    r <- idwint(m$x, m$y, m$v, tg$x, tg$y, p = a$p, k = a$k, rad = a$rad)
    want <- mapply(function(x, y) {
      ridw(m$x, m$y, m$v, x, y, a$p, a$k, a$rad, FALSE)
    }, tg$x, tg$y)
    expect_equal(r$val, unname(want), tolerance = 1e-13, info = a$p)
  }
  r <- idwint(m$x, m$y, m$v, tg$x, tg$y, k = 4, rad = 5)
  expect_identical(r$val[21], m$v[2])
  expect_identical(r$val[22], NA_real_)
  expect_identical(r$used[22], 0)
  # on the sphere, Rahman et al.'s rule: d^-2, four monitors within 50 km
  ll <- idwint(m$x / 10 - 118, m$y / 10 + 34, m$v, tg$x / 10 - 118,
               tg$y / 10 + 34, p = 2, k = 4, rad = 50, lonlat = TRUE)
  want <- mapply(function(x, y) {
    ridw(m$x / 10 - 118, m$y / 10 + 34, m$v, x, y, 2, 4, 50, TRUE)
  }, tg$x / 10 - 118, tg$y / 10 + 34)
  expect_equal(ll$val, unname(want), tolerance = 1e-13)
  # zonal means over the targets with a value
  zn <- rep(c("b", "a"), 11)
  z <- idwint(m$x, m$y, m$v, tg$x, tg$y, k = 4, rad = 5, zone = zn)
  expect_identical(z$zone, c("a", "b"))
  expect_equal(z$val, as.vector(tapply(r$val, zn, mean, na.rm = TRUE)),
               tolerance = 1e-14)
  expect_identical(z$n, as.numeric(tapply(!is.na(r$val), zn, sum)))
})

test_that("idwint() refusals name the function", {
  e <- function(m, ...) expect_error(idwint(...), m, fixed = TRUE)
  e("idwint: `mx`, `my` and `mv` must be numeric vectors of one length.",
    c(0, 1), c(0, 1), 1, 0, 0)
  e("idwint: `tx` and `ty` must have one length.", 0, 0, 1, c(0, 1), 0)
  e("idwint: `p` must be at least 0.", 0, 0, 1, 0, 0, p = -1)
  e("idwint: `k` must be a whole number.", 0, 0, 1, 0, 0, k = 1.5)
  e("idwint: `rad` must be a positive number or Inf.", 0, 0, 1, 0, 0,
    rad = -1)
  e("idwint: `zone` must give the zone of each target, with no NA.", 0, 0,
    1, c(0, 1), c(0, 1), zone = c(1, NA))
})

# the local fits of GTWR written out with lm(): coefficients, hat values,
# leave-one-out residuals and the variances of the coefficients
rgtw <- function(y, x, u, v, tt, h, tau) {
  n <- length(y)
  xm <- cbind(1, x)
  w <- function(i) {
    exp(-((u - u[i])^2 + (v - v[i])^2 + tau * (tt - tt[i])^2) / h^2)
  }
  loc <- lapply(seq_len(n), function(i) {
    wi <- w(i)
    f <- stats::lm(y ~ x, weights = wi)
    a <- solve(crossprod(xm, wi * xm))
    wo <- wi
    wo[i] <- 0
    fo <- stats::lm(y ~ x, weights = wo)
    list(b = stats::coef(f), s = drop(xm[i, , drop = FALSE] %*% a %*%
                                        t(wi * xm)),
         loo = y[i] - sum(stats::coef(fo) * xm[i, ]),
         vd = diag(a %*% crossprod(xm, wi^2 * xm) %*% a))
  })
  srow <- t(vapply(loc, `[[`, numeric(n), "s"))
  fit <- vapply(seq_len(n), function(i) sum(loc[[i]]$b * xm[i, ]), 0)
  rss <- sum((y - fit)^2)
  trs <- sum(diag(srow))
  s2 <- rss / (n - 2 * trs + sum(srow^2))
  list(b = t(vapply(loc, `[[`, numeric(2), "b")), fit = fit,
       hat = diag(srow), loo = vapply(loc, `[[`, 0, "loo"), rss = rss,
       trs = trs, s2 = s2,
       se = sqrt(s2 * t(vapply(loc, `[[`, numeric(2), "vd"))),
       aicc = n * log(rss / n) + n * log(2 * pi) +
         n * (n + trs) / (n - 2 - trs))
}

test_that("gtwreg() is lm() with kernel weights at each observation", {
  set.seed(52)
  d <- expand.grid(u = 1:4, v = 1:4, t = 1:3)
  d$x <- rnorm(nrow(d))
  d$y <- 1 + (0.5 + 0.2 * d$u) * d$x + rnorm(nrow(d), 0, 0.3)
  r <- gtwreg(d$y, d["x"], d$u, d$v, d$t, h = 2.5, tau = 0.5)
  g <- rgtw(d$y, d$x, d$u, d$v, d$t, 2.5, 0.5)
  expect_equal(unname(r$coef), unname(g$b), tolerance = 1e-10)
  expect_identical(colnames(r$coef), c("(const)", "x"))
  expect_equal(r$fit, g$fit, tolerance = 1e-10)
  expect_equal(r$hat, g$hat, tolerance = 1e-10)
  expect_equal(r$loo, g$loo, tolerance = 1e-9)
  expect_equal(r$stat[["cv"]], sum(g$loo^2), tolerance = 1e-9)
  expect_equal(r$stat[["trs"]], g$trs, tolerance = 1e-10)
  expect_equal(r$stat[["sigma2"]], g$s2, tolerance = 1e-10)
  expect_equal(unname(r$se), unname(g$se), tolerance = 1e-9)
  expect_equal(r$stat[["aicc"]], g$aicc, tolerance = 1e-10)
  expect_equal(r$resid, d$y - g$fit, tolerance = 1e-10)
  # a flat kernel is ordinary least squares at every observation
  f <- gtwreg(d$y, d["x"], d$u, d$v, h = 1e8)
  ols <- stats::coef(stats::lm(y ~ x, d))
  expect_equal(unname(f$coef[5, ]), unname(ols), tolerance = 1e-8)
  expect_equal(f$stat[["trs"]], 2, tolerance = 1e-6)
})

test_that("gtwreg() picks the pair with the smallest criterion", {
  set.seed(53)
  d <- expand.grid(u = 1:5, v = 1:5)
  d$x <- rnorm(25)
  d$y <- (1 + 0.3 * d$u) * d$x + rnorm(25, 0, 0.2)
  r <- gtwreg(d$y, d["x"], d$u, d$v, h = c(1, 2, 4, 50))
  cv <- vapply(c(1, 2, 4, 50), function(h) {
    sum(rgtw(d$y, d$x, d$u, d$v, 0 * d$u, h, 0)$loo^2)
  }, 0)
  expect_equal(r$grid$cv, cv, tolerance = 1e-8)
  expect_identical(r$stat[["h"]], c(1, 2, 4, 50)[which.min(cv)])
  a <- gtwreg(d$y, d["x"], d$u, d$v, h = c(1, 2, 4, 50), crit = "aicc")
  expect_identical(a$stat[["h"]], a$grid$h[which.min(a$grid$aicc)])
  expect_output(print(r), "geographically and temporally weighted")
  e <- function(m, ...) expect_error(gtwreg(...), m, fixed = TRUE)
  e("gtwreg: `h` must be positive.", d$y, d["x"], d$u, d$v, h = 0)
  e("gtwreg: `tau` must be at least 0.", d$y, d["x"], d$u, d$v, h = 1,
    tau = -1)
  e("gtwreg: `u` and `v` must have one value per value of `y`.", d$y,
    d["x"], 1, d$v, h = 1)
  e("gtwreg: the model needs more observations than coefficients plus one.",
    1:3, cbind(1:3), 1:3, 1:3, h = 1)
  e("gtwreg: `crit` must be one of", d$y, d["x"], d$u, d$v, h = 1,
    crit = "bic")
  e("gtwreg: a local fit is singular", d$y, cbind(z = rep(2, 25)), d$u,
    d$v, h = 1)
})

test_that("asrule() is the measures of table() written out", {
  set.seed(54)
  x <- factor(sample(c("vl", "l", "m"), 200, TRUE), c("vl", "l", "m"))
  y <- sample(c("p", "q"), 200, TRUE, prob = c(0.3, 0.7))
  r <- asrule(x, y, minsup = 0.05, minconf = 0.1, nm = c("pm", "ipv"))
  tb <- table(x, y) / 200
  px <- rowSums(tb)
  py <- colSums(tb)
  want <- do.call(rbind, c(
    lapply(names(px), function(a) {
      data.frame(lhs = "", rhs = paste0("pm=", a), support = px[[a]],
                 confidence = px[[a]], coverage = 1, lift = 1)
    }),
    lapply(names(py), function(b) {
      data.frame(lhs = "", rhs = paste0("ipv=", b), support = py[[b]],
                 confidence = py[[b]], coverage = 1, lift = 1)
    }),
    unlist(lapply(levels(x), function(a) {
      lapply(names(py), function(b) {
        s <- tb[a, b]
        rbind(data.frame(lhs = paste0("pm=", a), rhs = paste0("ipv=", b),
                         support = s, confidence = s / px[[a]],
                         coverage = px[[a]], lift = s / px[[a]] / py[[b]]),
              data.frame(lhs = paste0("ipv=", b), rhs = paste0("pm=", a),
                         support = s, confidence = s / py[[b]],
                         coverage = py[[b]], lift = s / py[[b]] / px[[a]]))
      })
    }), recursive = FALSE)))
  want <- want[want$support >= 0.05 & want$confidence >= 0.1, ]
  expect_identical(r$lhs, want$lhs)
  expect_identical(r$rhs, want$rhs)
  expect_equal(r$support, want$support, tolerance = 1e-14)
  expect_equal(r$confidence, want$confidence, tolerance = 1e-14)
  expect_equal(r$coverage, want$coverage, tolerance = 1e-14)
  expect_equal(r$lift, want$lift, tolerance = 1e-13)
  expect_equal(r$count, want$support * 200, tolerance = 1e-12)
  # without the empty left-hand side, and a pair never seen together
  q <- asrule(c("a", "a", "b"), c("p", "q", "q"), minsup = 0,
              minconf = 0, empty = FALSE)
  expect_identical(q$lhs, c("x=a", "y=p", "x=a", "y=q", "x=b", "y=q"))
  e <- function(m, ...) expect_error(asrule(...), m, fixed = TRUE)
  e("asrule: `x` and `y` must have one length and no NA.", 1:2, 1)
  e("asrule: `minsup` must be between 0 and 1.", 1, 1, minsup = 2)
  e("asrule: `nm` must be two names.", 1, 1, nm = "a")
})

test_that("gtwreg() takes x = NULL and refuses x of the wrong type", {
  v <- c(2.5, 1.3, 4.1, 3.3, 0.7, 5.2, 2.2, 1.9, 3.8, 4.4, 2.9, 1.1)
  g <- gtwreg(v, NULL, v, rev(v), h = 2)
  expect_identical(colnames(g$coef), "(const)")
  m <- "`x` must be finite numbers, one row per value of `y`"
  errs("gtwreg", list(list(m, v, function(x) x, v, v, h = 2),
                      list(m, v, globalenv(), v, v, h = 2),
                      list(m, v, list(1), v, v, h = 2)))
})
