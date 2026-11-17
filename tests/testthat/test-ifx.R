# A panel of 25 units over 10 years, the first five treated from year 7
# (three of them) or year 8, with one factor and one covariate.
ifdat <- function(seed = 71) {
  set.seed(seed)
  p <- expand.grid(id = 1:25, yr = 1:10)
  f <- cumsum(rnorm(10))
  l <- rnorm(25)
  p$x <- rnorm(nrow(p))
  a <- c(7, 7, 7, 8, 8, rep(Inf, 20))
  p$tr <- as.integer(p$yr >= a[p$id])
  p$y <- 0.5 * p$x + l[p$id] * f[p$yr] + p$id / 10 + 1.5 * p$tr +
    rnorm(nrow(p), 0, 0.2)
  p
}

# The counterfactuals written out: Bai's estimator on the controls with r
# factors (beta by alternating least squares on doubly demeaned data),
# each treated unit's intercept and loadings by lm() before adoption.
ifref <- function(p, r, x = TRUE) {
  co <- p[p$id > 5, ]
  co <- co[order(co$id, co$yr), ]
  nt <- 10
  dd <- function(m) {
    m <- matrix(m, nt)
    m - outer(rowMeans(m), colMeans(m), "+") + mean(m)
  }
  b <- 0
  if (x) b <- stats::coef(stats::lm(y ~ x + factor(id) + factor(yr), co))[[2]]
  for (i in 1:5000) {
    w <- dd(co$y - b * co$x)
    ev <- eigen(w %*% t(w), symmetric = TRUE)$vectors
    ff <- sqrt(nt) * ev[, seq_len(r), drop = FALSE]
    cm <- ff %*% t(ff) %*% w / nt
    if (!x) break
    b1 <- sum(dd(co$x) * (dd(co$y) - cm)) / sum(dd(co$x)^2)
    if (abs(b1 - b) < 1e-12) break
    b <- b1
  }
  wv <- matrix(co$y - b * co$x, nt)
  mu <- mean(wv)
  xi <- rowMeans(wv) - mu
  out <- NULL
  for (u in 1:5) {
    s <- p[p$id == u, ]
    s <- s[order(s$yr), ]
    a <- min(s$yr[s$tr == 1])
    z <- s$y - b * s$x - mu - xi
    pre <- s$yr < a
    g <- stats::lm.fit(cbind(1, ff)[pre, , drop = FALSE], z[pre])$coefficients
    y0 <- b * s$x + mu + xi + drop(cbind(1, ff) %*% g)
    out <- rbind(out, data.frame(unit = u, yr = s$yr, e = s$yr - a,
                                 y0 = y0, eff = s$y - y0))
  }
  list(b = b, cf = out)
}

test_that("ifesc() with no factors is the two-way model carried over", {
  p <- ifdat()
  r <- ifesc(p, "y", "tr", "id", "yr", x = "x", r = 0, se = "none")
  ref <- ifref(p, 0)
  expect_equal(unname(r$beta), unname(ref$b), tolerance = 1e-10)
  expect_equal(r$cf$y0, ref$cf$y0, tolerance = 1e-10)
  expect_equal(r$att$est, mean(ref$cf$eff[ref$cf$e >= 0]), tolerance = 1e-10)
  expect_identical(r$cf$unit, as.character(ref$cf$unit))
  expect_identical(r$r, 0L)
  expect_null(r$cv)
  expect_null(r$draws)
  expect_true(is.na(r$att$se))
})

test_that("ifesc() with one factor is Bai's estimator and lm()", {
  p <- ifdat()
  r <- ifesc(p, "y", "tr", "id", "yr", x = "x", r = 1, se = "none",
             tol = 1e-12)
  ref <- ifref(p, 1)
  expect_true(r$conv)
  expect_equal(unname(r$beta), unname(ref$b), tolerance = 1e-8)
  expect_equal(r$cf$y0, ref$cf$y0, tolerance = 1e-8)
  ev <- tapply(ref$cf$eff, ref$cf$e, mean)
  expect_equal(r$eff$est, as.vector(ev), tolerance = 1e-8)
  expect_identical(r$eff$e, as.numeric(names(ev)))
  expect_identical(r$eff$n, as.vector(table(ref$cf$e)))
  expect_equal(r$att$est, mean(ref$cf$eff[ref$cf$e >= 0]), tolerance = 1e-8)
  # the factor spans the eigenvector's direction, its sign as eigsym()
  # fixes it
  expect_identical(dim(r$factors), c(10L, 1L))
  expect_equal(sum(r$factors^2), 10, tolerance = 1e-10)
})

test_that("ifesc(r = NULL) picks r by leaving out pre-adoption periods", {
  p <- ifdat()
  r <- ifesc(p, "y", "tr", "id", "yr", rmax = 3, se = "none")
  # without covariates the controls' fit has no iteration; each treated
  # unit's period s is predicted from the others
  ms <- vapply(0:3, function(k) {
    co <- p[p$id > 5, ]
    co <- co[order(co$id, co$yr), ]
    w <- matrix(co$y, 10)
    mu <- mean(w)
    xi <- rowMeans(w) - mu
    w <- w - outer(rowMeans(w), colMeans(w), "+") + mu
    ev <- eigen(w %*% t(w), symmetric = TRUE)$vectors
    ff <- sqrt(10) * ev[, seq_len(k), drop = FALSE]
    e <- NULL
    for (u in 1:5) {
      s <- p[p$id == u, ]
      s <- s[order(s$yr), ]
      a <- min(s$yr[s$tr == 1])
      z <- s$y - mu - xi
      for (t in seq_len(a - 1)) {
        k1 <- setdiff(seq_len(a - 1), t)
        xx <- cbind(1, ff)
        g <- stats::lm.fit(xx[k1, , drop = FALSE], z[k1])$coefficients
        e <- c(e, z[t] - sum(xx[t, ] * g))
      }
    }
    mean(e^2)
  }, 0)
  expect_equal(r$cv$mspe, ms, tolerance = 1e-8)
  expect_identical(r$cv$r, 0:3)
  expect_identical(r$r, which.min(ms) - 1L)
})

test_that("ifesc()'s bootstrap draws units with R's generator", {
  p <- ifdat()
  a <- ifesc(p, "y", "tr", "id", "yr", r = 1, nboot = 20, seed = 4)
  b <- ifesc(p, "y", "tr", "id", "yr", r = 1, nboot = 20, seed = 4)
  expect_identical(a$draws, b$draws)
  expect_identical(dim(a$draws), c(20L, 1L + nrow(a$eff)))
  # a draw that could not be fitted is NA and left out
  expect_equal(a$att$se, stats::sd(a$draws[, 1L], na.rm = TRUE),
               tolerance = 1e-12)
  expect_identical(a$nboot, 20L)
})

test_that("ifesc() refuses what it cannot fit", {
  p <- ifdat()
  expect_error(ifesc(p[-1, ], "y", "tr", "id", "yr", r = 0),
               paste("ifesc: the panel must be balanced: each unit once at",
                     "each time."), fixed = TRUE)
  q <- p
  q$tr[q$id == 1 & q$yr == 9] <- 0L
  expect_error(ifesc(q, "y", "tr", "id", "yr", r = 0),
               "ifesc: the treatment must stay 1 once it turns 1.",
               fixed = TRUE)
  q <- p
  q$tr[q$id == 1] <- 1L
  expect_error(ifesc(q, "y", "tr", "id", "yr", r = 0),
               paste("ifesc: a treated unit is treated at the first time:",
                     "it has no periods before adoption."), fixed = TRUE)
  q <- p
  q$tr <- 0L
  expect_error(ifesc(q, "y", "tr", "id", "yr", r = 0),
               "ifesc: no unit is treated.", fixed = TRUE)
  q <- p
  q$tr[q$yr == 10] <- 1L
  expect_error(ifesc(q, "y", "tr", "id", "yr", r = 0),
               "ifesc: the panel needs two or more never-treated units.",
               fixed = TRUE)
  q <- p
  q$tr <- q$tr * 2
  expect_error(ifesc(q, "y", "tr", "id", "yr", r = 0),
               "ifesc: column `tr` must hold 0 or 1.", fixed = TRUE)
  expect_error(ifesc(p, "y", "tr", "id", "yr", r = 6, se = "none"),
               paste("ifesc: a treated unit has too few periods before",
                     "adoption for the factors, or a fit is singular."),
               fixed = TRUE)
})

test_that("ifesc() prints, and leaves out the draws it cannot fit", {
  p <- ifdat()
  r <- ifesc(p, "y", "tr", "id", "yr", x = "x", r = 1, se = "none",
             maxit = 1)
  # five treated units, twenty never treated, ten years; one alternation
  # of beta and the factors cannot meet the tolerance
  expect_false(r$conv)
  expect_output(print(r), paste("ifesc: generalized synthetic control, 5",
                                "treated and 20 control units, 10 times, 1",
                                "factor(s)\n"), fixed = TRUE)
  expect_output(print(r), "Bai's estimator did not converge.", fixed = TRUE)
  expect_identical(summary(r), r$eff)
  cv <- ifesc(p, "y", "tr", "id", "yr", rmax = 1, se = "none")
  expect_output(print(cv), "factor(s) (cross-validated)", fixed = TRUE)
  q <- p
  q$yr <- as.character(q$yr)
  expect_error(ifesc(q, "y", "tr", "id", "yr", r = 0),
               "ifesc: column `yr` must hold numbers or dates.", fixed = TRUE)
  # one treated unit and three controls: a draw can be fitted when it
  # holds the treated unit and two or more of the controls
  set.seed(5)
  s <- expand.grid(id = 1:4, yr = 1:6)
  s$tr <- as.integer(s$id == 1 & s$yr >= 4)
  s$y <- stats::rnorm(nrow(s)) + s$tr
  bad <- function(seed) {
    bm <- .ifdraws(4L, 3L, seed)
    sum(bm[1L, ] == 0L | colSums(bm[-1L, ] > 0L) < 2L)
  }
  expect_identical(bad(1), 1L)
  expect_warning(
    b <- ifesc(s, "y", "tr", "id", "yr", r = 0, nboot = 3, seed = 1),
    "ifesc: 1 of 3 bootstrap draws could not be fitted and are left out.",
    fixed = TRUE
  )
  expect_identical(b$nfail, 1L)
  expect_identical(bad(3), 2L)
  expect_error(
    suppressWarnings(ifesc(s, "y", "tr", "id", "yr", r = 0, nboot = 3,
                           seed = 3)),
    "ifesc: fewer than two bootstrap draws could be fitted.", fixed = TRUE
  )
})

test_that("ifesc() refuses counts past R's integers in words", {
  p <- ifdat()
  errs("ifesc", list(
    list("`maxit` must be between 1 and 2147483647", p, "y", "tr", "id",
         "yr", maxit = 2^53),
    list("`r` must be between 0 and 2147483647", p, "y", "tr", "id", "yr",
         r = 1e308)
  ))
})
