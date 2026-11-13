# cbatt(), smdbal() and evlx(): IPW, entropy balancing and matching
# against stored values of WeightIt 1.7.0 and MatchIt 4.7.2, balance
# against cobalt 4.6.3, E-values against EValue 4.1.4 (dev/ref-cba.R);
# g-computation and subclassification written out with lm() and cut();
# recovery of a known effect and the bootstrap error against the Monte
# Carlo spread; and the refusals.

z <- cbdata()

test_that("Hajek IPW and entropy balancing are WeightIt's", {
  want <- c(ate = 2.04916195169231, att = 2.15879019509982,
            atc = 1.92247223977355)
  for (t in names(want)) {
    r <- cbatt(z$y, z$d, z$x, tgt = t, reps = 0)
    expect_equal(r$res$est, unname(want[t]), tolerance = 1e-8, info = t)
  }
  e <- cbatt(z$y, z$d, z$x, est = "ebal", tgt = "att", reps = 0)
  expect_equal(e$res$est, 2.05766808929619, tolerance = 1e-8)
  # the weighted control means equal the treated means
  w <- e$w
  for (v in z$x) {
    expect_equal(weighted.mean(v[z$d == 0], w[z$d == 0]),
                 mean(v[z$d == 1]), tolerance = 1e-9)
  }
  expect_equal(sum(w[z$d == 0]), sum(z$d), tolerance = 1e-12)
  expect_true(all(is.na(e$ps)))
  # given propensities are used as they are
  p <- glm(z$d ~ z$x$a + z$x$b, family = binomial)$fitted.values
  g <- cbatt(z$y, z$d, z$x, ps = p, reps = 0)
  expect_equal(g$res$est, 2.04916195169231, tolerance = 1e-8)
  expect_equal(g$w[z$d == 1], 1 / p[z$d == 1], ignore_attr = TRUE)
})

test_that("matching is MatchIt's greedy caliper matching", {
  r <- cbatt(z$y, z$d, z$x, est = "psm", tgt = "att", calip = 0.2)
  expect_equal(r$res$est, 2.26923311008552, tolerance = 1e-10)
  expect_identical(2L * nrow(r$match), 240L)
  dif <- z$y[r$match$focal] - z$y[r$match$match]
  expect_equal(r$res$med, median(dif), tolerance = testthat_tolerance())
  expect_equal(r$res$se, sd(dif) / sqrt(length(dif)),
               tolerance = testthat_tolerance())
  # the ATC matches each control to treated units
  a <- cbatt(z$y, z$d, z$x, est = "psm", tgt = "atc", k = 2)
  expect_true(all(z$d[a$match$focal] == 0 & z$d[a$match$match] == 1))
  ft <- split(a$match$match, a$match$focal)
  dc <- vapply(names(ft), function(i) mean(z$y[ft[[i]]]) - z$y[as.integer(i)],
               0)
  expect_equal(a$res$est, mean(dc), tolerance = 1e-12)
  expect_output(print(a), "psm, ATC")
})

test_that("g-computation and subclassification are as written out", {
  d <- cbind(z$x, y = z$y)
  m1 <- lm(y ~ a + b, d[z$d == 1, ])
  m0 <- lm(y ~ a + b, d[z$d == 0, ])
  tau <- predict(m1, d) - predict(m0, d)
  for (t in c("ate", "att", "atc")) {
    k <- switch(t, ate = TRUE, att = z$d == 1, atc = z$d == 0)
    expect_equal(cbatt(z$y, z$d, z$x, est = "gcomp", tgt = t,
                       reps = 0)$res$est, mean(tau[k]), tolerance = 1e-10)
  }
  e <- glm(z$d ~ z$x$a + z$x$b, family = binomial)$fitted.values
  for (t in c("ate", "att")) {
    q <- if (t == "ate") e else e[z$d == 1]
    s <- findInterval(e, quantile(q, (1:4) / 5), left.open = TRUE)
    m <- tapply(z$y[z$d == 1], s[z$d == 1], mean) -
      tapply(z$y[z$d == 0], s[z$d == 0], mean)
    w <- if (t == "ate") table(s) else table(s[z$d == 1])
    expect_equal(cbatt(z$y, z$d, z$x, est = "sub", tgt = t,
                       reps = 0)$res$est, sum(m * w) / sum(w),
                 tolerance = 1e-8)
  }
})

test_that("the effect is recovered and the bootstrap error fits", {
  est <- se <- numeric(100)
  for (k in seq_along(est)) {
    s <- cbdata(100 + k, n = 200)
    r <- cbatt(s$y, s$d, s$x, reps = 99, seed = k)
    est[k] <- r$res$est
    se[k] <- r$res$se
  }
  expect_lt(abs(mean(est) - 2), 4 * sd(est) / sqrt(100))
  # 100 draws: the spread's own error is about 7%
  expect_lt(abs(mean(se) / sd(est) - 1), 0.2)
  a <- cbatt(z$y, z$d, z$x, est = "gcomp", reps = 30, seed = 2)
  b <- cbatt(z$y, z$d, z$x, est = "gcomp", reps = 30, seed = 2)
  expect_identical(a, b)
  set.seed(4)
  c1 <- cbatt(z$y, z$d, z$x, est = "sub", reps = 20)
  set.seed(4)
  c2 <- cbatt(z$y, z$d, z$x, est = "sub", reps = 20)
  expect_identical(c1, c2)
  expect_equal(c1$res$hi, c1$res$est + qnorm(0.975) * c1$res$se,
               tolerance = testthat_tolerance())
  e1 <- cbatt(z$y, z$d, z$x, est = "ebal", tgt = "att", reps = 20,
              seed = 1)
  expect_gt(e1$res$se, 0)
})

test_that("balance and overlap are cobalt's", {
  w <- cbatt(z$y, z$d, z$x, reps = 0)$w
  b <- smdbal(z$x, z$d, w = w)
  expect_equal(b$tab$smd, c(0.415415772942186, -0.460225491755417),
               tolerance = 1e-12)
  expect_equal(b$tab$vr, c(1.10227311997676, 1.24202198942649),
               tolerance = 1e-12)
  expect_equal(b$tab$smdw, c(0.0200126005626581, -0.0241650126413817),
               tolerance = 1e-10)
  expect_equal(b$tab$vrw, c(1.0586267943418, 1.28357498398007),
               tolerance = 1e-10)
  expect_true(b$bal)
  r <- smdbal(as.matrix(z$x), z$d)
  expect_false(r$bal)
  expect_identical(r$tab$var, c("a", "b"))
  e <- glm(z$d ~ z$x$a + z$x$b, family = binomial)$fitted.values
  o <- r$overlap
  expect_equal(o$min1, min(e[z$d == 1]), tolerance = 1e-9)
  lo <- max(min(e[z$d == 1]), min(e[z$d == 0]))
  hi <- min(max(e[z$d == 1]), max(e[z$d == 0]))
  expect_equal(o$out, mean(e < lo | e > hi), tolerance = 1e-12)
  expect_equal(o$wmax, max(ifelse(z$d == 1, 1 / e, 1 / (1 - e))),
               tolerance = 1e-8)
  expect_identical(r$ovl, o$out == 0)
  s <- smdbal(z$x$a, z$d, ps = rep(0.5, 300))
  expect_true(s$ovl)
  expect_identical(s$tab$var, "x1")
  expect_output(print(s), "smdbal: balanced")
})

test_that("E-values are EValue's", {
  expect_equal(evlx(3.9, 1.8, 8.7)[2:3], c(7.26303434416005, 3),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(evlx(0.8, 0.71, 0.91, type = "hr")[2:3],
               c(1.60889134837784, 1.3360826209926), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(evlx(2.5, 1.5, 4, type = "or")[2:3],
               c(2.53971129469801, 1.74939249466688), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(evlx(0.5, type = "smd", sd = 0.1)[2:3],
               c(2.52914198186127, 1.96803734698826), tolerance = 1e-12,
               ignore_attr = TRUE)
  # a rare outcome reads ratios as risk ratios; an interval over the null
  expect_identical(evlx(2, type = "or", rare = TRUE)[["rr"]], 2)
  expect_identical(evlx(2, 0.9, 4)[["evlim"]], 1)
  expect_identical(evlx(0.5, 0.3, 1.2)[["evlim"]], 1)
  expect_true(is.na(evlx(2)[["evlim"]]))
  expect_named(evlx(2), c("rr", "evalue", "evlim"))
})

test_that("cbatt(), smdbal() and evlx() refuse what they cannot do", {
  expect_error(cbatt(z$y, z$d, z$x, est = "ebal"), "ATT only")
  expect_error(cbatt(z$y, z$d, z$x, est = "psm"), "ATT or the ATC")
  expect_error(cbatt(z$y, c(1, rep(0, 299)), z$x), "two treated")
  expect_error(cbatt(z$y, z$d, z$x, ps = rep(1, 300)), "strictly between")
  sep <- data.frame(s = z$d)
  expect_error(cbatt(z$y, z$d, sep), "separates the arms")
  expect_error(smdbal(sep, z$d), "separates the arms")
  expect_error(cbatt(z$y, z$d, z$x, est = "psm", tgt = "att",
                     calip = 1e-9), "fewer than two")
  expect_error(cbatt(z$y, z$d, cbind(z$x$a, z$x$a), est = "gcomp"),
               "undefined")
  expect_error(cbatt(z$y, z$d, z$x, est = "sub", nsub = 150), "undefined")
  expect_error(cbatt(z$y, z$d, cbind(z$x$a, 2 * z$x$a), est = "ebal",
                     tgt = "att"), "undefined")
  expect_error(cbatt(z$y, z$d, matrix(1, 2, 2)), "one row per unit")
  expect_error(cbatt(z$y, z$d, matrix(NA_real_, 300, 2)), "finite numbers")
  expect_error(cbatt(z$y, z$d, list()), "one covariate")
  expect_error(cbatt(z$y, z$d, matrix(0, 300, 0)), "one covariate")
  expect_error(evlx(NULL), "one estimate")
  expect_error(evlx(2, 2.5, 3), "surround")
  expect_error(evlx(2, sd = 1), "standardised difference")
  expect_error(evlx(-1), "est")
})
