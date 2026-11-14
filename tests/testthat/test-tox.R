# txcal(), lodimp(), pmrrat(), etohlr() and mtxrel(): against weighted
# least squares, substitutions and rules done here in base R, and the
# refusals.

test_that("txcal() is the weighted least-squares line", {
  conc <- c(5, 10, 25, 50, 100, 250)
  resp <- c(0.051, 0.098, 0.26, 0.49, 1.02, 2.47)
  unk <- c(0.03, 0.6, NA, 0.004)
  for (w in c("1/x^2", "1/x", "none")) {
    wt <- switch(w, "1/x^2" = 1 / conc^2, "1/x" = 1 / conc, none = 1)
    g <- lm(resp ~ conc, weights = rep_len(wt, 6))
    s <- summary(g)
    b <- unname(coef(g))
    out <- txcal(conc, resp, w = w, unk = unk)
    expect_identical(out$w, w)
    expect_equal(unname(out$fit[c("int", "slope")]), b, tolerance = 1e-12)
    expect_equal(unname(out$fit["r2"]), s$r.squared, tolerance = 1e-12)
    expect_equal(unname(out$fit["sigma"]), s$sigma, tolerance = 1e-12)
    expect_equal(unname(out$fit["lod"]), 3.3 * s$sigma / b[2],
                 tolerance = 1e-12)
    expect_equal(unname(out$fit["loq"]), 10 * s$sigma / b[2],
                 tolerance = 1e-12)
    xh <- (unk - b[1]) / b[2]
    expect_equal(out$unk$conc, xh, tolerance = 1e-12)
    fl <- as.character(cut(xh, c(-Inf, c(3.3, 10) * s$sigma / b[2], Inf),
                           c("below_lod", "below_loq", "quantifiable"),
                           right = FALSE))
    expect_identical(as.character(out$unk$flag), fl)
    expect_identical(out$unk$resp, unk)
  }
  num <- txcal(conc, resp, w = 1 / conc, k = c(3, 9))
  g <- summary(lm(resp ~ conc, weights = 1 / conc))
  expect_identical(num$w, "numeric")
  expect_equal(unname(num$fit["lod"]), 3 * g$sigma / coef(g)[2, 1],
               tolerance = 1e-12)
  expect_null(num$unk)
})

test_that("lodimp() substitutes below the limit", {
  x <- c(0.4, 0.05, NA, 1.2, 0.08, 0.1)
  for (m in c("half", "sqrt2", "lod")) {
    f <- switch(m, half = 0.5, sqrt2 = 1 / sqrt(2), lod = 1)
    out <- lodimp(x, 0.1, m)
    cens <- !is.na(x) & x < 0.1
    expect_identical(out$tab$cens, cens)
    expect_identical(out$tab$imp, ifelse(cens, f * 0.1, x))
    expect_identical(out$frac, 2 / 5)
    expect_identical(out$nmiss, 1L)
  }
  l <- c(0.5, 0.01, 0.1, 0.1, 0.1, 0.2)
  out <- lodimp(x, l)
  expect_identical(out$tab$lod, l)
  expect_identical(out$tab$imp, ifelse(!is.na(x) & x < l, l / 2, x))
  expect_identical(lodimp(c(NA, NA), 1)$frac, NA_real_)
})

test_that("pmrrat() classes the central to peripheral ratio", {
  cen <- c(0.4, 1.5, 0.9, 2)
  per <- c(0.5, 1, 0.2, 1)
  out <- pmrrat(cen, per)
  expect_identical(out$ratio, cen / per)
  expect_identical(
    as.character(out$redist),
    c("minimal", "modest", "significant", "modest")
  )
  expect_identical(levels(out$redist), c("minimal", "modest", "significant"))
  o2 <- pmrrat(cen, per, cut = c(1.5, 4))
  expect_identical(
    as.character(o2$redist),
    c("minimal", "minimal", "significant", "modest")
  )
})

test_that("etohlr() reads ethyl glucuronide, then the higher alcohols", {
  eth <- c(0.8, 0.5, 0.3, 0.2, 0.6)
  prop <- c(0, 0.02, NA, NA, 0.01)
  but <- c(0.01, 0.03, 0.05, NA, 0)
  etg <- c(1.1, NA, 0, NA, NA)
  out <- etohlr(eth, prop, but, etg)
  hi <- suppressWarnings(pmax(prop, but, na.rm = TRUE))
  expect_identical(out$higher, hi)
  v <- vapply(seq_along(etg), function(i) {
    if (!is.na(etg[i]) && etg[i] > 0) return("antemortem")
    if (!is.na(hi[i]) && hi[i] > 0) "postmortem" else "indeterminate"
  }, "")
  expect_identical(as.character(out$verdict), v)
  o2 <- etohlr(eth, prop, but, etg, thr = 0.04)
  expect_identical(
    as.character(o2$verdict),
    c("antemortem", "indeterminate", "postmortem", "indeterminate",
      "indeterminate")
  )
  o3 <- etohlr(c(0.1, 0.2))
  expect_identical(o3$etg, c(NA_real_, NA_real_))
  expect_identical(as.character(o3$verdict), rep("indeterminate", 2))
})

test_that("mtxrel() ranks the matrices after the losses", {
  base <- c(vitreous = 0.95, liver = 0.8, femoral = 0.55, heart = 0.35,
            urine = 0.6)
  dil <- c(FALSE, FALSE, TRUE, TRUE, TRUE)
  out <- mtxrel(base, dil, days = 12, stage = 3)
  rel <- base * (1 - min(3 / 10, 0.6)) * ifelse(dil, 1 - min(12 / 30, 0.8), 1)
  o <- order(rel, decreasing = TRUE)
  expect_identical(out$matrix, names(base)[o])
  expect_equal(out$rel, unname(rel[o]), tolerance = 1e-15)
  expect_identical(out$base, unname(base[o]))
  expect_identical(out$dil, dil[o])
  expect_identical(out$rank, 1:5)
  u <- mtxrel(c(0.2, 0.9), c(TRUE, FALSE), days = 100, stage = 100,
              dmax = 10, dcap = 0.5, smax = 5, scap = 0.3)
  expect_identical(u$matrix, c("2", "1"))
  expect_equal(u$rel, c(0.9 * 0.7, 0.2 * 0.7 * 0.5), tolerance = 1e-15)
})

test_that("the toxicology functions refuse bad input", {
  errs("txcal", list(
    list("at least three calibrators", c(1, 2), c(1, 2)),
    list("`conc` must be greater than 0", c(0, 1, 2), 1:3),
    list("`w` must be one of", 1:3, 1:3, w = "1/y"),
    list("`w` must be greater than 0", 1:3, 1:3, w = c(1, 0, 1)),
    list("`k` must have length 2", 1:3, 1:3, k = 3),
    list("`unk` must be numbers or NA", 1:3, 1:3, unk = "a"),
    list("`unk` must be numbers or NA", 1:3, 1:3, unk = Inf),
    list("the calibrators give no line", 1:3, c(2, 2, 2))
  ))
  errs("lodimp", list(
    list("`lod` must be one limit, or one per value, 3", 1:3, c(1, 2)),
    list("`lod` must be greater than 0", 1:3, 0),
    list("`meth` must be one of", 1:3, 1, meth = "x"),
    list("`x` must be numbers or NA", character(0), 1)
  ))
  errs("pmrrat", list(
    list("`cut` must be two increasing ratios", 1, 1, cut = c(2, 1)),
    list("`per` must be greater than 0", 1, 0),
    list("`cen` must be finite", NA, 1)
  ))
  errs("etohlr", list(
    list("`prop` must be numbers or NA", c(1, 2), prop = c(1, 2, 3)),
    list("`thr` must be one number from 0", 1, thr = -1),
    list("`eth` must be at least 0", -1)
  ))
  errs("mtxrel", list(
    list("`dil` must be TRUE or FALSE", c(a = 1), 1),
    list("`days` must be one number from 0", 1, TRUE, days = -1),
    list("`dmax` must be one number strictly between 0", 1, TRUE, dmax = 0),
    list("`scap` must be one number from 0 to 1", 1, TRUE, scap = 2)
  ))
})
