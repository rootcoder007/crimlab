# Data description: the values of the reference implementation
# (helper-dsref.R: rmoriebricklayer 1.1.1 on these simulated inputs), and
# the refusals in crimlab's words.

cols <- function(x) unname(as.list(x))
utf8 <- function() {
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "not a UTF-8 session")
  skip_if_not(identical(Sys.getenv("CRIMLAB_ASCII_ONLY"), ""))
}

test_that("inlhist() draws the reference sketches", {
  utf8()
  set.seed(3)
  x <- round(rnorm(60), 3)
  expect_identical(inlhist(x), dsref$inl1)
  expect_identical(inlhist(c(x, 9), bins = 7), dsref$inl2)
  expect_identical(inlhist(rep(5, 3)), dsref$inl3)
  expect_identical(inlhist(c(1, 2, 2, 3, 3, 3), bins = 3), dsref$inl4)
  expect_identical(inlhist(NULL), strrep(" ", 10))
  expect_error(inlhist(c("a", "b")),
               "^inlhist: `x` must be numeric, got character\\.$")
  expect_error(inlhist(1:3, bins = 0),
               "^inlhist: `bins` must be a single positive integer\\.$")
})

test_that("freqtbl() counts as the reference does", {
  df <- data.frame(grade = c("a", "b", "b", "c", NA, "b"))
  f1 <- freqtbl(df, "grade")
  expect_s3_class(f1, "freqtbl")
  expect_identical(names(f1), c("value", "n", "pct", "pctvalid"))
  expect_equal(cols(f1), dsref$fq1)
  expect_type(f1$n, "integer")
  expect_equal(cols(freqtbl(c(1, 1, 2, 3, 3, 3))), dsref$fq2)
  expect_equal(cols(freqtbl(c("c", "a", "b", "a"), sort = FALSE)), dsref$fq3)
  set.seed(1)
  s <- sample(letters, 50, TRUE)
  expect_equal(cols(freqtbl(s, maxlvl = 5)), dsref$fq4)
  expect_identical(nrow(freqtbl(character(0))), 0L)
  expect_error(freqtbl(df),
               "^freqtbl: `col` is required when `data` is a data frame\\.$")
  expect_error(freqtbl(df, "zz"), "^freqtbl: no such column: zz\\.$")
  expect_error(freqtbl(1:3, maxlvl = 0), "`maxlvl` must be a single positive")
  expect_output(print(f1), "Frequency table")
})

cor_data <- function() {
  set.seed(2)
  d <- data.frame(a = round(rnorm(30), 2), b = round(rnorm(30), 2),
                  g = sample(letters[1:3], 30, TRUE))
  d$c <- round(d$a * 2 + rnorm(30, sd = 0.3), 2)
  d$e <- round(exp(d$a), 3)
  d$a[c(2, 5)] <- NA
  d$b[7] <- NA
  d
}

test_that("cortbl() and topcors() give the reference correlations", {
  d <- cor_data()
  c1 <- cortbl(d)
  expect_s3_class(c1, "cortbl")
  expect_identical(attr(c1, "method"), "spearman")
  expect_equal(cols(c1), dsref$ct1, tolerance = 1e-12)
  expect_equal(cols(cortbl(d, meth = "pearson", minprs = 29)), dsref$ct2,
               tolerance = 1e-12)
  expect_equal(cols(topcors(d)), dsref$tc1, tolerance = 1e-12)
  expect_equal(cols(topcors(d, n = 2, minabs = 0.5, meth = "pearson")),
               dsref$tc2, tolerance = 1e-12)
  expect_error(cortbl(data.frame(a = 1:5, b = letters[1:5])),
               "^cortbl: need at least two numeric columns to correlate\\.$")
  expect_error(cortbl(d, minprs = 1), "`minprs` must be a single integer")
  expect_error(topcors(d, minabs = 2),
               "^topcors: `minabs` must be a single value in \\[0, 1\\]\\.$")
  expect_error(topcors(d, meth = "kendall"), "`meth` must be one of")
  expect_output(print(topcors(d)), "Top correlations")
})

test_that("drpempty() and drpconst() drop what the reference drops", {
  de <- data.frame(keep = c(1, 2, NA), allna = c(NA, NA, NA),
                   constant = c(7, 7, 7), s = c("x", NA, "y"))
  e1 <- drpempty(de, quiet = TRUE)
  expect_equal(cols(e1), dsref$de1)
  expect_identical(attr(e1, "dropped"), "allna")
  expect_equal(cols(drpempty(data.frame(a = c(1, NA), b = c(2, NA)),
                             which = "rows", quiet = TRUE)), dsref$de2)
  c1 <- drpconst(de, quiet = TRUE)
  expect_equal(cols(c1), dsref$dc1)
  expect_identical(attr(c1, "dropped"), c("allna", "constant"))
  expect_equal(cols(drpconst(data.frame(v = c(NA, NA, 5)), naasval = TRUE,
                             quiet = TRUE)), dsref$dc2)
  expect_message(drpempty(de), "Empty rows and columns dropped")
  expect_error(drpempty(c(1, NA)),
               "^drpempty: `data` must be a data frame\\.$")
  expect_error(drpempty(de, which = "both"), "`which` must be one or more of")
})

test_that("clncolnms() cleans names as the reference does", {
  nm <- c("Total  Population (2021)", "% change", "Ville / City", "dup",
          "dup", "R&D #", "camelCaseName", "2021 count", "", "été")
  expect_identical(clncolnms(nm, quiet = TRUE), dsref$cn1)
  for (cs in c("lower_camel", "upper_camel", "screaming_snake", "none")) {
    expect_identical(
      clncolnms(c("Total Pop", "per capita", "Total Pop"), case = cs,
                sep = ".", quiet = TRUE),
      dsref[[cs]]
    )
  }
  df <- data.frame(`Total Pop` = 1:2, `% change` = 3:4, check.names = FALSE)
  cl <- clncolnms(df, quiet = TRUE)
  expect_identical(names(cl), c("total_pop", "pct_change"))
  expect_identical(attr(cl, "orignms"), c("Total Pop", "% change"))
  expect_error(clncolnms("a b", case = "kebab"), "`case` must be one of")
  expect_error(clncolnms(c("a", NA)), "without NA")
})

maha_data <- function() {
  set.seed(4)
  m <- data.frame(h = round(rnorm(25, 170, 10), 1))
  m$w <- round(m$h * 0.5 + rnorm(25, 0, 5), 1)
  m$z <- round(runif(25), 3)
  m[1, 1:2] <- c(155, 100)
  m$w[4] <- NA
  m
}

test_that("mahaoutl() gives the reference distances", {
  m <- maha_data()
  r1 <- mahaoutl(m)
  expect_s3_class(r1, "mahaoutl")
  expect_identical(names(r1),
                   c("row", "distance", "p_value", "logpvalue", "outlier"))
  expect_equal(cols(r1), dsref$mh1, tolerance = 1e-10)
  expect_equal(cols(mahaoutl(m, robust = FALSE, alpha = 0.05)), dsref$mh2,
               tolerance = 1e-10)
  dup <- m
  dup$h2 <- dup$h
  expect_error(mahaoutl(dup), "^mahaoutl: the columns are collinear")
  expect_error(mahaoutl(m, alpha = 1), "strictly inside \\(0, 1\\)")
  expect_error(mahaoutl(m[1:3, ]), "need more complete rows than columns")
  expect_output(print(r1), "Mahalanobis outliers \\(robust\\)")
})

test_that("missruns(), missgrid(), misspat() and misssmry() match", {
  mr <- data.frame(outage = c(1, 2, NA, NA, NA, NA, 7, 8),
                   scattered = c(1, NA, 3, 4, NA, 6, NA, NA),
                   b = c(NA, NA, NA, NA, 1, 2, 3, 4))
  expect_equal(cols(missruns(mr)), dsref$mr1)
  expect_equal(cols(missruns(mr, minrun = 1)), dsref$mr2)
  expect_identical(nrow(missruns(data.frame(x = 1:5))), 0L)
  expect_error(missruns(mr, minrun = 0), "`minrun` must be a single positive")
  mp <- data.frame(id = 1:10, a = c(rep(NA, 3), 4:10), b = c(1:7, rep(NA, 3)),
                   cc = c(NA, 2:9, NA))
  p1 <- misspat(mp)
  expect_equal(cols(p1), dsref$mp1)
  expect_identical(attr(p1, "columns"), c("id", "a", "b", "cc"))
  expect_equal(cols(misspat(mp, maxpat = 2)), dsref$mp2)
  expect_error(misspat(data.frame()), "^misspat: `data` has no columns\\.$")
  ms <- data.frame(a = c(1, NA, 3), b = c(NA, NA, 3), c = 1:3)
  expect_equal(unname(misssmry(ms)), dsref$ms1)
  expect_identical(
    names(misssmry(ms)),
    c("n_rows", "n_cols", "n_missing", "pct_missing", "n_complete_rows",
      "pct_complete_rows", "n_cols_any_missing", "n_cols_all_missing")
  )
  expect_equal(unname(misssmry(data.frame(x = 1:4, note = NA))), dsref$ms2)
  expect_error(misssmry(c(1, NA)), "^misssmry: `data` must be a data frame")
  utf8()
  set.seed(1)
  g <- data.frame(complete = 1:37, block = c(rep(NA, 11), 12:37),
                  scattered = ifelse(runif(37) < 0.3, NA, 1),
                  mostlygone = c(1:5, rep(NA, 32)))
  expect_identical(capture.output(r <- missgrid(g, height = 6, width = 5)),
                   dsref$mg1)
  expect_identical(r, dsref$mg1r)
  expect_error(missgrid(g, height = 0), "`height` must be a single positive")
})

test_that("mcartst() gives the reference statistics", {
  set.seed(5)
  n <- 80
  xx <- round(rnorm(n), 3)
  yy <- round(xx + rnorm(n), 3)
  zz <- round(rnorm(n), 3)
  mc <- data.frame(x = xx, y = yy, z = zz)
  mc$y[xx > 0.6] <- NA
  mc$z[c(3, 9, 30)] <- NA
  r <- mcartst(mc)
  expect_s3_class(r, "mcartst")
  ref <- dsref$mc1
  nm <- c("statistic", "df", "logpvalue", "p_value", "npat", "nused",
          "nvars", "iterations", "mu", "sigma", "note")
  expect_identical(names(r), nm)
  expect_equal(unname(unlist(r[1:8])), unlist(ref[1:8]), tolerance = 1e-8)
  expect_equal(unname(r$mu), ref[[9L]], tolerance = 1e-8)
  expect_equal(as.vector(r$sigma), ref[[10L]], tolerance = 1e-8)
  expect_identical(names(r$mu), c("x", "y", "z"))
  mc2 <- data.frame(x = xx, y = yy)
  mc2$y[sample(n, 20)] <- NA
  r2 <- mcartst(mc2)
  expect_equal(unname(unlist(r2[1:8])), unlist(dsref$mc2[1:8]),
               tolerance = 1e-8)
  r3 <- mcartst(data.frame(a = xx, b = yy))
  expect_equal(r3$statistic, 0, tolerance = 1e-12)
  expect_identical(r3$df, 0)
  expect_identical(r3$note, dsref$mc3[[11L]])
  dup <- mc
  dup$x2 <- dup$x
  expect_error(mcartst(dup), "^mcartst: the columns are collinear")
  expect_warning(mcartst(cbind(mc, g = "a")), "dropping non-numeric column")
  expect_error(mcartst(mc, maxit = 0), "`maxit` must be one whole number")
  expect_output(print(r), "Little's MCAR test")
})

test_that("duprows() groups the reference duplicates", {
  dr <- data.frame(id = c(1, 2, 2, 3, 3, 3, NA, NA),
                   value = c("a", "b", "b", "c", "d", "c", "e", "e"))
  d1 <- duprows(dr, quiet = TRUE)
  expect_identical(names(d1), c("id", "value", "dupe_count"))
  expect_equal(cols(d1), dsref$dr1)
  expect_equal(cols(duprows(dr, "id", quiet = TRUE)), dsref$dr2)
  expect_identical(nrow(duprows(data.frame(x = 1:3), quiet = TRUE)), 0L)
  expect_error(duprows(dr, "code"), "^duprows: no such column\\(s\\): code\\.$")
  expect_message(duprows(dr), "Duplicated rows")
})

test_that("profcols() profiles as the reference does", {
  utf8()
  set.seed(9)
  pc <- data.frame(
    clean = round(rnorm(40), 3), skewed = c(round(rnorm(39), 3), 500),
    grade = sample(c("a", "b", "c"), 40, TRUE),
    gappy = c(rep(NA, 5), round(runif(35), 3)),
    inf = c(Inf, -2, 0, 0, 5, rep(1, 35)), allna = NA_real_
  )
  p <- profcols(pc)
  expect_s3_class(p, "profcols")
  expect_identical(names(p), c(dsref$pc2[1:19], "q25", "q50", "q75"))
  expect_equal(cols(p), dsref$pc1, tolerance = 1e-12)
  expect_identical(names(profcols(pc, qnts = c(0.025, 0.5))), dsref$pc2)
  expect_error(profcols(pc[, 0]), "has no columns to profile")
  expect_error(profcols(pc, qnts = 2), "must be probabilities")
  expect_output(print(p), "Column profile")
})

test_that("dstntst() gives the reference test", {
  own <- data.frame(group = c("Aboriginal", "Torres Strait Islander", "Both"),
                    est = c(24.0, 31.0, 24.6), lo = c(22.1, 26.2, 22.8),
                    hi = c(25.9, 35.8, 26.4))
  d1 <- dstntst(own, "group", "est", "lo", "hi", pooled = "Both")
  expect_s3_class(d1, "dstntst")
  expect_identical(
    names(d1),
    c("groups", "pooled", "poose", "q", "df", "p_value", "i2", "leftout",
      "conflvl")
  )
  ref <- dsref$dt1
  expect_equal(cols(d1$groups), ref[[1L]], tolerance = 1e-12)
  expect_equal(unlist(d1[2:7]), unlist(ref[2:7]), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(d1$leftout, "Both")
  div <- data.frame(
    division = c("North", "Central", "South", "Total, all divisions"),
    rate = c(4.1, 5.3, 4.4, 4.6), lo = c(3.6, 4.7, 3.9, 4.3),
    hi = c(4.6, 5.9, 4.9, 4.9), q = c("", "E", "ok", NA)
  )
  d2 <- dstntst(div, "division", "rate", "lo", "hi", qual = "q",
                conflvl = 0.9)
  ref <- dsref$dt2
  expect_equal(cols(d2$groups), ref[[1L]], tolerance = 1e-12)
  expect_equal(unlist(d2[2:7]), unlist(ref[2:7]), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(d2$leftout, "Total, all divisions")
  expect_error(dstntst(own, "group", "est", "lo", "hi", pooled = "None"),
               "`pooled` names \"None\", not a group of `data`")
  expect_error(dstntst(own[1, ], "group", "est", "lo", "hi"),
               "needs two groups at least")
  expect_error(dstntst(own, "grp", "est", "lo", "hi"),
               "^dstntst: column `grp` not found in `data`\\.$")
  expect_output(print(d1), "Cochran's Q")
})
