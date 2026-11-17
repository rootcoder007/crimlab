# Data description: the values of the reference implementation
# (helper-dsref.R: rmoriebricklayer 1.1.1 on these simulated inputs), and
# the refusals in crimlab's words.

cols <- function(x) {
  x <- unname(as.list(x))
  attributes(x) <- NULL
  x
}
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
  df <- data.frame(grade = c("a", "b", "b", "c", NA, "b"),
                   stringsAsFactors = FALSE)
  f1 <- freqtbl(df, "grade")
  expect_s3_class(f1, "freqtbl")
  expect_named(f1, c("value", "n", "pct", "pctvalid"))
  expect_equal(cols(f1), dsref$fq1, tolerance = testthat_tolerance())
  expect_type(f1$n, "integer")
  expect_equal(cols(freqtbl(c(1, 1, 2, 3, 3, 3))), dsref$fq2,
               tolerance = testthat_tolerance())
  expect_equal(cols(freqtbl(c("c", "a", "b", "a"), sort = FALSE)), dsref$fq3,
               tolerance = testthat_tolerance())
  set.seed(1)
  s <- sample(letters, 50, TRUE)
  expect_equal(cols(freqtbl(s, maxlvl = 5)), dsref$fq4,
               tolerance = testthat_tolerance())
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
                   constant = c(7, 7, 7), s = c("x", NA, "y"),
                   stringsAsFactors = FALSE)
  e1 <- drpempty(de, quiet = TRUE)
  expect_equal(cols(e1), dsref$de1, tolerance = testthat_tolerance())
  expect_identical(attr(e1, "dropped"), "allna")
  expect_equal(cols(drpempty(data.frame(a = c(1, NA), b = c(2, NA)),
                             which = "rows", quiet = TRUE)), dsref$de2,
               tolerance = testthat_tolerance())
  c1 <- drpconst(de, quiet = TRUE)
  expect_equal(cols(c1), dsref$dc1, tolerance = testthat_tolerance())
  expect_identical(attr(c1, "dropped"), c("allna", "constant"))
  expect_equal(cols(drpconst(data.frame(v = c(NA, NA, 5)), naasval = TRUE,
                             quiet = TRUE)), dsref$dc2,
               tolerance = testthat_tolerance())
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
  expect_named(cl, c("total_pop", "pct_change"))
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
  expect_named(r1,
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
  expect_equal(cols(missruns(mr)), dsref$mr1, tolerance = testthat_tolerance())
  expect_equal(cols(missruns(mr, minrun = 1)), dsref$mr2,
               tolerance = testthat_tolerance())
  expect_identical(nrow(missruns(data.frame(x = 1:5))), 0L)
  expect_error(missruns(mr, minrun = 0), "`minrun` must be a single positive")
  mp <- data.frame(id = 1:10, a = c(rep(NA, 3), 4:10), b = c(1:7, rep(NA, 3)),
                   cc = c(NA, 2:9, NA))
  p1 <- misspat(mp)
  expect_equal(cols(p1), dsref$mp1, tolerance = testthat_tolerance())
  expect_identical(attr(p1, "columns"), c("id", "a", "b", "cc"))
  expect_equal(cols(misspat(mp, maxpat = 2)), dsref$mp2,
               tolerance = testthat_tolerance())
  expect_error(misspat(data.frame()), "^misspat: `data` has no columns\\.$")
  ms <- data.frame(a = c(1, NA, 3), b = c(NA, NA, 3), c = 1:3)
  expect_equal(unname(misssmry(ms)), dsref$ms1,
               tolerance = testthat_tolerance())
  expect_named(
    misssmry(ms),
    c("n_rows", "n_cols", "n_missing", "pct_missing", "n_complete_rows",
      "pct_complete_rows", "n_cols_any_missing", "n_cols_all_missing")
  )
  expect_equal(unname(misssmry(data.frame(x = 1:4, note = NA))), dsref$ms2,
               tolerance = testthat_tolerance())
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
  expect_named(r, nm)
  expect_equal(unname(unlist(r[1:8])), unlist(ref[1:8]), tolerance = 1e-8)
  expect_equal(unname(r$mu), ref[[9L]], tolerance = 1e-8)
  expect_equal(as.vector(r$sigma), ref[[10L]], tolerance = 1e-8)
  expect_named(r$mu, c("x", "y", "z"))
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
                   value = c("a", "b", "b", "c", "d", "c", "e", "e"),
                   stringsAsFactors = FALSE)
  d1 <- duprows(dr, quiet = TRUE)
  expect_named(d1, c("id", "value", "dupe_count"))
  expect_equal(cols(d1), dsref$dr1, tolerance = testthat_tolerance())
  expect_equal(cols(duprows(dr, "id", quiet = TRUE)), dsref$dr2,
               tolerance = testthat_tolerance())
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
  expect_named(p, c(dsref$pc2[1:19], "q25", "q50", "q75"))
  expect_equal(cols(p), dsref$pc1, tolerance = 1e-12)
  expect_named(profcols(pc, qnts = c(0.025, 0.5)), dsref$pc2)
  expect_error(profcols(pc[, 0]), "has no columns to profile")
  expect_error(profcols(pc, qnts = 2), "must be probabilities")
  expect_output(print(p), "Column profile")
})

test_that("dstntst() gives the reference test", {
  own <- data.frame(group = c("Aboriginal", "Torres Strait Islander", "Both"),
                    est = c(24.0, 31.0, 24.6), lo = c(22.1, 26.2, 22.8),
                    hi = c(25.9, 35.8, 26.4), stringsAsFactors = FALSE)
  d1 <- dstntst(own, "group", "est", "lo", "hi", pooled = "Both")
  expect_s3_class(d1, "dstntst")
  expect_named(
    d1,
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
    hi = c(4.6, 5.9, 4.9, 4.9), q = c("", "E", "ok", NA),
    stringsAsFactors = FALSE
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

test_that("the shared checks refuse each kind of bad argument", {
  d <- data.frame(a = c(1, 2, 3, 5), b = c(2, 1, 4, 3))
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  nn <- d
  names(nn) <- c("a", "")
  ef(cortbl(nn), paste("cortbl: every column of `data` must have a name",
                       "(names(data) has an empty one)."))
  dd <- d
  names(dd) <- c("a", "a")
  ef(topcors(dd), "topcors: `data` has two columns named a: rename one.")
  bad <- list(factor("a"), d, mean, list(1), 1i)
  what <- c("a factor", "a data frame", "a function", "list", "complex")
  for (i in seq_along(bad)) {
    ef(topcors(d, minabs = bad[[i]]),
       sprintf("topcors: `minabs` must be a number, not %s.", what[[i]]))
  }
  ef(topcors(d, minabs = "abc"),
     "topcors: `minabs` must be numeric, not text such as \"abc\".")
  # text that reads as a number is taken as that number
  expect_identical(nrow(topcors(d, minabs = "0.5")), 1L)
  ef(topcors(d, minabs = TRUE), "topcors: `minabs` must be a number, not TRUE.")
  for (v in c(2.5, Inf)) {
    ef(cortbl(d, minprs = v),
       sprintf("cortbl: `minprs` must be a whole number, not %s.", v))
  }
  ef(cortbl(d, minprs = 3e9),
     "cortbl: `minprs` is too large: 3e+09 (the limit is 2147483647).")
  ef(freqtbl(mean), "freqtbl: `data` must be data, not a function.")
  ef(freqtbl(globalenv()), "freqtbl: `data` must be data, not an environment.")
  # NULL takes the first choice, rows: the empty row goes, the column stays
  e <- drpempty(data.frame(a = c(NA, 1), b = NA), which = NULL, quiet = TRUE)
  expect_identical(dim(e), c(1L, 2L))
  ef(drpempty(d, which = 1),
     "drpempty: `which` must be one or more of \"rows\", \"cols\", not 1.")
  for (fn in c("cortbl", "topcors", "drpconst", "missruns", "missgrid",
               "misspat", "profcols", "duprows", "dstntst")) {
    a <- if (fn == "dstntst") list(1:3, "g", "e", "l", "u") else list(1:3)
    ef(do.call(fn, a), sprintf("%s: `data` must be a data frame.", fn))
  }
  ef(topcors(data.frame(a = 1:3, b = letters[1:3])),
     "topcors: need at least two numeric columns to correlate.")
  for (fn in c("missruns", "missgrid", "misssmry")) {
    ef(do.call(fn, list(data.frame())),
       sprintf("%s: `data` has no columns.", fn))
  }
  ef(missgrid(d, width = 0),
     "missgrid: `width` must be one whole number from 1 to 10000.")
  ef(clncolnms("a", sep = 1), "clncolnms: `sep` must be a single string.")
  ef(mahaoutl(1:3),
     "mahaoutl: `data` must be a data frame or a numeric matrix.")
  ef(mahaoutl(data.frame(a = letters)),
     "mahaoutl: `data` has no numeric columns.")
  ef(mcartst(d, tol = 0), "mcartst: `tol` must be one number from 1e-12 to 1.")
  ef(mcartst(1:3), "mcartst: `data` must be a data frame or a numeric matrix.")
  ef(mcartst(data.frame(a = 1:3)),
     "mcartst: Little's test needs at least two numeric columns.")
  ef(mcartst(data.frame(a = c(1, Inf, 3), b = c(1, 2, NA))),
     paste0("mcartst: `data` has infinite values, which Little's test ",
            "cannot use; set them to NA first."))
  # the all-missing column goes, and one column is left
  ef(mcartst(data.frame(a = c(1, 2), b = NA_real_)),
     "mcartst: too few complete columns or rows remain for Little's test.")
  dd <- stats::setNames(data.frame(1:2, 1:2), c("a", "a"))
  ef(duprows(dd, "a"),
     paste0("duprows: `data` has two columns named a: give the columns by ",
            "position, or rename one."))
})

test_that("empty and one-row tables keep their shape", {
  # one row: is.na() over the columns is a vector, not a matrix
  one <- data.frame(a = NA, b = 1)
  expect_identical(dim(drpempty(one, quiet = TRUE)), c(1L, 1L))
  expect_identical(dim(drpempty(data.frame(a = NA_real_, b = NA_real_),
                                which = "rows", quiet = TRUE)), c(0L, 2L))
  p <- misspat(one)
  expect_identical(cols(p), list("X.", 1L, 100, 1L, "a"))
  expect_identical(
    unname(misssmry(one)),
    c(1, 2, 1, 50, 0, 0, 1, 1)
  )
  expect_identical(capture.output(r <- missgrid(one)), c(r, ""))
  expect_identical(r[3L], if (.uniok()) {
    "       1 \u2588 "
  } else {
    "       1 # "
  })
  x <- drpconst(data.frame())
  expect_identical(attr(x, "dropped"), character(0))
  expect_identical(dim(x), c(0L, 0L))
  expect_identical(clncolnms(character(0)), character(0))
  expect_identical(clncolnms(data.frame()), data.frame())
  z <- missgrid(data.frame(a = 1)[0, , drop = FALSE])
  expect_identical(z[2L], "  (no rows)")
  expect_length(z, 3L)
  d0 <- duprows(data.frame(a = 1:2, b = 3:4)[0, ], quiet = TRUE)
  expect_named(d0, c("a", "b", "dupe_count"))
  expect_identical(nrow(d0), 0L)
  expect_identical(
    as.matrix(mahaoutl(as.matrix(data.frame(a = c(1, 2, 3, 4, 10),
                                            b = c(2, 1, 3, 5, 4))))[, 1:2]),
    as.matrix(mahaoutl(data.frame(a = c(1, 2, 3, 4, 10),
                                  b = c(2, 1, 3, 5, 4)))[, 1:2])
  )
})

test_that("the notes and printed forms say what was done", {
  expect_message(
    out <- drpconst(data.frame(a = c(1, 1), b = 1:2)),
    paste0("-- Constant columns dropped\n   columns 2 -> 1\n",
           "   columns dropped (one value): a"),
    fixed = TRUE
  )
  expect_named(out, "b")
  expect_message(drpconst(data.frame(b = 1:2)),
                 "-- Constant columns dropped\n   columns 1 -> 1\n",
                 fixed = TRUE)
  expect_message(
    expect_identical(clncolnms(c("Total Pop", "x")), c("total_pop", "x")),
    "-- Column names cleaned\n   Total Pop -> total_pop", fixed = TRUE
  )
  expect_message(clncolnms(c("a", "b")),
                 "-- Column names cleaned\n   no name changed", fixed = TRUE)
  f1 <- freqtbl(data.frame(g = c("a", "b", "b"), stringsAsFactors = FALSE), "g")
  # a subset of the columns prints as a plain table
  expect_identical(capture.output(print(f1[, c("pct", "pctvalid")])),
                   c("      pct pctvalid", " 66.66667 66.66667",
                     " 33.33333 33.33333"))
  d <- data.frame(a = c(1, 2, 3, 5), b = c(2, 1, 4, 3))
  expect_identical(
    capture.output(print(cortbl(d))),
    c("Correlations (spearman): 1 pairs", " x y correlation nprs",
      " a b         0.6    4")
  )
  expect_identical(capture.output(print(topcors(d, minabs = 0.9))),
                   c("Top correlations (spearman)",
                     "  (none above the threshold)"))
  mr <- data.frame(a = c(1, NA, NA, 4), b = c(NA, NA, NA, 1))
  expect_identical(
    capture.output(print(missruns(mr))),
    c("Runs of consecutive missing values: 2 ",
      " column start end length", "      b     1   3      3",
      "      a     2   3      2")
  )
  expect_identical(capture.output(print(missruns(data.frame(a = 1:3)))),
                   "Runs of consecutive missing values: 0 ")
  expect_identical(
    capture.output(print(misspat(data.frame(a = NA, b = 1)))),
    c("Missingness patterns, columns in pattern order: a, b ",
      " pattern n_rows pctrows n_missing columns",
      "      X.      1     100         1       a")
  )
})

test_that("the ASCII sketches are used when Unicode is turned off", {
  old <- Sys.getenv("CRIMLAB_ASCII_ONLY", NA)
  Sys.setenv(CRIMLAB_ASCII_ONLY = "1")
  on.exit(if (is.na(old)) {
    Sys.unsetenv("CRIMLAB_ASCII_ONLY")
  } else {
    Sys.setenv(CRIMLAB_ASCII_ONLY = old)
  }, add = TRUE)
  # counts 1, 2, 3 of a top of 3, on five levels: ceiling(5 c / 3)
  expect_identical(inlhist(c(1, 2, 2, 3, 3, 3), bins = 3), ":=#")
  # two bins, as hist(breaks = c(1, 2, 3)): counts 2 and 1
  expect_identical(inlhist(c(1, 2, 3), bins = 2), "#-")
  expect_identical(inlhist(rep(5, 3)), "     #    ")
  g <- data.frame(a = c(1, NA, 3, NA), b = c(NA, NA, 3, 4))
  out <- capture.output(r <- missgrid(g, height = 4))
  expect_identical(
    r,
    c(paste0("-- Missingness map ", strrep("-", 47)), "         ab",
      "       1  #", "       2 ##", "       3   ", "       4 # ", "",
      "  legend: ' ' none  '#' all missing", strrep("-", 66))
  )
})

test_that("BUG: missgrid() keeps a shorter name in its own column", {
  g <- data.frame(a = c(1, NA), bb = c(NA, 2))
  out <- capture.output(r <- missgrid(g))
  # the second letters: none for a, "b" for bb, under the second column
  expect_identical(r[3L], "          b")
})

test_that("Little's test is read out by its p-value", {
  set.seed(5)
  n <- 80
  xx <- round(rnorm(n), 3)
  yy <- round(xx + rnorm(n), 3)
  mc2 <- data.frame(x = xx, y = yy)
  mc2$y[sample(n, 20)] <- NA
  expect_output(print(mcartst(mc2)),
                "^Little's MCAR test: no departure from MCAR detected\n")
  out <- capture.output(print(mcartst(data.frame(a = xx, b = yy))))
  expect_identical(out[1L], "Little's MCAR test: no test possible")
  expect_identical(
    out[4L],
    paste0("  Note: df is 0: with a single missingness pattern there is ",
           "nothing to test. A complete table trivially satisfies MCAR.")
  )
})

test_that("Little's test without complete rows matches the reference", {
  # every row lacks one value: the rank check uses pairwise moments
  set.seed(11)
  n <- 30
  a <- rnorm(n)
  m3 <- data.frame(a = a, b = a + rnorm(n), c = rnorm(n))
  for (i in 1:n) m3[i, (i %% 3) + 1] <- NA
  r <- mcartst(m3)
  # rmoriebricklayer 1.1.1 mcartst() on these data
  expect_equal(unname(unlist(r[1:8])),
               c(2.36043941784393, 3, -0.691066760495222, 0.501041292820027,
                 3, 30, 3, 139),
               tolerance = 1e-10)
  # pairwise moments that are no covariance (an eigenvalue below zero):
  # the EM start is floored to a positive definite matrix
  set.seed(12)
  u <- rnorm(10)
  v <- rnorm(10)
  w <- rnorm(10)
  q <- data.frame(x = c(u, rep(NA, 10), w, rnorm(6)),
                  y = c(u, v, rep(NA, 10), rnorm(6)),
                  z = c(rep(NA, 10), v, -w, rnorm(6)))
  expect_lt(min(eigen(stats::cov(q, use = "pairwise"))$values), 0)
  r <- mcartst(q)
  # the start is floored on the scale of the correlations, the reference's
  # on that of the covariances: one optimum, reached to the tolerance of
  # the EM steps, in a step more or less (the reference took 96)
  expect_equal(unname(unlist(r[1:7])),
               c(3.64923455237147, 6, -0.322935774282878, 0.724020353601097,
                 4, 36, 3),
               tolerance = 1e-7)
  expect_lte(abs(r[[8]] - 96), 3)
})

test_that("near-collinear distances use the floored covariance", {
  set.seed(7)
  x <- round(rnorm(20), 3)
  d <- data.frame(a = x, b = x + 5e-4 * ((1:20) %% 3 - 1))
  # classical: centred by the mean, scaled by the range
  z <- scale(as.matrix(d), colMeans(d), apply(d, 2, function(v) diff(range(v))))
  e <- eigen(stats::cov(z), symmetric = TRUE)
  # the small eigenvalue, 5.2e-9, is above the rank check (1e-8 of the
  # largest) and below the floor (1e-8 times max(1, largest))
  expect_gt(e$values[2], 1e-8 * e$values[1])
  expect_lt(e$values[2], 1e-8)
  sf <- e$vectors %*% diag(pmax(e$values, 1e-8)) %*% t(e$vectors)
  m <- mahaoutl(d, robust = FALSE)
  m <- m[order(m$row), ]
  # a condition number of 2e7 leaves about 1e-9 of rounding
  expect_equal(m$distance, sqrt(rowSums((z %*% solve(sf)) * z)),
               tolerance = 1e-6, ignore_attr = TRUE)
  # robust: a column with MAD 0 is scaled by its sd; a constant one is
  # collinear with nothing to scale by
  r <- data.frame(a = c(1, 1, 1, 1, 2, 5, 3), b = c(2, 4, 1, 3, 6, 5, 7))
  # rmoriebricklayer 1.1.1 mahaoutl() on these data, rows 6, 3, 7
  expect_equal(mahaoutl(r)$distance[1:3],
               c(2.9757081022300254, 1.7459630857421335, 1.5080429530975228),
               tolerance = 1e-12)
  r$a <- 1
  expect_error(mahaoutl(r), "mahaoutl: the columns are collinear",
               fixed = TRUE)
})

test_that("profiles hold for constant and overflowing columns", {
  # 1.5e308 - (-1.5e308) overflows: the median is the weighted mean, 0
  p <- profcols(data.frame(x = c(-1.5e308, 1.5e308)))
  expect_identical(p$median, 0)
  expect_identical(p$mean, 0)
  k <- profcols(data.frame(k = c(3, 3, 3)))
  expect_identical(k$sd, 0)
  expect_identical(k$mad, 0)
})

test_that("dstntst() refuses and warns in its own words", {
  g <- data.frame(group = c("A", "B", "C"), est = c(1, 2, 3),
                  lo = c(0.5, 1.5, 2.5), hi = c(1.5, 2.5, 3.5),
                  q = c("ok", "Z", "E"), stringsAsFactors = FALSE)
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(dstntst(g, " ", "est", "lo", "hi"),
     "dstntst: `grp` must be a single non-empty string.")
  ef(dstntst(g, "group", "est", "lo", "hi", conflvl = 1),
     "dstntst: `conflvl` must be a single number strictly between 0 and 1.")
  ef(dstntst(g, "group", "est", "lo", "hi", pooled = 1),
     "dstntst: `pooled` must be the labels of the rows that pool other groups.")
  expect_warning(
    r <- dstntst(g, "group", "est", "lo", "hi", qual = "q"),
    "dstntst: quality flag(s) \"Z\" not recognised: treated as ok.",
    fixed = TRUE
  )
  expect_identical(r$groups$group, c("A", "B", "C"))
  g$est[2] <- 9
  expect_warning(
    r <- dstntst(g, "group", "est", "lo", "hi"),
    "dstntst: the interval does not contain the estimate for: B; not compared.",
    fixed = TRUE
  )
  expect_identical(r$groups$group, c("A", "C"))
  # B is still outside its interval: warned of before the refusal
  expect_warning(
    ef(dstntst(rbind(g, g[1, ]), "group", "est", "lo", "hi"),
       "dstntst: a group appears twice."),
    "dstntst: the interval does not contain the estimate for: B; not compared.",
    fixed = TRUE
  )
})

test_that("freqtbl() refuses a sort that is not TRUE or FALSE", {
  expect_error(freqtbl(c("a", "b"), sort = NA),
               "freqtbl: `sort` must be TRUE or FALSE.", fixed = TRUE)
})

test_that("mcartst() takes a numeric matrix as its data frame", {
  set.seed(5)
  m <- cbind(x = rnorm(60), y = rnorm(60))
  m[sample(60, 12), "y"] <- NA
  a <- mcartst(m)
  b <- mcartst(as.data.frame(m))
  expect_equal(a$statistic, b$statistic, tolerance = 1e-12)
  expect_identical(a$df, b$df)
})

test_that("inlhist() refuses more bins than it can draw", {
  expect_error(inlhist(1:5, bins = 1e9),
               "inlhist: `bins` must be at most 67108864.", fixed = TRUE)
})

test_that("a number given as invalid text is refused, not read", {
  bad <- rawToChar(as.raw(c(0xff, 0xfe, 0x41)))
  expect_error(missgrid(data.frame(a = c(1, NA)), height = bad),
               "^missgrid: `height` holds text that is not valid")
})

test_that("mcartst() takes columns in different units", {
  # money next to a rate: the rank check is made on the correlations, so
  # the data are not refused, and Little's statistic does not depend on
  # the units of a column
  set.seed(1)
  n <- 300
  d <- data.frame(income = rnorm(n, 5e4, 2e4), rate = rnorm(n, 0.01, 0.003),
                  age = rnorm(n, 40, 10))
  d$rate[d$age > 45] <- NA
  s <- data.frame(income = d$income / 1e4, rate = d$rate * 100, age = d$age)
  r <- mcartst(d)
  expect_equal(r$statistic, mcartst(s)$statistic, tolerance = 1e-5)
  # Little's d2 by hand from the returned estimates
  pt <- apply(!is.na(d), 1L, paste, collapse = "")
  d2 <- 0
  for (k in unique(pt)) {
    rk <- d[pt == k, , drop = FALSE]
    o <- which(!is.na(rk[1L, ]))
    df <- colMeans(rk[, o, drop = FALSE]) - r$mu[o]
    d2 <- d2 + nrow(rk) * drop(df %*% solve(r$sigma[o, o], df))
  }
  expect_equal(r$statistic, d2, tolerance = 1e-10)
  # a constant column is still refused
  d$k <- 0.1
  expect_error(mcartst(d), "^mcartst: the columns are collinear")
})

test_that("clncolnms() numbers repeats past the names already taken", {
  # x_2 is a name given: the repeat of x takes the next free number
  expect_identical(clncolnms(c("x", "x", "x_2"), quiet = TRUE),
                   c("x", "x_3", "x_2"))
  nm <- c("a", "a", "a_2", "a", "a_3", "a_2")
  out <- clncolnms(nm, quiet = TRUE)
  expect_identical(anyDuplicated(out), 0L)
  expect_identical(out, c("a", "a_4", "a_2", "a_5", "a_3", "a_2_2"))
  df <- data.frame(1, 2, 3)
  names(df) <- c("a", "a", "a_2")
  expect_named(clncolnms(df, quiet = TRUE), c("a", "a_3", "a_2"))
  # a tab or a line break is a word break, as a space is
  expect_identical(clncolnms(c("a\tb", "c\nd", "e\001f"), quiet = TRUE),
                   c("a_b", "c_d", "e_f"))
})

test_that("inlhist() bins values whose range overflows", {
  # -1e308 and 1e308 are more than the largest double apart: each end is
  # in its own bin, as the counts by hand give (one in the first of four
  # bins of width 5e307, one in the last)
  old <- Sys.getenv("CRIMLAB_ASCII_ONLY", NA)
  Sys.setenv(CRIMLAB_ASCII_ONLY = "1")
  on.exit(if (is.na(old)) Sys.unsetenv("CRIMLAB_ASCII_ONLY") else
    Sys.setenv(CRIMLAB_ASCII_ONLY = old))
  x <- c(-1e308, 1e308)
  br <- -1e308 + 0:4 * 5e307
  cnt <- tabulate(findInterval(x, br, rightmost.closed = TRUE), 4)
  expect_identical(cnt, c(1L, 0L, 0L, 1L))
  h <- strsplit(inlhist(x, 4), "", fixed = TRUE)[[1L]]
  expect_identical(h != " ", cnt > 0)
  h2 <- strsplit(inlhist(c(-1e308, 0, 1e308), 2), "",
                 fixed = TRUE)[[1L]]
  expect_identical(h2 != " ", c(TRUE, TRUE))
  expect_identical(inlhist(c(-1e308, 1e308), 4),
                   inlhist(c(-1, 1), 4))
})
