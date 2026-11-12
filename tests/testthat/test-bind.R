# The compiled binding (src/bind.cpp, src/rchk.cpp) reproduces what the R
# wrappers did, including R's own rules for odd inputs; each test names
# the R behaviour it pins down.

test_that("classed numbers are read as their data, dates are refused", {
  # is.numeric() is TRUE for a class without its own method, FALSE for
  # Date, POSIXt and difftime (base's methods)
  x <- structure(c(3, 4), class = "foo")
  expect_identical(vodtst(x, c(1, 1)), vodtst(c(3, 4), c(1, 1)))
  expect_error(vodtst(as.Date("2020-01-01") + 0:1, c(1, 1)), "`day`")
  expect_error(vodtst(as.POSIXct("2020-01-01", tz = "UTC") + 0:1, c(1, 1)),
               "`day`")
  expect_error(vodtst(as.difftime(c(1, 2), units = "days"), c(1, 1)), "`day`")
  expect_error(vodtst(factor(c(1, 2)), c(1, 1)), "`day`")
  # as.vector() drops the class of a numeric matrix
  d <- structure(matrix(1, 2, 2), class = "difftime", units = "days")
  expect_identical(dyadrt(d, c(1, 1)), dyadrt(matrix(1, 2, 2), c(1, 1)))
  expect_identical(dyadrt(matrix(1L, 2, 2), c(1, 1)),
                   dyadrt(matrix(1, 2, 2), c(1, 1)))
  expect_error(dyadrt(matrix(factor(1:4), 2, 2), c(1, 1)), "`count`")
})

test_that("q is NA in every form all(is.na(q)) accepts", {
  ref <- gtmp(c(0.2, 0.3), hom = TRUE)
  for (q in list(c(NA, NA), c(NA_character_, NA), factor(c(NA, NA)),
                 as.Date(c(NA, NA)), list(NA, NA_real_))) {
    expect_identical(gtmp(c(0.2, 0.3), q, hom = TRUE), ref)
  }
  expect_error(gtmp(c(0.2, 0.3), list(NA, 1), hom = TRUE), "`q`")
  expect_error(gtmp(c(0.2, 0.3), factor(c(NA, "a")), hom = TRUE), "`q`")
  expect_error(gtmp(c(0.2, 0.3), as.Date(c(NA, 1)), hom = TRUE), "`q`")
  # is.na() warns on a type that is not a vector, and all() of a length-0
  # answer is TRUE: the length check then refuses q
  expect_warning(expect_error(gtmp(0.2, new.env(), hom = TRUE), "`q`"),
                 "type 'environment'")
  expect_error(gtmp(c(0.2, 0.3), c(NA, 1i), hom = TRUE), "`q`")
  expect_identical(gtmp(c(0.2, 0.3), c(NA, NA_complex_), hom = TRUE), ref)
  expect_error(gtmp(c(0.2, 0.3), as.raw(1:2), hom = TRUE), "`q`")
  expect_error(gtmp(0.2, 1L), "`q`")
  expect_identical(gtmp(0.2, 0L, hom = TRUE), gtmp(0.2, hom = TRUE))
})

test_that("eceplt takes the limits of calibrated ratios at extreme odds", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  s <- c(2, 3)
  d <- 0.5
  # the ratios separate perfectly: calibrated, same-source ratios are Inf
  # and different-source ratios 0, and their cost is 0 at any odds (it was
  # NaN, or an error, at odds of 10^-400 and 10^400)
  e <- eceplt(s, d, lo = c(-400, 0, 400))
  expect_identical(e[["min"]], c(0, 0, 0))
  # at even odds the cost is cllr; far out it tends to 0
  cl <- 0.5 * (mean(log2(1 + 1 / s)) + mean(log2(1 + d)))
  expect_equal(e[["ece"]], c(0, cl, 0), tolerance = 1e-12)
  expect_equal(e[["ref"]], c(0, 1, 0), tolerance = 1e-12)
  # s * odds underflows: log(1 + 1 / x) is taken as -log(s) - lo * log(10)
  want <- (1e-100 * 400 * log(10) + 1e-100) / log(2)
  expect_equal(eceplt(1e-300, 1, lo = -100)[["ece"]], want, tolerance = 1e-12)
})

test_that("a genotype of two values becomes a one-row table", {
  r <- kinidx(c(1, 2), c(1, 2), c(0.1, 0.2))
  # matrix(x, 1L) gives a factor's labels and drops a date's class
  expect_identical(kinidx(factor(c("1", "2")), c(1, 2), c(0.1, 0.2)), r)
  dt <- as.Date("2020-01-01") + 0:1
  expect_identical(kinidx(dt, dt, c(0.1, 0.2)), r)
  expect_identical(kinidx(c(1L, 2L), c("1", "2"), c(0.1, 0.2)), r)
  expect_error(kinidx(c(TRUE, FALSE), c(1, 2), c(0.1, 0.2)), "`g1`")
  expect_error(kinidx(list(1, 2), c(1, 2), c(0.1, 0.2)), "`g1`")
  # alleles are compared as as.character() writes them: 15 digits
  expect_identical(kinidx(c(0.1 + 0.2, 0.3), c(0.3, 0.3), c(0.1, 0.1)),
                   kinidx(c(0.3, 0.3), c(0.3, 0.3), c(0.1, 0.1)))
})

test_that("many alleles are coded alike by value and by label", {
  # whole numbers are coded by value, other labels as text; both must agree
  # on every locus, however many alleles share the hash table
  g <- cbind(1:3000, 3001:6000)
  h <- cbind(3001:6000, 1:3000)
  f <- matrix(0.001, 3000, 2)
  t <- function(x) matrix(as.character(x), ncol = 2)
  expect_identical(kinidx(g, h, f), kinidx(t(g), t(h), f))
  expect_identical(patidx(g, h, f, g), patidx(t(g), t(h), f, t(g)))
  expect_error(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), NA_character_), "`rel`")
})

test_that("loci are named as names(g[, 1L]) names them", {
  f <- c(0.1, 0.2)
  m <- function(...) matrix(1:2, 1, dimnames = list(...))
  loc <- function(g) kinidx(g, c(1, 2), f)[["locus"]][["locus"]]
  # a one-row table keeps a name only when exactly one of its row and
  # column names exists
  expect_identical(loc(m("r", NULL)), "r")
  expect_identical(loc(m(NULL, c("x", "y"))), "x")
  expect_identical(loc(m("r", c("x", "y"))), "1")
  g <- rbind(D3 = c(1, 2), vWA = c(3, 3))
  r <- kinidx(g, g, rbind(c(0.1, 0.2), c(0.3, 0.3)))
  expect_identical(r[["locus"]][["locus"]], c("D3", "vWA"))
})

test_that("kernel errors keep the call the R helper made", {
  e <- tryCatch(kinidx(matrix(1, 0, 2), matrix(1, 0, 2), matrix(0.1, 0, 2)),
                error = identity)
  expect_identical(conditionMessage(e),
                   "crimlab: a size or an index is out of range")
  expect_identical(
    conditionCall(e),
    quote(.idx(g1, .Call(C_clab_kinidx, .code(g1, u), .code(g2, u), f, k)))
  )
  z <- matrix(1, 0, 2)
  e <- tryCatch(patidx(z, z, matrix(0.1, 0, 2), z), error = identity)
  expect_identical(conditionMessage(e),
                   "crimlab: a size or an index is out of range")
})

test_that("t0 + s past the largest integer is refused, without a warning", {
  # the old R code added the two as integers: NA with a warning, then an
  # error from if (); the sum is now formed in double and checked
  expect_no_warning(
    expect_error(incapr(c(2, 1), .Machine$integer.max, 1),
                 "`lam` must have length at least t0 + s + 1.", fixed = TRUE)
  )
  expect_error(incapr(c(2, 1), 1, .Machine$integer.max),
               "`lam` must have length at least t0 + s + 1.", fixed = TRUE)
})

test_that("z must be two labels, ordered by code point", {
  expect_identical(sntbnd(c(0, 1), c("b", "a"))[["term"]],
                   c("a", "b", "effect"))
  # upper case sorts before lower case on every machine
  expect_identical(sntbnd(c(0, 1), c("a", "B"))[["term"]],
                   c("B", "a", "effect"))
  expect_identical(sntbnd(c(0, 1, 1, 0), factor(c("x", "y", "x", "y"))),
                   sntbnd(c(0, 1, 1, 0), c("x", "y", "x", "y")))
  expect_identical(sntbnd(c(0, 1), c(2, 10))[["term"]],
                   c("10", "2", "effect"))
  expect_error(sntbnd(c(0, 1, 1), c("a", "b", "c")), "exactly two")
  expect_error(sntbnd(c(0, 1), as.Date("2020-01-01") + 0:1), "labels")
  expect_error(sntbnd(c(0, 1), as.difftime(c(1, 2), units = "days")),
               "labels")
  e <- new.env()
  assign("a", 1, e)
  assign("b", 2, e)
  expect_error(sntbnd(c(0, 1), e), "labels")
  expect_error(sntbnd(c(0, 1), list("a", NA)), "`z`")
  expect_identical(sntbnd(c(0, 1), list("a", "b")),
                   sntbnd(c(0, 1), c("a", "b")))
  expect_error(sntbnd(c(0, 1), c("a", "b"), treat = list("a")), "`treat`")
})

test_that("a long message is kept whole, as stop() keeps it", {
  long <- strrep("x", 300)
  expect_error(proof(c(long, "y")), paste0("no proof recorded for: ", long,
                                           ", y."), fixed = TRUE)
})

test_that("names are matched as R matches strings", {
  p <- stats::setNames(c(0.1, 0.2), c("a", NA))
  expect_identical(mixlr(p, "a", 1, "a", 1)[["lr"]], 1)
  b <- "\xe9"
  Encoding(b) <- "bytes"
  u <- "é"
  # a "bytes" string equals only a "bytes" string
  expect_error(vodtst(stats::setNames(c(1, 1), c(b, "x")), c(1, 1), ref = u),
               "names no group")
  expect_identical(
    vodtst(stats::setNames(c(1, 1), c(b, "x")), c(1, 1), ref = b)[["ratio"]],
    c(1, 1)
  )
  l1 <- iconv(u, "UTF-8", "latin1")
  expect_identical(
    vodtst(stats::setNames(c(1, 2), c("x", l1)), c(1, 1), ref = u),
    vodtst(stats::setNames(c(1, 2), c("x", l1)), c(1, 1), ref = 2)
  )
})

test_that("the kernel entry points still give the exported results", {
  v <- vodtst(c(60, 40), c(50, 50), shift = c(1, 2))
  r <- .Call(C_clab_vod, c(60, 40), c(50, 50), 0L, c(1, 2))
  expect_identical(r[["adj"]], v[["adj"]])
  expect_identical(.Call(C_clab_orrisk, 3, 0.4, 0.5),
                   unname(orrisk(3, 0.4, 0.5)[1:2]))
  expect_identical(.Call(C_clab_outtst, c(0.1, 0.9), c(1, 1), 0.3),
                   unname(outtst(c(0.1, 0.9), 0.3)))
  expect_identical(.Call(C_clab_dspbm, c(9, 1), c(9, 3), c(9, 6), 0L)[["rd"]],
                   dspbm(c(9, 1), c(9, 3), c(9, 6))[["rd"]])
  m <- matrix(c(9, 2, 1, 3), 2)
  expect_identical(
    .Call(C_clab_dyadrt, c(9, 2, 1, 3), c(0L, 1L, 0L, 1L), c(0L, 0L, 1L, 1L),
          c(8, 2))[["ratio"]],
    dyadrt(m, c(8, 2))[["ratio"]]
  )
  expect_identical(.Call(C_clab_gtmp, 0.2, 0, TRUE, 0)[["log10mp"]],
                   gtmp(0.2, hom = TRUE)[["log10mp"]])
  expect_identical(.Call(C_clab_cllr, 2, 0.5), cllr(2, 0.5))
  expect_identical(.Call(C_clab_lrcal, c(0.5, 2), c(FALSE, TRUE))[["min"]],
                   cllrmin(2, 0.5)[["min"]])
  expect_identical(.Call(C_clab_mixp, c(0.1, 0.2), c(TRUE, FALSE), 1L),
                   mixlr(c(a = 0.1, b = 0.2), "a", 1, "a", 1)[["num"]])
  f <- c(0.1, 0.2)
  expect_identical(.Call(C_clab_kinidx, 1:2, 1:2, f, c(0, 1, 0)),
                   kinidx(c(1, 2), c(1, 2), f)[["locus"]][["idx"]])
  expect_identical(.Call(C_clab_patidx, 1:2, c(1L, 3L), 1:2, f),
                   patidx(c(1, 2), c(1, 2), f, c(1, 3))[["locus"]][["idx"]])
  cs <- clrsmp(c(5, 3), c(10, 5), c(4, 3), c(6, 7))
  expect_identical(.Call(C_clab_clrsmp, c(5, 3), c(10, 5), c(4, 3),
                         c(6, 7))[["tot"]], unname(cs[["total"]][1:5]))
  expect_identical(.Call(C_clab_expbnd, c(2, 3), c(1, 1), 0L, 1.5)[["hi"]],
                   expbnd(c(2, 3), c(1, 1), 1.5)[["hi"]])
  expect_identical(.Call(C_clab_misrd, 0.3, 0.1, 0.4, 0.9, 0.8),
                   unname(misrd(0.3, 0.1, 0.4, 0.9, 0.8)))
  expect_identical(.Call(C_clab_rtmhs, c(5, 1), c(4, 2), c(1, 1), 4),
                   unname(rtmhs(c(5, 1), c(4, 2), 4)))
  expect_identical(.Call(C_clab_incap, c(1, 4), c(0.3, 0.6), c(2, 2),
                         c(0.7, 0.3))[["agg"]],
                   unlist(incap(c(1, 4), c(0.3, 0.6), 2,
                                share = c(0.7, 0.3))[3, c(2, 5, 6)],
                          use.names = FALSE))
  expect_identical(.Call(C_clab_incapr, c(3, 2, 1), 0L, 1L, 0.2),
                   unname(incapr(c(3, 2, 1), 0, 1, 0.2)))
  expect_identical(.Call(C_clab_sellab, c(1, 0), c(TRUE, TRUE),
                         c(TRUE, FALSE), c(1, 1)),
                   unname(sellab(c(1, 0), c(TRUE, TRUE), c(TRUE, FALSE))))
  expect_identical(.Call(C_clab_losbnd, c(1, 2), 0.5),
                   unname(losbnd(c(1, 2), 0.5)))
  expect_identical(.Call(C_clab_sntbnd, c(0, 1), c(0L, 1L), c(1, 1))[5],
                   sntbnd(c(0, 1), c("a", "b"))[["lower"]][1])
  expect_identical(.Call(C_clab_hzdsel, 0.3, 0.4, 0.1, c(0.9, 0.8),
                         c(0.95, 0.9)),
                   unname(hzdsel(0.3, 0.4, 0.1, c(0.9, 0.8), c(0.95, 0.9))))
})

test_that("lengths past INT_MAX are refused, never truncated", {
  # no test can allocate such a vector: C_clab_tlen runs the one length
  # check of the binding and of the kernel entry points on a given length
  expect_identical(.Call(C_clab_tlen, 5, "x"), 5L)
  expect_identical(.Call(C_clab_tlen, 2147483647, FALSE), 2147483647L)
  expect_error(.Call(C_clab_tlen, 2147483648, "x"),
               "`x` is too long for crimlab: more than 2147483647 values.",
               fixed = TRUE)
  expect_error(.Call(C_clab_tlen, 3e9, FALSE),
               "crimlab: a size or an index is out of range", fixed = TRUE)
  # a length is printed as length() gives it: an integer, a double above
  expect_identical(.Call(C_clab_tlen, 7, NULL), as.character(7L))
  expect_identical(.Call(C_clab_tlen, 3e9, NULL), as.character(3e9))
})

test_that("cllrmin forms the calibration loss from the full cost", {
  # cal = cllr - min is formed in the binding; cllr is checked against its
  # closed form, min against the cost of the calibrated ratios
  s <- c(3, 0.7, 20, 5)
  d <- c(0.2, 1, 4, 0.05)
  a <- cllrmin(s, d)
  full <- 0.5 * (mean(log2(1 + 1 / s)) + mean(log2(1 + d)))
  expect_equal(a[["cllr"]], full, tolerance = 1e-12)
  expect_equal(a[["cal"]], full - a[["min"]], tolerance = 1e-12)
  expect_gt(a[["cal"]], 0)
})
