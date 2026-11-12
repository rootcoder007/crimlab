# Values of the forensic kernels are checked in C++ (src/t_frns.cpp); here
# what R forms itself from them, the plots, the labels and the messages.

test_that("kinidx sums the loci in R and labels them", {
  p <- 0.1
  q <- 0.2
  g1 <- matrix(c(1, 1, 2, 1), 2, dimnames = list(c("D3", "vWA"), NULL))
  r <- kinidx(g1, g1, matrix(c(p, p, q, p), 2))
  expect_identical(r[["locus"]][["locus"]], c("D3", "vWA"))
  expect_equal(r[["log10idx"]], log10((p + q) / (4 * p * q)) + log10(1 / p),
               tolerance = 1e-10)
  expect_equal(r[["post"]], 1 / (1 + 10^-r[["log10idx"]]), tolerance = 1e-10)
})

test_that("mixlr forms the ratio of its two sides", {
  p <- c(a = 0.1, b = 0.2, c = 0.05, d = 0.15)
  r <- mixlr(p, character(0), 0, c("c", "d"), 1)
  expect_equal(unname(r), c(1, 2 * 0.05 * 0.15, 1 / (2 * 0.05 * 0.15)),
               tolerance = 1e-10)
})

test_that("dbsrch bounds the exact chance by np", {
  r <- dbsrch(1e-3, 100)
  expect_equal(r[["np"]], 0.1)
  expect_equal(r[["any"]], 1 - (1 - 1e-3)^100, tolerance = 1e-10)
  expect_lte(r[["any"]], r[["np"]])
  expect_equal(dbsrch(0.5, 10)[["np"]], 1, tolerance = 1e-10)
  expect_equal(dbsrch(1e-6, 1e4, pop = 1e6)[["post"]],
               1 / (1 + (1e6 - 1e4) * 1e-6), tolerance = 1e-10)
  expect_equal(dbsrch(1e-6, 1e6, pop = 1e6)[["post"]], 1, tolerance = 1e-10)
  expect_equal(dbsrch(1e-20, 10)[["any"]], 1e-19, tolerance = 1e-10)
})

test_that("cllrmin and lrcal give the calibrated ratios", {
  s <- c(50, 200, 8)
  d <- c(0.01, 0.2, 1.5)
  r <- cllrmin(s, d)
  expect_equal(r[["cal"]], r[["cllr"]], tolerance = 1e-10)
  s <- c(3, 0.7, 20, 5)
  d <- c(0.2, 1, 4, 0.05)
  a <- cllrmin(s, d)
  cal <- lrcal(s, d)
  expect_equal(cal[["same"]], c(1, 1, Inf, Inf), tolerance = 1e-10)
  expect_equal(cal[["diff"]], c(0, 1, 1, 0), tolerance = 1e-10)
  # calibrated ratios are their own minimum
  fin <- function(x) pmin(pmax(x, 1e-300), 1e300)
  expect_equal(cllr(fin(cal[["same"]]), fin(cal[["diff"]])), a[["min"]],
               tolerance = 1e-10)
})

test_that("the plots return what they draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  s <- c(50, 200, 8, 0.5)
  d <- c(0.01, 0.2, 1.5, 0.001)
  t <- tipplt(s, d)
  expect_equal(t[["same"]][["share"]], c(1, 0.75, 0.5, 0.25))
  expect_equal(t[["same"]][["log10lr"]], sort(log10(s)), tolerance = 1e-10)
  e <- eceplt(s, d, lo = c(-1, 0, 1))
  m <- cllrmin(s, d)
  expect_equal(e[["ece"]][2], m[["cllr"]], tolerance = 1e-10)
  expect_equal(e[["min"]][2], m[["min"]], tolerance = 1e-10)
  expect_equal(e[["ref"]][2], 1, tolerance = 1e-10)
  expect_true(all(e[["min"]] <= e[["ece"]] + 1e-12))
  expect_true(all(e[["min"]] <= e[["ref"]] + 1e-12))
})

test_that("bad input to the second forensic set stops in R", {
  errs("kinidx", list(
    list("`g1`", 1:3, 1:2, c(0.1, 0.2)),
    list("1 rows", c(1, 2), matrix(1, 2, 2), c(0.1, 0.2)),
    list("`g2`", c(1, 2), c(1, NA), c(0.1, 0.2)),
    list("`freq`", c(1, 2), c(1, 2), c(0.1, 1.2)),
    list("homozygote", c(1, 1), c(1, 2), c(0.1, 0.2)),
    list("`rel`", c(1, 2), c(1, 2), c(0.1, 0.2), "twin"),
    list("`rel`", c(1, 2), c(1, 2), c(0.1, 0.2), c(0.5, 0.2, 0.2))
  ))
  errs("patidx", list(
    list("`mother`", c(1, 2), c(1, 2), c(0.1, 0.2), mother = 1)
  ))
  errs("mixlr", list(
    list("named", c(0.1, 0.2), "a", 1, "a", 1),
    list("at most 1", c(a = 0.7, b = 0.6), "a", 1, "a", 1),
    list("`np`", c(a = 0.1), "z", 1, "a", 1),
    list("`xd`", c(a = 0.1), "a", 1, "a", 1.5)
  ))
  errs("dbsrch", list(
    list("at most 1", 2, 10),
    list("`n`", 0.1, 0),
    list("`pop`", 0.1, 10, pop = 5)
  ))
  errs("cllrmin", list(list("`same`", c(1, -1), 1)))
  errs("lrcal", list(list("`diff`", 1, numeric(0))))
  errs("tipplt", list(list("`same`", 0, 1)))
  errs("eceplt", list(list("`lo`", 1, 1, lo = NA)))
})
