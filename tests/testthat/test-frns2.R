test_that("coancestry can make a common heterozygote rarer", {
  expect_lt(gtmp(0.5, 0.5, theta = 0.1)[["locus"]][["mp"]], 0.5)
  expect_gte(gtmp(0.5, hom = TRUE, theta = 0.1)[["locus"]][["mp"]], 0.25)
})

test_that("kinidx gives the textbook indices", {
  p <- 0.1
  q <- 0.2
  k <- function(g1, g2, f, rel) kinidx(g1, g2, f, rel)[["locus"]][["idx"]]
  expect_equal(k(c(1, 2), c(1, 2), c(p, q), "parent"), (p + q) / (4 * p * q))
  expect_equal(k(c(1, 1), c(1, 1), c(p, p), "parent"), 1 / p)
  expect_equal(k(c(1, 1), c(1, 2), c(p, p), "parent"), 1 / (2 * p))
  expect_equal(k(c(1, 2), c(1, 3), c(p, q), "parent"), 1 / (4 * p))
  expect_equal(k(c(1, 2), c(3, 4), c(p, q), "parent"), 0)
  expect_equal(k(c(1, 2), c(3, 4), c(p, q), "sibling"), 0.25)
  expect_equal(k(c(1, 2), c(2, 1), c(p, q), "sibling"),
               0.25 + 0.5 * (p + q) / (4 * p * q) + 0.25 / (2 * p * q))
  expect_equal(k(c(1, 1), c(1, 1), c(p, p), c(0, 0, 1)), 1 / p^2)
  expect_equal(k(c("a", "b"), c("b", "c"), c(p, q), "half"),
               0.5 + 0.5 / (4 * q))
  expect_equal(k(c(1, 2), c(1, 2), c(p, q), c(1, 0, 0)), 1)
  g1 <- matrix(c(1, 1, 2, 1), 2, dimnames = list(c("D3", "vWA"), NULL))
  r <- kinidx(g1, g1, matrix(c(p, p, q, p), 2))
  expect_equal(r[["locus"]][["locus"]], c("D3", "vWA"))
  expect_equal(r[["log10idx"]], log10((p + q) / (4 * p * q)) + log10(1 / p))
  expect_equal(r[["post"]], 1 / (1 + 10^-r[["log10idx"]]))
})

test_that("patidx gives the textbook indices", {
  p <- 0.1
  q <- 0.2
  pi <- function(...) patidx(...)[["locus"]][["idx"]]
  expect_equal(pi(c(1, 2), c(2, 2), c(p, q), mother = c(1, 1)), 1 / q)
  expect_equal(pi(c(1, 2), c(2, 3), c(p, q), mother = c(1, 1)), 1 / (2 * q))
  expect_equal(pi(c(1, 2), c(1, 2), c(p, q), mother = c(1, 2)), 1 / (p + q))
  expect_equal(pi(c(1, 1), c(1, 3), c(p, p), mother = c(1, 2)), 1 / (2 * p))
  expect_equal(pi(c(1, 2), c(3, 3), c(p, q), mother = c(1, 1)), 0)
  expect_true(is.nan(pi(c(1, 2), c(1, 2), c(p, q), mother = c(3, 3))))
  # no mother: the parent-child kinship index
  expect_equal(pi(c(1, 2), c(1, 2), c(p, q)), (p + q) / (4 * p * q))
  # brute force: every trio drawn from three alleles
  f <- c(0.2, 0.3, 0.5)
  gs <- rbind(c(1, 1), c(1, 2), c(1, 3), c(2, 2), c(2, 3), c(3, 3))
  trans <- function(m, a, ch) {
    mean(outer(m, a, function(x, y) {
      pmin(x, y) == min(ch) & pmax(x, y) == max(ch)
    }))
  }
  for (ic in 1:6) for (im in 1:6) for (ia in 1:6) {
    ch <- gs[ic, ]
    den <- sum(vapply(1:6, function(j) {
      g <- gs[j, ]
      pg <- if (g[1] == g[2]) f[g[1]]^2 else 2 * f[g[1]] * f[g[2]]
      pg * trans(gs[im, ], g, ch)
    }, numeric(1)))
    num <- trans(gs[im, ], gs[ia, ], ch)
    got <- pi(ch, gs[ia, ], f[ch], mother = gs[im, ])
    if (den == 0) expect_true(is.nan(got)) else expect_equal(got, num / den)
  }
})

test_that("mixlr follows inclusion-exclusion", {
  p <- c(a = 0.1, b = 0.2, c = 0.05, d = 0.15)
  r <- mixlr(p, character(0), 0, c("c", "d"), 1)
  expect_equal(unname(r), c(1, 2 * 0.05 * 0.15, 1 / (2 * 0.05 * 0.15)))
  # nothing to explain, one unknown: random man not excluded
  expect_equal(mixlr(p, character(0), 1, character(0), 1)[["num"]], 0.5^2)
  # two unknowns carry all four alleles: 4! ordered assignments
  expect_equal(mixlr(p, names(p), 2, names(p), 2)[["num"]], 24 * prod(p))
  # three alleles cannot sit in one person
  expect_equal(mixlr(p, c("a", "b", "c"), 1, "a", 1)[["num"]], 0)
  # brute force, two unknowns must show c
  s <- sum(p)
  expect_equal(mixlr(p, "c", 2, "c", 2)[["num"]], s^4 - (s - 0.05)^4)
  expect_equal(mixlr(p, "a", 0, "a", 1)[["num"]], 0)
})

test_that("dbsrch bounds the exact chance by np", {
  r <- dbsrch(1e-3, 100)
  expect_equal(r[["np"]], 0.1)
  expect_equal(r[["any"]], 1 - (1 - 1e-3)^100)
  expect_lte(r[["any"]], r[["np"]])
  expect_equal(dbsrch(0.5, 10)[["np"]], 1)
  expect_equal(dbsrch(1e-6, 1e4, pop = 1e6)[["post"]],
               1 / (1 + (1e6 - 1e4) * 1e-6))
  expect_equal(dbsrch(1e-6, 1e6, pop = 1e6)[["post"]], 1)
  expect_equal(dbsrch(1e-20, 10)[["any"]], 1e-19)
})

test_that("cllrmin recalibrates without changing the order", {
  s <- c(50, 200, 8)
  d <- c(0.01, 0.2, 1.5)
  r <- cllrmin(s, d)
  expect_equal(r[["cllr"]], cllr(s, d))
  expect_equal(r[["min"]], 0)  # the two sets do not overlap
  expect_equal(r[["cal"]], r[["cllr"]])
  # scaling every ratio changes the cost, not the minimum
  s <- c(3, 0.7, 20, 5)
  d <- c(0.2, 1, 4, 0.05)
  a <- cllrmin(s, d)
  b <- cllrmin(s * 1000, d * 1000)
  expect_equal(a[["min"]], b[["min"]])
  expect_gt(b[["cllr"]], a[["cllr"]])
  expect_true(a[["min"]] <= a[["cllr"]] && a[["min"]] > 0)
  # by hand: order 0.05d 0.2d 0.7s 1d 3s 4d 5s 20s
  # blocks: d d | s d | s d | s s  with same-source share 0, 0.5, 0.5, 1
  expect_equal(a[["min"]], 0.5 * (2 * log2(2) / 4 + 2 * log2(2) / 4))
  cal <- lrcal(s, d)
  expect_equal(cal[["same"]], c(1, 1, Inf, Inf))
  expect_equal(cal[["diff"]], c(0, 1, 1, 0))
  # calibrated ratios are their own minimum
  fin <- function(x) pmin(pmax(x, 1e-300), 1e300)
  expect_equal(cllr(fin(cal[["same"]]), fin(cal[["diff"]])), a[["min"]])
  # ties across the two sets are pooled
  expect_equal(cllrmin(c(2, 2), c(2, 2))[["min"]], 1)
  # an uninformative reversed system
  expect_equal(cllrmin(c(0.1, 0.2), c(5, 9))[["min"]], 1)
})

test_that("the plots return what they draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  s <- c(50, 200, 8, 0.5)
  d <- c(0.01, 0.2, 1.5, 0.001)
  t <- tipplt(s, d)
  expect_equal(t[["same"]][["share"]], c(1, 0.75, 0.5, 0.25))
  expect_equal(t[["same"]][["log10lr"]], sort(log10(s)))
  e <- eceplt(s, d, lo = c(-1, 0, 1))
  m <- cllrmin(s, d)
  expect_equal(e[["ece"]][2], m[["cllr"]])
  expect_equal(e[["min"]][2], m[["min"]])
  expect_equal(e[["ref"]][2], 1)
  expect_true(all(e[["min"]] <= e[["ece"]] + 1e-12))
  expect_true(all(e[["min"]] <= e[["ref"]] + 1e-12))
})

test_that("bad input to the second forensic set stops in R", {
  expect_error(kinidx(1:3, 1:2, c(0.1, 0.2)), "`g1`")
  expect_error(kinidx(c(1, 2), matrix(1, 2, 2), c(0.1, 0.2)), "1 rows")
  expect_error(kinidx(c(1, 2), c(1, NA), c(0.1, 0.2)), "`g2`")
  expect_error(kinidx(c(1, 2), c(1, 2), c(0.1, 1.2)), "`freq`")
  expect_error(kinidx(c(1, 1), c(1, 2), c(0.1, 0.2)), "homozygote")
  expect_error(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), "twin"), "`rel`")
  expect_error(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), c(0.5, 0.2, 0.2)), "`rel`")
  expect_error(patidx(c(1, 2), c(1, 2), c(0.1, 0.2), mother = 1), "`mother`")
  expect_error(mixlr(c(0.1, 0.2), "a", 1, "a", 1), "named")
  expect_error(mixlr(c(a = 0.7, b = 0.6), "a", 1, "a", 1), "at most 1")
  expect_error(mixlr(c(a = 0.1), "z", 1, "a", 1), "`np`")
  expect_error(mixlr(c(a = 0.1), "a", 1, "a", 1.5), "`xd`")
  expect_error(dbsrch(2, 10), "at most 1")
  expect_error(dbsrch(0.1, 0), "`n`")
  expect_error(dbsrch(0.1, 10, pop = 5), "`pop`")
  expect_error(cllrmin(c(1, -1), 1), "`same`")
  expect_error(lrcal(1, numeric(0)), "`diff`")
  expect_error(tipplt(0, 1), "`same`")
  expect_error(eceplt(1, 1, lo = NA), "`lo`")
})
