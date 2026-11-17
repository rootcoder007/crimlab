test_that("cemtch() keeps strata with both groups, weighted to the treated", {
  set.seed(31)
  n <- 200
  d <- data.frame(pov = rbinom(n, 1, 0.4), age = runif(n, 14, 18),
                  race = sample(c("a", "b", "c"), n, TRUE))
  d$lead <- rbinom(n, 1, plogis(-1 + 1.5 * d$pov))
  m <- cemtch(d, "lead", c("pov", "age", "race"), cut = list(age = 3))
  # written out: age in three equal bins, the strata by paste()
  br <- seq(min(d$age), max(d$age), length.out = 4)
  ab <- findInterval(d$age, br, rightmost.closed = TRUE, all.inside = TRUE)
  s <- paste(d$pov, ab, d$race)
  nt <- as.vector(tapply(d$lead, s, sum)[s])
  nc <- as.vector(tapply(1 - d$lead, s, sum)[s])
  keep <- unname(nt > 0 & nc > 0)
  mt <- sum(d$lead[keep])
  mc <- sum(1 - d$lead[keep])
  w <- nt / nc * mc / mt
  w[d$lead == 1] <- 1
  w[!keep] <- 0
  expect_identical(m$keep, keep)
  expect_equal(m$w, unname(w), tolerance = 1e-14)
  expect_identical(m$cut$age, br)
  expect_identical(as.vector(m$tab["matched", ]), as.integer(c(mc, mt)))
  expect_identical(m$nstrata[["all"]], length(unique(s)))
  expect_identical(m$nstrata[["kept"]], length(unique(s[keep])))
  # strata numbered from 1 in order of first appearance among kept rows
  expect_identical(m$strata[keep], match(s[keep], unique(s[keep])))
  # the weighted controls match the treated in every kept stratum
  sw <- tapply(m$w[keep] * (1 - d$lead[keep]), s[keep], sum) / mc
  st <- tapply(d$lead[keep], s[keep], sum) / mt
  expect_equal(sw, st, tolerance = 1e-14)
})

test_that("cemtch()'s L1 is half the summed share differences", {
  set.seed(32)
  n <- 150
  d <- data.frame(x = rnorm(n), g = sample(1:2, n, TRUE))
  d$t <- rbinom(n, 1, plogis(d$x))
  m <- cemtch(d, "t", c("x", "g"), cut = list(x = 2),
              lcut = list(x = c(-10, -1, 0, 1, 10)))
  cell <- paste(findInterval(d$x, c(-10, -1, 0, 1, 10)), d$g)
  l1 <- function(w) {
    f <- tapply(w * d$t, cell, sum) / sum(w * d$t)
    g <- tapply(w * (1 - d$t), cell, sum) / sum(w * (1 - d$t))
    sum(abs(f - g)) / 2
  }
  expect_equal(m$l1[["before"]], l1(rep(1, n)), tolerance = 1e-14)
  expect_equal(m$l1[["after"]], l1(m$w), tolerance = 1e-14)
  expect_error(cemtch(d, "t", "x", cut = list(x = c(0, 1))),
               paste("cemtch: `cut$x` must be increasing breaks covering the",
                     "column."), fixed = TRUE)
  expect_error(cemtch(d, "t", "x", cut = list(z = 2)),
               "cemtch: `cut` must be a list named by columns of `x`.",
               fixed = TRUE)
  expect_error(cemtch(d, "x", "g"), "cemtch: column `x` must hold 0 or 1.",
               fixed = TRUE)
})

test_that("cemtch() prints its strata and refuses bad columns and bins", {
  set.seed(33)
  n <- 120
  d <- data.frame(x = rnorm(n), k = 2, g = sample(c("a", "b"), n, TRUE))
  d$t <- rbinom(n, 1, 0.5)
  m <- cemtch(d, "t", c("x", "k", "g"), cut = list(x = 3))
  # a constant column is one bin of width 1 about its value
  expect_identical(m$cut$k, c(1.5, 2.5))
  br <- seq(min(d$x), max(d$x), length.out = 4)
  s <- paste(findInterval(d$x, br, rightmost.closed = TRUE,
                          all.inside = TRUE), d$g)
  both <- tapply(d$t, s, function(v) any(v == 1) && any(v == 0))
  head <- sprintf("cemtch: coarsened exact matching, %d of %d strata kept",
                  sum(both), length(both))
  out <- utils::capture.output(v <- print(m))
  expect_identical(out[1L], head)
  expect_match(out[length(out)], "^L1 imbalance: ")
  expect_identical(v, m)
  q <- d
  q$x[3] <- Inf
  expect_error(cemtch(q, "t", "x"), "cemtch: column `x` must be finite.",
               fixed = TRUE)
  expect_error(cemtch(d, "t", "x", cut = list(x = 0)),
               "cemtch: `cut$x` must be a whole number of bins from 1.",
               fixed = TRUE)
  expect_error(cemtch(d, "t", "x", cut = list(x = 2.5)),
               "cemtch: `cut$x` must be a whole number of bins from 1.",
               fixed = TRUE)
})

test_that("cemtch() refuses a table of one row in words", {
  expect_error(cemtch(data.frame(t = 1, x = 2), "t", "x"),
               "^cemtch: `data` needs two rows at least")
})
