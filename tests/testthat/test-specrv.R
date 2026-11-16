sc_dat <- function(n = 150, eff = 0.4, seed = 61) {
  set.seed(seed)
  d <- data.frame(t = rbinom(n, 1, 0.5), age = runif(n, 18, 60),
                  a = rnorm(n), b = rpois(n, 2))
  d$y1 <- eff * d$t + 0.3 * d$a + 0.01 * d$age + rnorm(n)
  d$y2 <- eff * d$t + 0.2 * d$b + rnorm(n)
  d
}

test_that("specrv() fits every outcome with every subset of the options", {
  d <- sc_dat()
  r <- specrv(d, c("y1", "y2"), "t", x = "age", opt = c("a", "b"),
              se = "iid", nsim = 20, seed = 1)
  sp <- r$spec
  expect_identical(sp$outcome, rep(c("y1", "y2"), each = 4))
  expect_identical(sp$ctrl, rep(c("none", "a", "b", "a + b"), 2))
  ctl <- list(character(0), "a", "b", c("a", "b"))
  for (q in seq_len(8)) {
    v <- c(sp$outcome[q], "age", ctl[[(q - 1L) %% 4L + 1L]], "t")
    fit <- stats::lm(d[v])
    cf <- summary(fit)$coefficients["t", ]
    expect_equal(sp$est[q], cf[[1L]], tolerance = 1e-10)
    expect_equal(sp$se[q], cf[[2L]], tolerance = 1e-10)
    expect_equal(sp$p[q], cf[[4L]], tolerance = 1e-8)
    expect_equal(c(sp$lo[q], sp$hi[q]),
                 unname(stats::confint(fit)["t", ]), tolerance = 1e-10)
  }
  # the three summaries of the curve, written out
  sg <- if (stats::median(sp$est) >= 0) 1 else -1
  expect_identical(r$sign, sg)
  expect_equal(r$test$obs,
               c(stats::median(sp$est), sum(sp$p < 0.05 & sp$est * sg > 0),
                 sum(sp$stat) / sqrt(8)), tolerance = 1e-12)
  expect_identical(r$test$stat, c("median", "nsig", "z"))
  expect_identical(sp$sig, sp$p < 0.05)
})

test_that("specrv() puts a strong effect beyond every permuted curve", {
  d <- sc_dat(eff = 3, seed = 62)
  r <- specrv(d, "y1", "t", opt = c("a", "b"), nsim = 99, seed = 2)
  expect_identical(r$nok, 99)
  # the median and the summed t are far beyond the permuted ones; the
  # count of significant ones can be matched by chance
  expect_identical(r$test$p[c(1L, 3L)], rep(1 / 100, 2))
  expect_identical(dim(r$dist), c(99L, 3L))
  expect_identical(colnames(r$dist), c("median", "nsig", "z"))
  b <- specrv(d, "y1", "t", opt = "a", null = "boot", nsim = 99, seed = 2)
  expect_identical(b$test$p[c(1L, 3L)], rep(1 / (b$nok + 1), 2))
  expect_gte(b$nok, 90)
  # under the bootstrap null each estimate has its own removed
  expect_lt(abs(stats::median(b$dist[, "median"], na.rm = TRUE)), 0.2)
})

test_that("specrv()'s p-values count the permuted curves at least as extreme", {
  d <- sc_dat(eff = 0.1, seed = 63)
  r <- specrv(d, c("y1", "y2"), "t", opt = "a", nsim = 200, seed = 3)
  sg <- r$sign
  o <- r$test$obs
  k <- c(sum(sg * r$dist[, 1L] >= sg * o[1L] - 1e-9 * abs(o[1L])),
         sum(r$dist[, 2L] >= o[2L]),
         sum(sg * r$dist[, 3L] >= sg * o[3L] - 1e-9 * abs(o[3L])))
  expect_equal(r$test$p, (1 + k) / 201, tolerance = 1e-15)
  expect_identical(specrv(d, c("y1", "y2"), "t", opt = "a", nsim = 200,
                          seed = 3), r)
  set.seed(4)
  a <- specrv(d, "y1", "t", nsim = 10)
  expect_identical(specrv(d, "y1", "t", nsim = 10, seed = a$seed), a)
  out <- utils::capture.output(v <- print(r))
  expect_match(out[1L], "^specrv: 4 specifications, [0-9]+ significant; ")
  expect_identical(v, r)
})

test_that("specrv() refuses specifications it cannot fit", {
  d <- sc_dat()
  expect_error(specrv(d, "y1", "t", x = "a", opt = "a"),
               "specrv: `y`, `d`, `x` and `opt` must name different columns.",
               fixed = TRUE)
  w <- as.data.frame(matrix(rnorm(150 * 11), 150))
  w$t <- d$t
  w$y <- d$y1
  expect_error(specrv(w, "y", "t", opt = paste0("V", 1:11)),
               "specrv: `opt` may name at most ten columns.", fixed = TRUE)
  d$z <- 1
  expect_error(specrv(d, "y1", "z"), "specrv: column `z` must vary.",
               fixed = TRUE)
  expect_error(specrv(d, "y1", "t", opt = "z", seed = 1),
               paste("specrv: the covariates are collinear, or `d` does not",
                     "vary once they are taken out."), fixed = TRUE)
  expect_error(specrv(d, "y1", "t", nsim = 0),
               "specrv: `nsim` must be between 1 and 1e+05.", fixed = TRUE)
})
