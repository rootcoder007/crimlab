fc_dat <- function(n = 120, seed = 41) {
  set.seed(seed)
  d <- data.frame(age = runif(n, 18, 60), k = rbinom(n, 1, 0.5),
                  u = rnorm(n))
  d$t <- rbinom(n, 1, plogis(0.8 * d$u))
  d$y <- 0.4 * d$t + 0.5 * d$u + 0.01 * d$age + rnorm(n)
  d$c1 <- 0.5 * d$u + rnorm(n)
  d$c2 <- rnorm(n)
  d$c3 <- 2 * d$t + rnorm(n, sd = 0.1)
  d
}

test_that("falsctl() gives lm()'s coefficient, classical error and interval", {
  d <- fc_dat()
  r <- falsctl(d, "y", "t", c("c1", "c2"), x = c("age", "k"), se = "iid")
  tab <- r$tab
  expect_identical(tab$outcome, c("y", "c1", "c2"))
  expect_identical(tab$role, c("primary", "control", "control"))
  for (j in 1:3) {
    v <- tab$outcome[j]
    fit <- stats::lm(d[[v]] ~ d$age + d$k + d$t)
    cf <- summary(fit)$coefficients["d$t", ]
    ci <- stats::confint(fit)["d$t", ]
    expect_equal(tab$est[j], cf[[1L]], tolerance = 1e-10)
    expect_equal(tab$se[j], cf[[2L]], tolerance = 1e-10)
    expect_equal(tab$stat[j], cf[[3L]], tolerance = 1e-10)
    expect_equal(tab$p[j], cf[[4L]], tolerance = 1e-8)
    expect_equal(c(tab$lo[j], tab$hi[j]), unname(ci), tolerance = 1e-10)
  }
  expect_equal(tab$padj[2:3], stats::p.adjust(tab$p[2:3], "holm"),
               tolerance = 1e-15)
  expect_true(is.na(tab$padj[1L]))
  expect_true(all(is.na(tab$pperm)))
  expect_identical(r$nflag, sum(tab$padj <= 0.05, na.rm = TRUE))
  expect_identical(r$seed, NA_real_)
})

test_that("falsctl()'s HC1 error is the sandwich written out", {
  d <- fc_dat(seed = 42)
  r <- falsctl(d, NULL, "t", c("c1", "c2"), x = "age", adj = "BH")
  X <- cbind(1, d$age, d$t)
  br <- solve(crossprod(X))
  for (j in 1:2) {
    v <- c("c1", "c2")[j]
    e <- d[[v]] - X %*% (br %*% crossprod(X, d[[v]]))
    meat <- crossprod(X * as.vector(e))
    vc <- nrow(X) / (nrow(X) - 3) * br %*% meat %*% br
    expect_equal(r$tab$se[j], sqrt(vc[3, 3]), tolerance = 1e-10)
  }
  expect_identical(r$tab$role, c("control", "control"))
  expect_equal(r$tab$padj, stats::p.adjust(r$tab$p, "BH"), tolerance = 1e-15)
})

test_that("falsctl() flags a control the treatment drives", {
  d <- fc_dat(seed = 43)
  r <- falsctl(d, "y", "t", c("c2", "c3"), nperm = 199, seed = 6)
  # c3 is twice the treatment plus little noise: no permutation reaches it
  expect_identical(r$tab$pperm[3L], 1 / 200)
  expect_true(all(r$tab$pperm >= 1 / 200 & r$tab$pperm <= 1))
  expect_equal(r$tab$padj[2:3], stats::p.adjust(r$tab$pperm[2:3], "holm"),
               tolerance = 1e-15)
  expect_gte(r$nflag, 1L)
  expect_identical(falsctl(d, "y", "t", c("c2", "c3"), nperm = 199,
                           seed = 6), r)
  set.seed(1)
  a <- falsctl(d, "y", "t", "c2", nperm = 20)
  expect_identical(falsctl(d, "y", "t", "c2", nperm = 20, seed = a$seed), a)
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], paste0(
    "falsctl: ", r$nflag, " of 2 negative control(s) associated with the ",
    "treatment at alpha 0.05 (holm)"
  ))
  expect_identical(v, r)
})

test_that("falsctl() refuses columns and settings it cannot use", {
  d <- fc_dat()
  expect_error(falsctl(d, "y", "t", "y"),
               "falsctl: `y`, `d`, `nc` and `x` must name different columns.",
               fixed = TRUE)
  d$z <- 1
  expect_error(falsctl(d, "y", "z", "c1"), "falsctl: column `z` must vary.",
               fixed = TRUE)
  expect_error(falsctl(d, "y", "t", "c1", x = c("age", "age2")),
               "falsctl: `x` names columns not in `data`: age2.",
               fixed = TRUE)
  d$age2 <- 2 * d$age
  expect_error(falsctl(d, "y", "t", "c1", x = c("age", "age2")),
               paste("falsctl: the covariates are collinear, or `d` does not",
                     "vary once they are taken out."), fixed = TRUE)
  expect_error(falsctl(d[1:3, ], "y", "t", "c1", x = "age"),
               "falsctl: `data` needs more rows than covariates plus two.",
               fixed = TRUE)
  expect_error(falsctl(d, "y", "t", "c1", nperm = -1),
               "falsctl: `nperm` must be between 0 and 1e+06.", fixed = TRUE)
  expect_error(falsctl(d, "y", "t", "c1", adj = "fdr"),
               paste("falsctl: `adj` must be one of \"holm\", \"bonferroni\",",
                     "\"hochberg\", \"BH\", \"BY\", \"none\", not \"fdr\"."),
               fixed = TRUE)
})
