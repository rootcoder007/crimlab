# Extracted from test-dsc.R:530

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "crimlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
cols <- function(x) {
  x <- unname(as.list(x))
  attributes(x) <- NULL
  x
}
utf8 <- function() {
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "not a UTF-8 session")
  skip_if_not(identical(Sys.getenv("CRIMLAB_ASCII_ONLY"), ""))
}
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
maha_data <- function() {
  set.seed(4)
  m <- data.frame(h = round(rnorm(25, 170, 10), 1))
  m$w <- round(m$h * 0.5 + rnorm(25, 0, 5), 1)
  m$z <- round(runif(25), 3)
  m[1, 1:2] <- c(155, 100)
  m$w[4] <- NA
  m
}

# test -------------------------------------------------------------------------
set.seed(11)
n <- 30
a <- rnorm(n)
m3 <- data.frame(a = a, b = a + rnorm(n), c = rnorm(n))
for (i in 1:n) m3[i, (i %% 3) + 1] <- NA
r <- mcartst(m3)
expect_equal(unname(unlist(r[1:8])),
               c(2.36043941784393, 3, -0.691066760495222, 0.501041292820027,
                 3, 30, 3, 139),
               tolerance = 1e-10)
set.seed(12)
u <- rnorm(10)
v <- rnorm(10)
w <- rnorm(10)
q <- data.frame(x = c(u, rep(NA, 10), w, rnorm(6)),
                  y = c(u, v, rep(NA, 10), rnorm(6)),
                  z = c(rep(NA, 10), v, -w, rnorm(6)))
expect_lt(min(eigen(stats::cov(q, use = "pairwise"))$values), 0)
r <- mcartst(q)
expect_equal(unname(unlist(r[1:8])),
               c(3.64923455237147, 6, -0.322935774282878, 0.724020353601097,
                 4, 36, 3, 96),
               # the start is floored on the scale of the correlations, the
               # reference's on that of the covariances: one optimum, reached
               # to the tolerance of the EM steps
               tolerance = 1e-7)
