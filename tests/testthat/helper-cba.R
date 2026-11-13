# Units for the design-based causal tests: two covariates, a treatment
# whose probability rises with the first and falls with the second, and an
# effect of 2 on the outcome.
cbdata <- function(seed = 1, n = 300) {
  set.seed(seed)
  x <- data.frame(a = rnorm(n), b = rnorm(n))
  d <- rbinom(n, 1, plogis(0.5 * x$a - 0.3 * x$b))
  y <- 2 * d + x$a - x$b + rnorm(n)
  list(x = x, d = d, y = y)
}
