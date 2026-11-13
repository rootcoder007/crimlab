# A staggered adoption for the policy-evaluation tests: 30 towns over ten
# years, cohorts first treated in 2018, 2020 and 2022 and towns never
# treated, an effect that grows with time since adoption, a covariate
# that bends the trends, counts, and coarser clusters (ten regions).
pepanel <- function(seed = 2) {
  set.seed(seed)
  d <- expand.grid(year = 2015:2024, town = 1:30)
  d$first <- rep(c(2018, 2020, 2022, NA, NA), 6)[d$town]
  d$on <- as.integer(!is.na(d$first) & d$year >= d$first)
  d$x <- rnorm(30)[d$town] + 0.1 * rnorm(300)
  e <- ifelse(d$on == 1, d$year - d$first, 0)
  d$rate <- 30 + d$town / 3 + 0.5 * (d$year - 2015) +
    d$x * (d$year - 2015) / 5 - (2 + e) * d$on + rnorm(300)
  d$cnt <- rpois(300, exp(1 + d$town / 50 + 0.05 * (d$year - 2015) +
    0.3 * d$on))
  d$pop <- 1000 * (1 + d$town %% 3)
  d$region <- (d$town - 1) %/% 3
  d$w <- 1 + d$town %% 4
  d
}
