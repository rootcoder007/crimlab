# A stratified cluster sample for the survey tests: three regions, six
# blocks of ten respondents in each, unequal weights.
svdata <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(region = rep(c("n", "s", "e"), each = 60),
                  block = rep(1:18, each = 10))
  d$v <- rbinom(180, 1, 0.3)
  d$y <- rpois(180, 3)
  d$x <- rpois(180, 5) + 1
  d$w <- runif(180, 50, 150)
  d$N <- c(n = 40, s = 30, e = 50)[d$region]
  d
}
