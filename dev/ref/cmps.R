# Reference values from cmpsR's extract_feature_cmps() at its defaults
# (segments of 50, peaks 5, 3, 1, tolerance 25) on simulated land
# signatures: for tests/testthat/test-cmpsc.R, written to
# tests/testthat/ref/cmps.txt with the signatures. Pairs 1 to 4 share a
# source (the reference is a shifted, noisier copy); pairs 5 to 8 do not.
# From the package root, with cmpsR in the library:
#   Rscript dev/ref/cmps.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- "cmpsR"
source(file.path("dev", "ref", "put.R"))
set.seed(52)
# a land signature: striae as a smoothed random walk plus fine marks,
# missing points at the ends as an extraction leaves them
sig <- function(n) {
  w <- stats::filter(cumsum(rnorm(n + 40, sd = 0.6)), rep(1 / 9, 9), sides = 2)
  as.numeric(w[21:(n + 20)] - mean(w, na.rm = TRUE)) + rnorm(n, sd = 0.3)
}
cases <- list()
for (k in 1:8) {
  n <- sample(900:1100, 1)
  a <- sig(n)
  if (k <= 4) {
    sh <- sample(-60:60, 1)
    i <- seq(max(1, 1 + sh), min(n, n + sh))
    b <- a[i] + rnorm(length(i), sd = 0.35)
  } else {
    b <- sig(sample(900:1100, 1))
  }
  a[1:sample(0:15, 1)] <- NA
  b[length(b) - sample(0:15, 1):0] <- NA
  r <- cmpsR::extract_feature_cmps(a, b, include = "full_result")
  cases[[k]] <- list(
    a = a,
    b = b,
    score = r$CMPS_score,
    nseg = r$nseg,
    pos = r$congruent_pos
  )
}
names(cases) <- paste0("p", seq_along(cases))
put(cases, "cmps", need)
