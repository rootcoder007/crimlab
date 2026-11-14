# Reference values from cmpsR's extract_feature_cmps() at its defaults
# (segments of 50, peaks 5, 3, 1, tolerance 25) on simulated land
# signatures: for tests/testthat/test-cmpsc.R, written to
# tests/testthat/ref/cmps.txt with the signatures: the score, the number of
# segments, the congruent position and every segment's candidates. Pairs 1
# to 4 share a source (the reference is a shifted, noisier copy); pairs 5
# to 8 do not; then edge cases (different lengths, missing runs at the
# ends and inside, a flat stretch, ties, one level of peaks) and cases
# cmpsR refuses (err TRUE).
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
run <- function(a, b, npeaks = c(5, 3, 1)) {
  r <- tryCatch(cmpsR::extract_feature_cmps(a, b, npeaks_set = npeaks, include = "full_result"),
                error = function(e) NULL)
  if (is.null(r)) return(list(a = a, b = b, npeaks = npeaks, err = TRUE))
  list(a = a, b = b, npeaks = npeaks, err = FALSE, score = r$CMPS_score, nseg = r$nseg,
       pos = r$congruent_pos, ccp = lapply(r$ccp_list, function(z) if (is.null(z)) numeric(0) else z))
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
  cases[[k]] <- run(a, b)
}
names(cases) <- paste0("p", seq_along(cases))
a <- sig(700)
cases$len <- run(a, c(a[101:700], sig(600)))
a <- sig(1000)
b <- a[21:1000] + rnorm(980, sd = 0.3)
a[c(1:25, 961:1000)] <- NA
b[c(1:30, 971:980)] <- NA
cases$ends <- run(a, b)
a <- sig(900)
b <- a[11:900]
a[300:419] <- 0.5
cases$flat <- run(a, b)
a <- sig(900)
b <- a[31:900] + rnorm(870, sd = 0.2)
a[401:500] <- NA
cases$hole <- run(a, b)
pat <- sin(seq_len(50) / 4) + rnorm(50, sd = 0.3)
a <- rep(pat, 12)
cases$tie <- run(a, rep(pat, 14)[-(1:7)])
a <- sig(800)
cases$one <- run(a, a[41:800] + rnorm(760, sd = 0.3), npeaks = 5)
cases$short <- run(sig(40), sig(800))
cases$shortb <- run(sig(800), sig(150))
put(cases, "cmps", need)
