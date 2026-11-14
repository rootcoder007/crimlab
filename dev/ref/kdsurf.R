# Reference values of tests/testthat/test-kdsurf.R: MASS::bandwidth.nrd()
# of the points the test builds. From the package root, with MASS in the
# library:
#   Rscript dev/ref/kdsurf.R
source(file.path("dev", "ref", "put.R"))

set.seed(61)
p <- cbind(x = runif(30, 0, 4), y = runif(30, 0, 3))
put(list(
  inp = sum(p),
  bw = c(MASS::bandwidth.nrd(p[, 1]), MASS::bandwidth.nrd(p[, 2]))
), "kdsurf", "MASS")
