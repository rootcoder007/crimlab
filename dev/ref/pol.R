# Reference values of tests/testthat/test-pol.R: the rmorie functions the
# policing functions were ported from, on the inputs the test builds.
# From the package root, with rmorie in the library:
#   Rscript dev/ref/pol.R
source(file.path("dev", "ref", "put.R"))

o <- rmorie::morie_disparity_exposure_bounds(
  c(A = 300, B = 100, C = 50), c(A = 1000, B = 1000, C = 400),
  gamma = 1.5, reference = "B"
)
exb <- list(rate = o[["proxy_rate"]], lo = o[["rate_lower"]],
            hi = o[["rate_upper"]], rlo = o[["ratio_lower"]],
            rhi = o[["ratio_upper"]], dir = o[["direction_identified"]])
o <- rmorie::morie_misclassified_rd(0.3, 0.1, 0.4, 0.8, 0.9)
mis <- list(measured = o[["measured_rd"]], factor = o[["factor"]])
set.seed(5)
mu <- rgamma(300, 2, 0.5)
x1 <- rpois(300, mu)
x2 <- rpois(300, mu)
cut <- quantile(x1, 0.8, names = FALSE)
o <- rmorie::morie_regression_to_mean(x1, x2, cut)
rtm <- list(inp = c(sum(mu), sum(x1), sum(x2)),
            change = o[["selected_change"]], mirror = o[["mirror_change"]],
            sym = o[["symmetrised_change"]],
            low = o[["low_selected_change"]],
            lowsym = o[["low_symmetrised_change"]])
put(list(exb = exb, mis = mis, rtm = rtm), "pol", "rmorie")
