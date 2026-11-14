# Reference values of tests/testthat/test-disp.R: the rmorie functions the
# disparity functions were ported from, on the inputs the test builds.
# From the package root, with rmorie in the library:
#   Rscript dev/ref/disp.R
source(file.path("dev", "ref", "put.R"))

m <- rmorie::morie_veil_of_darkness(c(w = 600, b = 400), c(w = 700, b = 300))
vod <- list(ratio = m[["bias_ratio"]], shday = m[["share_day"]])
orr <- unname(rmorie::morie_relative_risk_from_or(3, 0.4, 0.5)[["risks"]])
g <- list(A = list(q = c(0.3, 0.9), w = c(1, 1), t = 0.5),
          B = list(q = 0.9, w = 1, t = 0.2))
m <- rmorie::morie_outcome_test(g)
out <- list(stop = m[["stop_rate"]], hit = m[["hit_rate"]])
m <- rmorie::morie_disparity_benchmark(c(a = 9000, b = 1000),
                                       c(a = 900, b = 300), c(a = 9, b = 6),
                                       reference = "a")
bm <- list(cd = m[["contact_disparity"]],
           fd = m[["force_given_contact_disparity"]],
           rd = m[["resident_disparity"]])
m <- rmorie::morie_interracial_rates(
  c(A_on_A = 900, B_on_A = 200, A_on_B = 120, B_on_B = 300),
  c(A = 80000, B = 20000)
)
dy <- list(off = m[["offender"]], vic = m[["victim"]],
           rtoff = m[["rate_per_offender_group"]],
           rtpair = m[["rate_per_pair_exposure"]],
           nulloff = m[["null_rate_per_offender_group"]],
           ratio = m[["ratio_to_null"]])
put(list(vod = vod, orr = orr, out = out, bm = bm, dy = dy), "disp",
    "rmorie")
