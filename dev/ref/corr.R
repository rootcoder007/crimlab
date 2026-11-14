# Reference values of tests/testthat/test-corr.R: the rmorie functions the
# corrections functions were ported from, on the inputs the test builds.
# From the package root, with rmorie in the library:
#   Rscript dev/ref/corr.R
source(file.path("dev", "ref", "put.R"))

o <- rmorie::morie_incapacitation(c(2, 10), 0.1, 1, shares = c(0.8, 0.2))
inc <- list(
  rate = o[["rate"]],
  prev = o[["prevented_share"]],
  marg = o[["marginal_prevention"]],
  agg = unname(attr(o, "aggregate")[["prevented_share"]])
)
lam <- c(12, 10, 8, 6, 5, 4, 3, 2, 2, 1)
o <- rmorie::morie_incapacitation_career(lam, 2, 3, 0.25)
car <- c(o[["prevented"]], o[["prevented_net"]], o[["upper"]], o[["lower"]],
         o[["later"]])
car3 <- rmorie::morie_incapacitation_career(c(5, 4, 3), 0, 2)[["later"]]

set.seed(4)
risk <- runif(500)
y <- rbinom(500, 1, risk)
o <- rmorie::morie_selective_labels(y, risk < 0.7, risk < 0.8)
sel <- list(
  inp = c(sum(risk), sum(y)),
  obs = o[["observed_rate"]],
  unobs = o[["rule_share_unobserved"]],
  bounds = unname(o[["bounds"]]),
  naive = o[["naive_rate"]]
)

o <- rmorie::morie_backlog_censoring(c(30, 45, 60, 90, 120), c(100, 150, 200))
los <- list(lower = o[["lower_bound"]], bias = o[["bias_lower"]])

set.seed(1)
z <- sample(c("community", "custody"), 500, TRUE, c(0.6, 0.4))
y <- rbinom(500, 1, ifelse(z == "custody", 0.55, 0.35))
o <- rmorie::morie_sentence_effect_bounds(y, z)
snt <- list(
  inp = c(sum(z == "custody"), sum(y)),
  lower = o[["ate_bounds"]][["lower"]],
  upper = o[["ate_bounds"]][["upper"]],
  naive = o[["naive_difference"]]
)

o <- rmorie::morie_hazard_selection(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9))
hz <- list(haz = unname(o[["period2_hazard"]]),
           ratio = o[["period2_hazard_ratio"]])

put(list(inc = inc, car = car, car3 = car3, sel = sel, los = los,
         snt = snt, hz = hz), "corr", "rmorie")
