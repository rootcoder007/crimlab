# Reference values of the OTIS pair presets (test-clprep.R): rmorie's
# treatment-outcome builders morie_otis_make_pair_alert_to_volatility_*()
# and morie_otis_make_pair_a/b/c() on the simulated table opre of
# tests/testthat/helper-pre.R, written to tests/testthat/ref/preset.txt.
# From the package root, with rmorie in the library:
#   Rscript dev/ref/preset.R
need <- "rmorie"
source(file.path("dev", "ref", "put.R"))
sys.source(file.path("tests", "testthat", "helper-pre.R"),
           envir = environment())
f <- function(k) getExportedValue("rmorie", paste0("morie_otis_make_pair_", k))
ru <- f("alert_to_volatility_ruhela")(opre)
nv <- f("alert_to_volatility_naive")(opre)
# the two wrappers add nothing to the builders they call
both <- f("alert_to_volatility_all")(opre)
stopifnot(identical(both, list(ruhela = ru, naive = nv)),
          identical(f("alert_to_volatility_a01")(opre), ru))
# one pair: its treatment, outcome and covariate columns
one <- function(p) {
  c(list(t = p$data[[p$T]], y = p$data[[p$Y]]), as.list(p$data[p$covariates]))
}
put(list(ruhela = one(ru), naive = one(nv), pa = one(f("a")(opre)),
         pb = one(f("b")(opre)), pc = one(f("c")(opre))),
    "preset", need)
