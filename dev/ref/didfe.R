# Reference values of tests/testthat/test-didfe.R: MASS::glm.nb() on the
# overdispersed panel the tests build from pepanel() (helper-pe.R). From
# the package root, with MASS in the library:
#   Rscript dev/ref/didfe.R
source(file.path("dev", "ref", "put.R"))
sys.source(file.path("tests", "testthat", "helper-pe.R"),
           envir = environment())

pd <- pepanel()
set.seed(7)
od <- pd
od$cnt <- rnbinom(nrow(od), mu = exp(2 + od$town / 50 + 0.3 * od$on),
                  size = 3)
nb <- MASS::glm.nb(
  cnt ~ on + factor(town) + factor(year) + offset(log(pop)), od
)
put(list(
  inp = sum(od$cnt),
  est = unname(coef(nb)[2]),
  theta = nb$theta
), "didfe", "MASS")
