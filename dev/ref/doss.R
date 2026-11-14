# Reference values of tests/testthat/test-doss.R: spdep::localG() with the
# cell in its own queen neighbourhood, on the grid of counts the test
# builds. From the package root, with spdep in the library:
#   Rscript dev/ref/doss.R
source(file.path("dev", "ref", "put.R"))

set.seed(28)
m <- matrix(rpois(48, 3), 6, 8)
nb <- spdep::include.self(spdep::cell2nb(6, 8, type = "queen"))
lg <- spdep::localG(as.vector(t(m)), spdep::nb2listw(nb, style = "B"))
put(list(inp = sum(m), lg = as.vector(lg)), "doss", "spdep")
