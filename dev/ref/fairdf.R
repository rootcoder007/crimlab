# Reference values of tests/testthat/test-fairdf.R: the disparate impact
# ratios of rmorie::morie_fairness_disparate_impact() on the inputs of
# dev/ref-fair.R (fairin() in the test). From the package root, with
# rmorie in the library:
#   Rscript dev/ref/fairdf.R
source(file.path("dev", "ref", "put.R"))

set.seed(2)
n <- 300
g <- sample(c("A", "B", "C"), n, TRUE)
y <- rbinom(n, 1, 0.4)
p <- rbinom(n, 1, c(A = 0.5, B = 0.35, C = 0.2)[g])
di <- suppressWarnings(
  rmorie::morie_fairness_disparate_impact(p, g, privileged = "A"))
put(list(inp = c(sum(y), sum(p)), di = unlist(di$ratios)), "fairdf",
    "rmorie")
