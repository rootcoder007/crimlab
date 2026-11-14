# Reference values of tests/testthat/test-ctabx.R: DescTools::CramerV(),
# with its interval and bias-corrected, on the table the tests build.
# From the package root, with DescTools in the library:
#   Rscript dev/ref/ctabx.R
source(file.path("dev", "ref", "put.R"))

set.seed(7)
ct <- data.frame(
  a = sample(c("x", "y", "z"), 120, TRUE),
  b = sample(c("p", "q", "r", "s"), 120, TRUE)
)
tb <- table(ct$a, ct$b)
put(list(
  inp = as.vector(tb),
  v = unname(DescTools::CramerV(tb, conf.level = 0.95)),
  vbc = DescTools::CramerV(tb, correct = TRUE)
), "ctabx", "DescTools")
