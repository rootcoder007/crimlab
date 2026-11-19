library(testthat)
library(crimlab)

# The whole suite takes minutes. A check of CRAN runs one file of it,
# test-crimlab.R: the compiled self-test of every kernel, and one valid
# call of every exported function with the shape of its result. With
# NOT_CRAN set (the project's own checks on every push, a developer) the
# whole suite runs.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  test_check("crimlab")
} else {
  test_check("crimlab", filter = "^crimlab$")
}
