#' @srrstats {G2.2} A matrix with more than one row and column is refused
#'   for each kind of vector argument; a one-column matrix is read as its
#'   vector.
test_that("a vector argument refuses a matrix", {
  m <- matrix(c(6, 4, 5, 5), 2)
  v <- "must be a vector, not a matrix"
  errs("vodtst", list(list(v, m, c(5, 5, 4, 6))))
  errs("sellab", list(list(v, c(1, 0), m > 4, c(TRUE, TRUE)),
                      list(v, matrix(c(1, 0, 1, 0), 2), rep(TRUE, 4),
                           rep(TRUE, 4))))
  errs("gtmp", list(list(v, c(0.1, 0.2, 0.1, 0.2), m / 20),
                    list(v, c(0.1, 0.2, 0.1, 0.2), hom = matrix(TRUE, 2, 2))))
  errs("mixlr", list(list(v, c(a = 0.1, b = 0.2), matrix("a", 2, 2), 1,
                          "a", 2)))
  errs("sntbnd", list(list(v, c(1, 0, 1, 0), matrix(c("a", "b"), 2, 2))))
  expect_identical(vodtst(matrix(c(6, 4)), c(5, 5)), vodtst(c(6, 4), c(5, 5)))
})
