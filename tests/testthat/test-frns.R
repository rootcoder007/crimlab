# Values of the forensic kernels are checked in C++ (src/t_frns.cpp); here
# what R forms itself, the labels and the messages.

test_that("gtmp labels its loci and gives the ratio as the inverse", {
  r <- gtmp(c(a = 0.1, b = 0.2), c(0.3, NA), hom = c(FALSE, TRUE))
  expect_equal(r[["log10lr"]], -r[["log10mp"]], tolerance = 1e-10)
  expect_identical(r[["locus"]][["locus"]], c("a", "b"))
})

test_that("lrpost applies Bayes rule", {
  expect_equal(lrpost(1e6, 1 / 5e6), (1e6 / (5e6 - 1)) / (1 + 1e6 / (5e6 - 1)),
               tolerance = 1e-10)
  expect_equal(lrpost(1, 0.3), 0.3)
  expect_equal(lrpost(c(1, 3), c(0.5, 0.25)), c(0.5, 0.5))
  expect_equal(lrpost(1e308, 0.5), 1, tolerance = 1e-10)
})

test_that("bad forensic input stops with the argument named", {
  errs("gtmp", list(
    list("`p`", 0),
    list("at most 1", 1.2, hom = TRUE),
    list("`q`", 0.2),
    list("`q`", 0.6, 0.6),
    list("`q`", 0.2, "x"),
    list("`q`", c(0.2, 0.2), 0.1),
    list("`hom`", 0.2, hom = NA),
    list("`theta`", 0.2, hom = TRUE, theta = 1)
  ))
  errs("cllr", list(list("`same`", c(1, 0), 1), list("`diff`", 1, -1)))
  errs("lrpost", list(
    list("`lr`", 0, 0.5),
    list("`prior`", 2, 1),
    list("`prior`", 2, c(0.1, 0.2))
  ))
})
