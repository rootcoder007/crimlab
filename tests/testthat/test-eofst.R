test_that("eofst() matches the singular value decomposition", {
  set.seed(4)
  z <- matrix(stats::rpois(60, 8), 12, 5) + outer(sin(1:12), 1:5)
  colnames(z) <- paste0("p", 1:5)
  e <- eofst(z)
  a <- sweep(z, 2, colMeans(z))
  s <- svd(a)
  sg <- apply(s$v, 2, function(v) sign(v[which.max(abs(v))]))
  v <- sweep(s$v, 2, sg, `*`)
  expect_identical(dim(e$eofs), c(5L, 5L))
  expect_identical(dimnames(e$eofs), list(colnames(z), paste0("EOF", 1:5)))
  expect_equal(unname(e$eofs), v, tolerance = 1e-8)
  expect_equal(unname(e$pcs), unname(a %*% v), tolerance = 1e-8)
  expect_equal(e$explained, s$d^2 / sum(s$d^2), tolerance = 1e-10)
  expect_equal(e$sdev, stats::prcomp(z)$sdev, tolerance = 1e-10)
  expect_equal(e$center, colMeans(z), tolerance = 1e-14)
  # more places than periods: the decomposition of A A', rank 4 at most
  w <- t(z)
  f <- eofst(w)
  aw <- sweep(w, 2, colMeans(w))
  sw <- svd(aw)
  expect_identical(ncol(f$eofs), 4L)
  vw <- sw$v[, 1:4]
  vw <- sweep(vw, 2, apply(vw, 2, function(u) sign(u[which.max(abs(u))])),
              `*`)
  expect_equal(unname(f$eofs), vw, tolerance = 1e-8)
  expect_equal(f$explained, sw$d[1:4]^2 / sum(sw$d^2), tolerance = 1e-10)
  expect_identical(ncol(eofst(z, k = 2)$pcs), 2L)
  expect_identical(ncol(eofst(as.data.frame(z), k = 99)$eofs), 5L)
  expect_output(print(e), "eofst: 5 components of 12 periods by 5 places")
})

test_that("eofst() finds one pattern written by hand", {
  # columns (1, 3) and (2, 2): centred (-1, 1) and (0, 0), one component
  e <- eofst(matrix(c(1, 3, 2, 2), 2))
  expect_equal(unname(e$eofs[, 1]), c(1, 0), tolerance = 1e-14)
  expect_equal(unname(e$pcs[, 1]), c(-1, 1), tolerance = 1e-14)
  expect_equal(e$explained, 1, tolerance = 1e-14)
  expect_equal(e$sdev, sqrt(2), tolerance = 1e-14)
  expect_null(names(e$center))
})

test_that("eofst() refuses what it cannot decompose", {
  msg <- paste("eofst: `z` must be a numeric matrix of at least two rows:",
               "periods by places.")
  expect_error(eofst(1:3), msg, fixed = TRUE)
  expect_error(eofst(matrix(1:3, 1)), msg, fixed = TRUE)
  expect_error(eofst(matrix("a", 2, 2)), msg, fixed = TRUE)
  expect_error(eofst(data.frame(a = c("x", "y"))), msg, fixed = TRUE)
  expect_error(eofst(matrix(c(1, NA, 2, 3), 2)), "eofst: `z` must be finite.",
               fixed = TRUE)
  expect_error(eofst(matrix(1:4, 2), k = 0),
               "eofst: `k` must be NULL or one whole number from 1.",
               fixed = TRUE)
  expect_error(eofst(matrix(1:4, 2), k = 1.5),
               "eofst: `k` must be NULL or one whole number", fixed = TRUE)
  expect_error(eofst(matrix(1, 3, 2)),
               "eofst: `z` has no variation: every column is constant.",
               fixed = TRUE)
})
