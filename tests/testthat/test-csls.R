# The kernels check their own input, so a caller of the C interface that
# skips the R layer gets an error code, never a read out of bounds or a
# meaningless number. These calls bypass the R checks on purpose.

earg <- "csls: a size or an index is out of range"
edom <- "csls: a value is outside the domain of the method"

test_that("the disparity kernels refuse what the R layer refuses", {
  expect_error(.Call(C_csls_vod, c(1, 1), c(1, 1), 2L, NULL), earg)
  expect_error(.Call(C_csls_vod, c(1, 1), c(1, 1), -1L, NULL), earg)
  expect_error(.Call(C_csls_vod, c(-1, 1), c(1, 1), 1L, NULL), edom)
  expect_error(.Call(C_csls_vod, c(1, 1), c(1, 0), 0L, NULL), edom)
  expect_error(.Call(C_csls_vod, c(1, 1), c(1, 1), 0L, c(1, 0)), edom)
  expect_error(.Call(C_csls_vod, c(0, 1), c(1, 1), 0L, NULL), edom)
  expect_error(.Call(C_csls_orrisk, 0, 0.5, 0.5), edom)
  expect_error(.Call(C_csls_orrisk, 2, 1, 0.5), edom)
  expect_error(.Call(C_csls_orrisk, 2, 0.5, 0), edom)
  expect_error(.Call(C_csls_outtst, numeric(0), numeric(0), 0.5), earg)
  expect_error(.Call(C_csls_outtst, 0.5, 1, Inf), edom)
  expect_error(.Call(C_csls_outtst, 1.5, 1, 0.5), edom)
  expect_error(.Call(C_csls_outtst, 0.5, -1, 0.5), edom)
  expect_error(.Call(C_csls_outtst, 0.5, 0, 0.5), edom)
  expect_error(.Call(C_csls_dspbm, 1, 1, 1, 1L), earg)
  expect_error(.Call(C_csls_dspbm, 1, 0, 1, 0L), edom)
  one <- 0L
  expect_error(.Call(C_csls_dyadrt, numeric(0), integer(0), integer(0), 1),
               earg)
  expect_error(.Call(C_csls_dyadrt, 1, one, one, numeric(0)), earg)
  expect_error(.Call(C_csls_dyadrt, 1, one, one, 0), edom)
  expect_error(.Call(C_csls_dyadrt, 1, 1L, one, 1), earg)
  expect_error(.Call(C_csls_dyadrt, 1, one, -1L, 1), earg)
  expect_error(.Call(C_csls_dyadrt, -1, one, one, 1), edom)
})

test_that("the forensic kernels refuse what the R layer refuses", {
  expect_error(.Call(C_csls_gtmp, numeric(0), numeric(0), logical(0), 0),
               earg)
  expect_error(.Call(C_csls_gtmp, 0.5, 0, TRUE, 1), edom)
  expect_error(.Call(C_csls_gtmp, 0, 0, TRUE, 0), edom)
  expect_error(.Call(C_csls_gtmp, 0.6, 0.6, FALSE, 0), edom)
  expect_error(.Call(C_csls_gtmp, 0.5, 0, FALSE, 0), edom)
  expect_error(.Call(C_csls_cllr, numeric(0), 1), earg)
  expect_error(.Call(C_csls_cllr, 0, 1), edom)
  expect_error(.Call(C_csls_cllr, 1, -1), edom)
  expect_error(.Call(C_csls_lrcal, 1, TRUE), earg)
  expect_error(.Call(C_csls_lrcal, c(2, 1), c(TRUE, FALSE)), edom)
  expect_error(.Call(C_csls_lrcal, c(0, 1), c(TRUE, FALSE)), edom)
  expect_error(.Call(C_csls_lrcal, c(1, 2), c(TRUE, TRUE)), edom)
  expect_error(.Call(C_csls_lrcal, c(1, 2), c(FALSE, FALSE)), edom)
  expect_error(.Call(C_csls_mixp, numeric(0), logical(0), 1L), earg)
  expect_error(.Call(C_csls_mixp, 0.1, TRUE, -1L), earg)
  expect_error(.Call(C_csls_mixp, rep(0.01, 21), rep(TRUE, 21), 11L), earg)
  expect_error(.Call(C_csls_mixp, 0, TRUE, 1L), edom)
  expect_error(.Call(C_csls_mixp, c(0.6, 0.6), c(TRUE, FALSE), 1L), edom)
  g <- c(0L, 1L)
  f <- c(0.1, 0.2)
  k <- c(0, 1, 0)
  expect_error(.Call(C_csls_kinidx, integer(0), integer(0), numeric(0), k),
               earg)
  expect_error(.Call(C_csls_kinidx, g, g, f, c(-1, 1, 1)), edom)
  expect_error(.Call(C_csls_kinidx, g, g, f, c(0.5, 0.2, 0.2)), edom)
  expect_error(.Call(C_csls_kinidx, g, g, c(0, 0.2), k), edom)
  expect_error(.Call(C_csls_kinidx, c(0L, 0L), g, f, k), edom)
  expect_error(
    .Call(C_csls_patidx, integer(0), integer(0), integer(0), numeric(0)),
    earg
  )
  expect_error(.Call(C_csls_patidx, g, g, g, c(0.1, 1.5)), edom)
  expect_error(.Call(C_csls_patidx, c(0L, 0L), g, g, f), edom)
})
