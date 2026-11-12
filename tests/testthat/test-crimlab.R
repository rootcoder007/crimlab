# The kernels check their own input, so a caller of the C interface that
# skips the R layer gets an error code, never a read out of bounds or a
# meaningless number. These calls bypass the R checks on purpose.

earg <- "crimlab: a size or an index is out of range"
edom <- "crimlab: a value is outside the domain of the method"

test_that("the disparity kernels refuse what the R layer refuses", {
  expect_error(.Call(C_clab_vod, c(1, 1), c(1, 1), 2L, NULL), earg)
  expect_error(.Call(C_clab_vod, c(1, 1), c(1, 1), -1L, NULL), earg)
  expect_error(.Call(C_clab_vod, c(-1, 1), c(1, 1), 1L, NULL), edom)
  expect_error(.Call(C_clab_vod, c(1, 1), c(1, 0), 0L, NULL), edom)
  expect_error(.Call(C_clab_vod, c(1, 1), c(1, 1), 0L, c(1, 0)), edom)
  expect_error(.Call(C_clab_vod, c(0, 1), c(1, 1), 0L, NULL), edom)
  expect_error(.Call(C_clab_orrisk, 0, 0.5, 0.5), edom)
  expect_error(.Call(C_clab_orrisk, 2, 1, 0.5), edom)
  expect_error(.Call(C_clab_orrisk, 2, 0.5, 0), edom)
  expect_error(.Call(C_clab_outtst, numeric(0), numeric(0), 0.5), earg)
  expect_error(.Call(C_clab_outtst, 0.5, 1, Inf), edom)
  expect_error(.Call(C_clab_outtst, 1.5, 1, 0.5), edom)
  expect_error(.Call(C_clab_outtst, 0.5, -1, 0.5), edom)
  expect_error(.Call(C_clab_outtst, 0.5, 0, 0.5), edom)
  expect_error(.Call(C_clab_dspbm, 1, 1, 1, 1L), earg)
  expect_error(.Call(C_clab_dspbm, 1, 0, 1, 0L), edom)
  one <- 0L
  expect_error(.Call(C_clab_dyadrt, numeric(0), integer(0), integer(0), 1),
               earg)
  expect_error(.Call(C_clab_dyadrt, 1, one, one, numeric(0)), earg)
  expect_error(.Call(C_clab_dyadrt, 1, one, one, 0), edom)
  expect_error(.Call(C_clab_dyadrt, 1, 1L, one, 1), earg)
  expect_error(.Call(C_clab_dyadrt, 1, one, -1L, 1), earg)
  expect_error(.Call(C_clab_dyadrt, -1, one, one, 1), edom)
})

test_that("the forensic kernels refuse what the R layer refuses", {
  expect_error(.Call(C_clab_gtmp, numeric(0), numeric(0), logical(0), 0),
               earg)
  expect_error(.Call(C_clab_gtmp, 0.5, 0, TRUE, 1), edom)
  expect_error(.Call(C_clab_gtmp, 0, 0, TRUE, 0), edom)
  expect_error(.Call(C_clab_gtmp, 0.6, 0.6, FALSE, 0), edom)
  expect_error(.Call(C_clab_gtmp, 0.5, 0, FALSE, 0), edom)
  expect_error(.Call(C_clab_cllr, numeric(0), 1), earg)
  expect_error(.Call(C_clab_cllr, 0, 1), edom)
  expect_error(.Call(C_clab_cllr, 1, -1), edom)
  expect_error(.Call(C_clab_lrcal, 1, TRUE), earg)
  expect_error(.Call(C_clab_lrcal, c(2, 1), c(TRUE, FALSE)), edom)
  expect_error(.Call(C_clab_lrcal, c(0, 1), c(TRUE, FALSE)), edom)
  expect_error(.Call(C_clab_lrcal, c(1, 2), c(TRUE, TRUE)), edom)
  expect_error(.Call(C_clab_lrcal, c(1, 2), c(FALSE, FALSE)), edom)
  expect_error(.Call(C_clab_mixp, numeric(0), logical(0), 1L), earg)
  expect_error(.Call(C_clab_mixp, 0.1, TRUE, -1L), earg)
  expect_error(.Call(C_clab_mixp, rep(0.01, 21), rep(TRUE, 21), 11L), earg)
  expect_error(.Call(C_clab_mixp, 0, TRUE, 1L), edom)
  expect_error(.Call(C_clab_mixp, c(0.6, 0.6), c(TRUE, FALSE), 1L), edom)
  g <- c(0L, 1L)
  f <- c(0.1, 0.2)
  k <- c(0, 1, 0)
  expect_error(.Call(C_clab_kinidx, integer(0), integer(0), numeric(0), k),
               earg)
  expect_error(.Call(C_clab_kinidx, g, g, f, c(-1, 1, 1)), edom)
  expect_error(.Call(C_clab_kinidx, g, g, f, c(0.5, 0.2, 0.2)), edom)
  expect_error(.Call(C_clab_kinidx, g, g, c(0, 0.2), k), edom)
  expect_error(.Call(C_clab_kinidx, c(0L, 0L), g, f, k), edom)
  expect_error(
    .Call(C_clab_patidx, integer(0), integer(0), integer(0), numeric(0)),
    earg
  )
  expect_error(.Call(C_clab_patidx, g, g, g, c(0.1, 1.5)), edom)
  expect_error(.Call(C_clab_patidx, c(0L, 0L), g, g, f), edom)
})

test_that("the policing and corrections kernels refuse bad values", {
  expect_error(.Call(C_clab_clrsmp, -1, 1, 1, 1), edom)
  expect_error(.Call(C_clab_clrsmp, 1, 0, 1, 1), edom)
  expect_error(.Call(C_clab_misrd, 1.5, 0.2, 0.5, 0.9, 0.9), edom)
  expect_error(.Call(C_clab_misrd, 0.5, 0.2, 1, 0.9, 0.9), edom)
  expect_error(.Call(C_clab_misrd, 0.5, 0.2, 0.5, 0.9, -0.1), edom)
  expect_error(.Call(C_clab_rtmhs, c(1, NaN), c(1, 2), c(1, 1), 0.5), edom)
  expect_error(.Call(C_clab_rtmhs, c(1, 2), c(1, 2), c(1, -1), 0.5), edom)
  expect_error(.Call(C_clab_incap, 0, 0.5, 1, NULL), edom)
  expect_error(.Call(C_clab_incap, 2, 1.5, 1, NULL), edom)
  expect_error(.Call(C_clab_incap, 2, 0.5, -1, NULL), edom)
  expect_error(.Call(C_clab_sntbnd, c(0, 0.5), c(0L, 1L), c(1, 1)), edom)
  expect_error(.Call(C_clab_sntbnd, c(0, 1), c(0L, 2L), c(1, 1)), edom)
  expect_error(.Call(C_clab_sntbnd, c(0, 1), c(0L, 1L), c(1, -1)), edom)
  expect_error(.Call(C_clab_hzdsel, 0.5, 0.4, 0.2, c(0.9, 0), c(0.9, 0.8)),
               edom)
  expect_error(.Call(C_clab_hzdsel, 0.5, 0.4, 0.2, c(0.9, 0.8), c(1.2, 0.8)),
               edom)
  expect_error(hzdsel(0.5, 0.4, 0.2, c(0.9, 0.8), c(0.9, 1.2)), "`surl`")
  expect_error(hzdsel(0.5, 0.4, 0.2, c(1.2, 0.8), c(0.9, 0.8)), "`surh`")
})

# The kernels' own tests live in C++ (src/t_*.cpp): hand formulas,
# brute-force enumeration and seeded simulation, checked through crimlab.h.
# This runs them all; a family that stops running drops the count below the
# floor and fails here. The floor sits a little below the count on x86-64
# (13778): a few checks are conditional on the random cases, and the
# smallest family has more than 300.

test_that("the compiled self-test passes every check", {
  r <- .Call(C_clab_selftest)
  expect_identical(r[["report"]], "")
  expect_identical(r[["nfail"]], 0L)
  expect_gte(r[["n"]], 880L)
})

test_that("every exported function answers one valid call in its shape", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  calls <- list(
    list(quote(vodtst(c(w = 6, b = 4), c(5, 5), shift = c(1, 2))),
                  c("group", "share_day", "share_night", "ratio", "sign",
                    "adj")),
    list(quote(orrisk(3, 0.4, 0.5)),
                  c("exposed", "unexposed", "rr", "or")),
    list(quote(outtst(c(0.1, 0.5, 0.9), 0.3)),
                  c("stop", "hit")),
    list(quote(dspbm(c(9, 1), c(3, 1), c(2, 1), ref = 2)),
                 c("group", "cd", "fd", "rd")),
    list(quote(dyadrt(matrix(1, 2, 2), c(a = 1, b = 2))),
                  c("off", "vic", "count", "rate_off", "rate_pair",
                    "null_off", "ratio")),
    list(quote(gtmp(c(0.1, 0.2), c(0.3, NA), c(FALSE, TRUE), 0.01)),
                c("locus", "log10mp", "log10lr")),
    list(quote(cllr(c(5, 9), 0.5)), NULL),
    list(quote(lrpost(c(1, 3), 0.25)), NULL),
    list(quote(kinidx(c(1, 2), c(2, 3), c(0.1, 0.2), "sibling")),
                  c("locus", "log10idx", "post")),
    list(quote(patidx(c(1, 2), c(2, 2), c(0.1, 0.2), c(1, 1))),
                  c("locus", "log10idx", "post")),
    list(quote(patidx(c(1, 2), c(2, 2), c(0.1, 0.2))),
                  c("locus", "log10idx", "post")),
    list(quote(mixlr(c(a = 0.1, b = 0.2), "a", 1, "b", 1)),
                 c("num", "den", "lr")),
    list(quote(dbsrch(1e-6, 1e4, pop = 1e6)),
                  c("np", "any", "post")),
    list(quote(cllrmin(c(3, 0.7), c(0.2, 1))),
                   c("cllr", "min", "cal")),
    list(quote(lrcal(c(3, 0.7), c(0.2, 1))), c("same", "diff")),
    list(quote(tipplt(c(3, 0.7), c(0.2, 1))), c("same", "diff")),
    list(quote(eceplt(c(3, 0.7), c(0.2, 1), lo = 0)),
                  c("lo", "ece", "min", "ref")),
    list(quote(clrsmp(c(9, 2), c(10, 5), c(1, 6), c(2, 14))),
                  c("parts", "total")),
    list(quote(expbnd(c(1, 3), c(10, 10), 1.5, ref = 2)),
                  c("group", "rate", "lo", "hi", "ratio", "rlo", "rhi",
                    "sign")),
    list(quote(misrd(0.3, 0.1, 0.4, 0.8, 0.9)),
                 c("true", "measured", "factor", "d1", "d0")),
    list(quote(rtmhs(c(10, 2, 7), c(6, 3, 8), 6, w = c(2, 1, 1))),
                 c("nsel", "change", "mirror", "sym", "excess", "low",
                   "lowsym")),
    list(quote(incap(c(2, 10), 0.1, 1, share = c(0.8, 0.2))),
                 c("group", "lam", "q", "s", "rate", "prev", "marg")),
    list(quote(incapr(c(5, 4, 3), 0, 2, r = 0.5)),
                  c("prev", "net", "upper", "lower", "later")),
    list(quote(sellab(c(1, NA), c(TRUE, FALSE), c(TRUE, TRUE),
                               w = c(1, 2))),
                  c("obs", "unobs", "lower", "upper", "width", "naive")),
    list(quote(losbnd(c(30, 45), 100)),
                  c("mean", "age", "lower", "bias")),
    list(quote(sntbnd(c(1, 0, 1), c("a", "b", "b"), "a",
                               c(1, 1, 2))),
                  c("term", "share", "lower", "upper", "observed")),
    list(quote(hzdsel(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9))),
                  c("share0", "share1", "haz0", "haz1", "ratio")),
    list(quote(proof("cllr")), c("fn", "theorem", "file"))
  )
  fns <- vapply(calls, function(x) as.character(x[[1L]][[1L]]), "")
  expect_setequal(unique(fns), getNamespaceExports("crimlab"))
  for (k in seq_along(calls)) {
    expect_named(eval(calls[[k]][[1L]]), calls[[k]][[2L]], info = fns[k])
  }
})
