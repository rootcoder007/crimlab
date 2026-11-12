# Old-versus-new comparison of the exported functions: the same calls, on
# random and edge inputs, run against two installed builds; values,
# attributes (in order), error and warning messages, classes and calls.
#   Rscript dev/cmp.R run LIB OUT.rds    evaluate every case with crimlab
#                                        from LIB, save the outcomes
#   Rscript dev/cmp.R diff OLD.rds NEW.rds
# Prints one line per case that is not identical, then a summary line
# "CMP identical N, within 1e-12 M, different K".
a <- commandArgs(TRUE)

cases <- function() {
  set.seed(20261112)
  r <- function(n, lo = 0.01, hi = 1) stats::runif(n, lo, hi)
  big <- r(1000, 1e-3, 1e3)
  lst <- list(
    # disparity
    quote(vodtst(c(w = 60, b = 40), c(w = 50, b = 50))),
    quote(vodtst(c(w = 60, b = 40), c(w = 50, b = 50), shift = c(1, 2 / 3))),
    quote(vodtst(c(w = 60, b = 40), c(w = 50, b = 50), ref = "b")),
    quote(vodtst(c(5L, 5L), c(7L, 7L))),
    quote(vodtst(c(1, NA), c(1, 1))), quote(vodtst(c(1, 1), c(1, 0))),
    quote(vodtst(c(1, 1), c(1, 1, 1))), quote(vodtst(c(0, 1), c(1, 1))),
    quote(vodtst(c(1, 1), c(1, 1), ref = 3)),
    quote(vodtst(c(1, 1), c(1, 1), ref = 1.5)),
    quote(vodtst(c(1, 1), c(1, 1), ref = NA_character_)),
    quote(vodtst(c(1, 1), c(1, 1), ref = TRUE)),
    quote(vodtst(c(1, 1), c(1, 1), ref = 2L)),
    quote(vodtst(c(a = 1, b = 1), c(1, 1), ref = "z")),
    quote(vodtst(c(a = 1, b = 1), c(1, 1), ref = c("a", "b"))),
    quote(vodtst(c(1, 1), c(1, 1), shift = c(1, -1))),
    quote(vodtst(factor(c(1, 2)), c(1, 1))),
    quote(vodtst(Sys.Date() + 0:1, c(1, 1))),
    quote(vodtst(structure(c(3, 4), class = "foo"), c(1, 1))),
    quote(vodtst(array(c(3, 4), 2, list(c("x", "y"))), c(1, 1))),
    quote(vodtst("a", 1)), quote(vodtst(numeric(0), numeric(0))),
    quote(vodtst(c(1, Inf), c(1, 1))), quote(vodtst(c(1, NaN), c(1, 1))),
    quote(vodtst(list(1, 2), c(1, 1))), quote(vodtst(NULL, 1)),
    bquote(vodtst(.(r(50, 1, 100)), .(r(50, 1, 100)), shift = .(r(50)))),
    quote(orrisk(3, 0.4, 0.5)), quote(orrisk(0, 0.5, 0.5)),
    quote(orrisk(2, 1, 0.5)), quote(orrisk(2, 0.5, 0)),
    quote(orrisk(1e-320, 0.5, 0.5)), quote(orrisk(c(a = 2), 0.3, 0.2)),
    quote(orrisk(2L, 0.3, 0.2)), quote(orrisk(c(2, 3), 0.3, 0.2)),
    quote(orrisk(2, NA, 0.2)), quote(orrisk(2, 0.3, c(0.2, 0.1))),
    quote(orrisk(2, "a", 0.2)), quote(orrisk(1e-17, 0.99, 0.98)),
    quote(outtst(c(0.1, 0.5, 0.9), 0.3)),
    quote(outtst(c(0.1, 0.5, 0.9), 0.3, w = c(1, 3, 1))),
    quote(outtst(c(0.1, 0.2), 0.5)), quote(outtst(c(0.5, 1.5), 0.2)),
    quote(outtst(c(0.5, 0.6), 0.2, w = c(0, 0))), quote(outtst("a", 0.2)),
    quote(outtst(0.5, c(1, 2))), quote(outtst(0.5, 0.2, w = c(1, 2))),
    bquote(outtst(.(r(300)), 0.4, w = .(r(300)))),
    quote(dspbm(c(a = 9000, b = 1000), c(a = 900, b = 300), c(a = 9, b = 6))),
    quote(dspbm(c(1, 1), c(1, 0), c(1, 1))),
    quote(dspbm(c(1, 1), c(1, 1), c(1, 1), ref = 0)),
    quote(dyadrt(outer(c(a = 800, b = 150, c = 50), c(800, 150, 50)) * 0.002,
                 c(a = 800, b = 150, c = 50))),
    quote(dyadrt(matrix(c(900, 200, 120, 300), 2), c(A = 80000, B = 20000))),
    quote(dyadrt(matrix(1, 2, 3), c(1, 1))),
    quote(dyadrt(matrix(0, 2, 2), c(1, 1))),
    quote(dyadrt(matrix(-1, 2, 2), c(1, 1))),
    quote(dyadrt(matrix("a", 2, 2), c(1, 1))),
    quote(dyadrt(matrix(1L, 2, 2), c(1, 1))),
    quote(dyadrt(structure(matrix(1, 2, 2), class = "table"), c(1, 1))),
    quote(dyadrt(as.difftime(matrix(1, 2, 2), units = "days"), c(1, 1))),
    quote(dyadrt(structure(matrix(1, 2, 2), class = "difftime",
                           units = "days"), c(1, 1))),
    quote(dyadrt(structure(matrix(1L, 2, 2), class = "foo"), c(1, 1))),
    quote(gtmp(0.2, new.env(), hom = TRUE)),
    quote(gtmp(0.2, function() 1, hom = TRUE)),
    quote(gtmp(c(0.2, 0.3), c(NA, 1i), hom = TRUE)),
    quote(gtmp(c(0.2, 0.3), c(NA, NA_complex_), hom = TRUE)),
    quote(gtmp(c(0.2, 0.3), as.raw(1:2), hom = TRUE)),
    quote(gtmp(c(0.2, 0.3), factor(c(NA, NA)), hom = TRUE)),
    quote(gtmp(c(0.2, 0.3), as.Date(c(NA, NA)), hom = TRUE)),
    quote(gtmp(c(0.2, 0.3), list(NA, 1), hom = TRUE)),
    quote(eceplt(2, 0.5, lo = -400)),
    quote(kinidx(factor(c("1", "2")), c(1, 2), c(0.1, 0.2))),
    quote(kinidx(c(TRUE, FALSE), c(1, 2), c(0.1, 0.2))),
    quote(sntbnd(c(0, 1), c("a", "B"))), quote(sntbnd(c(0, 1), c(2, 10))),
    quote(sntbnd(c(0, 1), list("a", "b"))),
    quote(sntbnd(c(0, 1), list("a", NA))),
    quote(mixlr(stats::setNames(c(0.1, 0.2), c("a", NA)), "a", 1, "a", 1)),
    quote(dyadrt(matrix(factor(1:4), 2, 2), c(1, 1))),
    quote(dyadrt(data.frame(a = 1:2, b = 1:2), c(1, 1))),
    quote(dyadrt(1:4, c(1, 1))),
    # forensics
    quote(gtmp(c(a = 0.1, b = 0.2), c(0.3, NA), hom = c(FALSE, TRUE))),
    quote(gtmp(c(0.1, 0.1), c(NA, 0.3), hom = c(TRUE, FALSE),
               theta = 0.03)),
    quote(gtmp(rep(1e-3, 400), hom = TRUE)), quote(gtmp(0.5, 0.5)),
    quote(gtmp(0)), quote(gtmp(1.2, hom = TRUE)), quote(gtmp(0.2)),
    quote(gtmp(0.6, 0.6)), quote(gtmp(0.2, "x")), quote(gtmp(c(0.2, 0.2), 0.1)),
    quote(gtmp(0.2, hom = NA)), quote(gtmp(0.2, hom = TRUE, theta = 1)),
    quote(gtmp(0.2, hom = TRUE, theta = -0.1)),
    quote(gtmp(0.2, hom = TRUE, theta = c(0.1, 0.2))),
    quote(gtmp(0.2, NA, hom = TRUE)), quote(gtmp(0.2, NA_character_, TRUE)),
    quote(gtmp(c(0.2, 0.3), list(NA, NA), TRUE)),
    quote(gtmp(c(0.2, 0.3), c(0.1, 0.2), hom = c(TRUE, FALSE, TRUE))),
    quote(gtmp(c(0.2, 0.3), c(0.1, 0.2), hom = 1)),
    quote(gtmp(0.2, Inf)), quote(gtmp(0.2, 2L)),
    quote(gtmp(c(0.2, 0.3), c(0.1, NaN), hom = c(FALSE, TRUE))),
    quote(cllr(c(50, 200, 8), c(0.01, 0.2, 1.5))), quote(cllr(c(1, 0), 1)),
    quote(cllr(1, -1)), quote(cllr(1e-310, 1)), quote(cllr(1, 1e300)),
    bquote(cllr(.(big), .(1 / big))),
    quote(lrpost(1e6, 1 / 5e6)), quote(lrpost(c(1, 3), c(0.5, 0.25))),
    quote(lrpost(1e308, 0.5)), quote(lrpost(0, 0.5)), quote(lrpost(2, 1)),
    quote(lrpost(2, c(0.1, 0.2))), quote(lrpost(c(a = 2, b = 3), 0.1)),
    quote(lrpost(1e300, 1e-300)),
    quote(mixlr(c(a = 0.1, b = 0.2, c = 0.3), c("a", "b"), 1, "a", 2)),
    quote(mixlr(c(a = 0.1, b = 0.2), c("a", "b"), 0, c("a", "a"), 1)),
    quote(mixlr(c(0.1, 0.2), "a", 1, "a", 1)),
    quote(mixlr(c(a = 0.7, b = 0.6), "a", 1, "a", 1)),
    quote(mixlr(c(a = 0.1, a = 0.2), "a", 1, "a", 1)),
    quote(mixlr(c(a = 0.1), "z", 1, "a", 1)),
    quote(mixlr(c(a = 0.1), "a", 1, "a", 1.5)),
    quote(mixlr(c(a = 0.1), "a", 51, "a", 1)),
    quote(mixlr(c(a = 0.1), NA_character_, 1, "a", 1)),
    quote(mixlr(c(a = 0.1), factor("a"), 1, "a", 1)),
    quote(mixlr(c(a = 0.5, b = 0.5 + 1e-9), "a", 1, "a", 1)),
    quote(mixlr(setNames(rep(0.01, 21), letters[1:21]), letters[1:21], 11,
                "a", 1)),
    quote(dbsrch(1e-6, 1e5)), quote(dbsrch(1e-6, 1e5, pop = 3e7)),
    quote(dbsrch(0.5, 10)), quote(dbsrch(2, 10)), quote(dbsrch(0.1, 0)),
    quote(dbsrch(0.1, 10, pop = 5)), quote(dbsrch(0.1, 100000, pop = 5)),
    quote(dbsrch(0.1, 1.5, pop = 1)), quote(dbsrch(0.1, 4 / 3, pop = 1)),
    quote(dbsrch(0.1, 123456789012, pop = 5)),
    quote(cllrmin(c(50, 200, 8), c(0.01, 0.2, 1.5))),
    quote(cllrmin(c(1, -1), 1)), quote(lrcal(1, numeric(0))),
    quote(lrcal(c(5, 0.5, 5, 2), c(0.1, 5, 1, 1 + 2e-16))),
    bquote(lrcal(.(big[1:500]), .(big[501:1000]))),
    bquote(cllrmin(.(big[1:500]), .(big[501:1000]))),
    quote(lrcal(c(1, 1, 1), c(1, 1))),
    # kinship
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2))),
    quote(kinidx(rbind(D3 = c(14, 15), vWA = c(16, 16)),
                 rbind(c(14, 16), c(16, 17)), rbind(c(0.1, 0.2), c(0.3, 0.3)),
                 "sibling")),
    quote(kinidx(c("a", "b"), c("a", "c"), c(0.1, 0.2), "half")),
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), c(0.2, 0.3, 0.5))),
    quote(kinidx(1:3, 1:2, c(0.1, 0.2))),
    quote(kinidx(c(1, 2), matrix(1, 2, 2), c(0.1, 0.2))),
    quote(kinidx(c(1, 2), c(1, NA), c(0.1, 0.2))),
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 1.2))),
    quote(kinidx(c(1, 1), c(1, 2), c(0.1, 0.2))),
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), "twin")),
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), c(0.5, 0.2, 0.2))),
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), NA_character_)),
    quote(kinidx(factor(c("a", "b")), c("a", "b"), c(0.1, 0.2))),
    quote(kinidx(Sys.Date() + 0:1, Sys.Date() + 0:1, c(0.1, 0.2))),
    quote(kinidx(c(0.1 + 0.2, 0.3), c(0.3, 0.3), c(0.1, 0.2))),
    quote(kinidx(matrix(1:2, 1, dimnames = list("r", NULL)), c(1, 2),
                 c(0.1, 0.2))),
    quote(kinidx(matrix(1:2, 1, dimnames = list(NULL, c("x", "y"))),
                 c(1, 2), c(0.1, 0.2))),
    quote(kinidx(matrix(1:2, 1, dimnames = list("r", c("x", "y"))),
                 c(1, 2), c(0.1, 0.2))),
    quote(kinidx(matrix(1, 0, 2), matrix(1, 0, 2), matrix(0.1, 0, 2))),
    quote(kinidx(matrix(TRUE, 1, 2), c(1, 2), c(0.1, 0.2))),
    quote(kinidx(list(1, 2), c(1, 2), c(0.1, 0.2))),
    quote(kinidx(c(1, 2), c(1, 2), matrix(c(0.1, 0.2), 2, 1))),
    quote(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2), 1:3)),
    quote(patidx(c(1, 2), c(1, 2), c(0.1, 0.2))),
    quote(patidx(c(1, 2), c(1, 2), c(0.1, 0.2), mother = c(1, 3))),
    quote(patidx(c(1, 2), c(1, 2), c(0.1, 0.2), mother = c(3, 3))),
    quote(patidx(c(1, 2), c(1, 2), c(0.1, 0.2), mother = 1)),
    quote(patidx(c(1, 2), c(1, 2), c(0.1, 1.01), mother = c(1, 1))),
    quote(patidx(c(1, 2), 1:3, c(0.1, 0.2))),
    quote(patidx(matrix(1, 0, 2), matrix(1, 0, 2), matrix(0.1, 0, 2),
                 matrix(1, 0, 2))),
    quote(tipplt(c(50, 200, 8), c(0.01, 0.2, 1.5))), quote(tipplt(0, 1)),
    quote(eceplt(c(50, 200, 8), c(0.01, 0.2, 1.5))),
    quote(eceplt(1e-310, 1, lo = 0)), quote(eceplt(1, 1, lo = NA)),
    quote(eceplt(c(2, 3), 0.5, lo = c(-400, 2, 400))),
    quote(eceplt(2, 0.5, lo = 400)), quote(eceplt(c(2, 3), 0.5, lo = -400)),
    quote(eceplt(c(0.5, 2), c(3, 0.1), lo = c(-400, 400))),
    quote(eceplt(c(0.5, 2, 9), c(3, 0.1), lo = 400)),
    bquote(eceplt(.(big[1:300]), .(big[301:600]))),
    # policing
    quote(clrsmp(c(minor = 50, serious = 30), c(100, 50), c(40, 35),
                 c(60, 70))),
    quote(clrsmp(c(1, 1), c(1, 0), c(1, 1), c(1, 1))),
    quote(clrsmp(c(1, 1), c(1, 1), c(1, 1), "a")),
    quote(expbnd(c(a = 2, b = 3), c(1, 1), 1.5)),
    quote(expbnd(c(2, 3), c(1, 1), 1, ref = 2)),
    quote(expbnd(c(1, 1), c(1, 1), 0.9)), quote(expbnd(c(1, 1), c(1, 1), c(1, 2))),
    quote(misrd(0.3, 0.1, 0.4, 0.9, 0.8)), quote(misrd(0.1, 0.1, 0.5, 0, 1)),
    quote(misrd(1.1, 0.1, 0.5, 0.9, 0.9)), quote(misrd(0.1, NA, 0.5, 0.9, 0.9)),
    quote(misrd(0.1, 0.1, 1, 0.9, 0.9)), quote(misrd(0.1, 0.1, 0.5, 0.9, c(1, 1))),
    quote(rtmhs(c(5, 1, 7, 3), c(4, 2, 5, 3), 4)),
    quote(rtmhs(c(5, 1), c(1, 1), 4)), quote(rtmhs(c(1, 2), c(1, 2), 5)),
    quote(rtmhs(c(1, 2), c(1, 2), 1, w = c(1, 0))),
    quote(rtmhs(c(1, Inf), c(1, 2), 0)), quote(rtmhs(c(1, 2), c(1, 2), NA)),
    bquote(rtmhs(.(r(500, 0, 10)), .(r(500, 0, 10)), 5, w = .(r(500)))),
    # corrections
    quote(incap(c(a = 1, b = 2, c = 3), 0.5, c(1, 2, 3))),
    quote(incap(c(low = 1, high = 4), c(0.3, 0.6), 2, share = c(0.7, 0.3))),
    quote(incap(c(1L, 2L), 0.5, 1)), quote(incap(0, 0.1, 1)),
    quote(incap(1, 1.5, 1)), quote(incap(c(1, 2, 3), c(0.1, 0.2), 1)),
    quote(incap(c(1, 2, 3), 0.1, c(1, 2))),
    quote(incap(c(1, 2), 0.1, 1, share = c(0.5, 0.6))),
    quote(incap(c(1, 2), 0.1, 1, share = c(1.5, -0.5))),
    quote(incapr(c(3, 2, 1, 0.5), 1, 2, r = 0.2)),
    quote(incapr(c(1, 2), 0, 1)), quote(incapr(c(2, 1), 0, 0)),
    quote(incapr(c(2, 1), -1, 1)), quote(incapr(c(2, 1), 0.5, 1)),
    quote(incapr(c(2, 1), 1, 1)), quote(incapr(c(2, 1), 0, 1, r = 2)),
    quote(incapr(c(2, 1), 2147483647, 1)), quote(incapr(c(2, 1), 0, 3e9)),
    quote(sellab(c(1, 0, 1, 0), c(TRUE, TRUE, FALSE, TRUE),
                 c(TRUE, FALSE, TRUE, TRUE))),
    quote(sellab(c(1, 0), c(TRUE, NA), c(TRUE, TRUE))),
    quote(sellab(c(1, 0), c(TRUE, TRUE), c(1, 0))),
    quote(sellab(c(1, 2), c(TRUE, TRUE), c(TRUE, TRUE))),
    quote(sellab(c(NA, 1), c(FALSE, TRUE), c(TRUE, TRUE))),
    quote(sellab("a", TRUE, TRUE)),
    quote(sellab(c(1, 0), c(FALSE, FALSE), c(TRUE, TRUE))),
    quote(sellab(c(1, 0), c(TRUE, TRUE), c(FALSE, FALSE))),
    quote(losbnd(c(1, 2, 3), c(0.5, 4))), quote(losbnd(numeric(0), 1)),
    quote(sntbnd(c(0, 1, 1, 0), c("a", "b", "a", "b"))),
    quote(sntbnd(c(0, 1, 1, 0), c("a", "b", "a", "b"), treat = "a")),
    quote(sntbnd(c(0, 1, 1, 0), factor(c("x", "y", "x", "y")))),
    quote(sntbnd(c(0, 1, 1, 0), c(1, 2, 1, 2), treat = 1)),
    quote(sntbnd(c(0, 1, 1, 0), c(TRUE, FALSE, TRUE, FALSE))),
    quote(sntbnd(c(0, 2), c("a", "b"))), quote(sntbnd(c(0, 1), c("a", "a"))),
    quote(sntbnd(c(0, 1, 1), c("a", "b", "c"))),
    quote(sntbnd(c(0, 1), c("a", NA))), quote(sntbnd(c(0, 1), "a")),
    quote(sntbnd(c(0, 1), c("a", "b"), treat = "c")),
    quote(sntbnd(c(0, 1), c("a", "b"), treat = list("a"))),
    quote(sntbnd(c(0, 1), c("a", "b"), w = c(0, 0))),
    quote(hzdsel(0.3, 0.4, 0.1, c(0.9, 0.8), c(0.95, 0.9))),
    quote(hzdsel(1, 0.5, 0.1, c(1, 1), c(1, 1))),
    quote(hzdsel(0.5, 0.5, 0.6, c(1, 1), c(1, 1))),
    quote(hzdsel(0.5, 0.5, 0.1, 1, c(1, 1))),
    quote(hzdsel(0.5, 0.5, 0.1, c(1, 1), c(1, 1.2))),
    quote(proof()), quote(proof(c("cllr", "dbsrch"))),
    quote(proof(character(0))), quote(proof(c(a = "gtmp"))),
    quote(proof(NA_character_)), quote(proof(1)),
    quote(proof(c("cllr", "nope", "zzz", "nope")))
  )
  lst
}

run <- function(e) {
  w <- list()
  v <- withCallingHandlers(
    tryCatch(eval(e, globalenv()), error = function(err) {
      structure(list(cls = class(err), msg = conditionMessage(err),
                     call = conditionCall(err)), class = "cmperr")
    }),
    warning = function(x) {
      w[[length(w) + 1L]] <<- list(msg = conditionMessage(x),
                                   call = conditionCall(x))
      invokeRestart("muffleWarning")
    }
  )
  list(value = v, warn = w)
}

near <- function(x, y) {
  if (is.numeric(x) && is.numeric(y) && length(x) == length(y) &&
        identical(attributes(x), attributes(y)) &&
        identical(is.na(x), is.na(y))) {
    ok <- !is.na(x)
    return(all(x[ok] == y[ok] | abs(x[ok] - y[ok]) <=
                 1e-12 * pmax(abs(x[ok]), abs(y[ok]))))
  }
  if (is.list(x) && is.list(y) && length(x) == length(y) &&
        identical(attributes(x), attributes(y))) {
    return(all(mapply(near, x, y)))
  }
  identical(x, y)
}

if (a[1] == "run") {
  .libPaths(c(a[2], .libPaths()))
  library(crimlab)
  grDevices::pdf(NULL)
  cs <- cases()
  out <- lapply(cs, run)
  saveRDS(list(cases = cs, out = out), a[3])
} else {
  o <- readRDS(a[2])
  n <- readRDS(a[3])
  k <- c(same = 0, near = 0, diff = 0)
  for (i in seq_along(o$out)) {
    x <- o$out[[i]]
    y <- n$out[[i]]
    s <- if (identical(x, y, attrib.as.set = FALSE)) {
      "same"
    } else if (near(x, y)) {
      "near"
    } else {
      "diff"
    }
    k[s] <- k[s] + 1
    if (s != "same") {
      cat(toupper(s), deparse(o$cases[[i]], width.cutoff = 500L)[1], "\n")
      if (s == "diff") {
        utils::str(x, give.attr = TRUE, vec.len = 3, nchar.max = 300)
        utils::str(y, give.attr = TRUE, vec.len = 3, nchar.max = 300)
      }
    }
  }
  cat(sprintf("CMP identical %d, within 1e-12 %d, different %d\n", k[["same"]],
              k[["near"]], k[["diff"]]))
}
