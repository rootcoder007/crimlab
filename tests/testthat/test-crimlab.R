# The kernels' own tests live in C++ (src/t_*.cpp): hand formulas,
# brute-force enumeration and seeded simulation, checked through crimlab.h.
# This runs them all; a family that stops running drops the count below the
# floor and fails here. The floor sits a little below the count on x86-64
# (13787): a few checks are conditional on the random cases, and the
# smallest family has more than 300.

test_that("the compiled self-test passes every check", {
  r <- .Call(C_clab_selftest)
  expect_identical(r[["report"]], "")
  expect_identical(r[["nfail"]], 0L)
  expect_gte(r[["n"]], 13700L)
})

test_that("every exported function answers one valid call in its shape", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  calls <- list(
    list(quote(vodtst(c(w = 6, b = 4), c(5, 5), shift = c(1, 2))),
                  c("group", "ref", "shday", "shnight", "ratio",
                    "sign", "adj")),
    list(quote(orrisk(3, 0.4, 0.5)),
                  c("exposed", "unexposed", "rr", "or")),
    list(quote(outtst(c(0.1, 0.5, 0.9), 0.3)),
                  c("stop", "hit")),
    list(quote(dspbm(c(9, 1), c(3, 1), c(2, 1), ref = 2)),
                 c("group", "ref", "cd", "fd", "rd")),
    list(quote(dyadrt(matrix(1, 2, 2), c(a = 1, b = 2))),
                  c("off", "vic", "count", "rtoff", "rtpair",
                    "nulloff", "ratio")),
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
                  c("group", "ref", "rate", "lo", "hi", "ratio", "rlo",
                    "rhi", "sign")),
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
    list(quote(coofg(c(1, 1), c("a", "b"))), c("edges", "nodes")),
    list(quote(ntdeg("a", "b")), c("node", "deg", "wdeg")),
    list(quote(ntbtw("a", "b")), c("node", "btw")),
    list(quote(ntcls("a", "b")), c("node", "reach", "cls", "ncls", "hrm")),
    list(quote(ntevc("a", "b")), c("cent", "value", "iter", "conv")),
    list(quote(ntcore("a", "b")), c("node", "core")),
    list(quote(ntcmp("a", "b")), c("node", "weak", "strong")),
    list(quote(ntlvn("a", "b")), c("memb", "q", "levels")),
    list(quote(ntmod("a", "b", c(1, 1))), NULL),
    list(quote(ntbrk("a", "b", c(1, 2))),
                 c("node", "grp", "coord", "itin", "repr", "gate", "liais",
                   "total")),
    list(quote(ntchg("a", "b")), c("bound", "set", "iter", "conv")),
    list(quote(proof("cllr")), c("fn", "theorem", "file")),
    list(quote(relchk(matrix(1:4, 2), 2)),
         c("verdict", "ok", "value", "tol", "digits", "type", "claim",
           "reads")),
    list(quote(catchk(c("a", "b"))), c("ok", "nmin", "vars", "levels")),
    list(quote(cntreg(c(3, 5, 9, 4), c(1, 1, 2, 2))),
         c("irr", "coef", "fit", "ref", "n", "ncl", "level")),
    list(quote(grprr(c(3, 5), c(10, 20))),
         c("rates", "cmp", "ref", "type", "per", "level")),
    list(quote(dastd(c(3, 5), c(10, 20), c("a", "b"), std = c(1, 1))),
         c("rates", "cmp", "ref", "per", "level")),
    list(quote(sirex(c(3, 5), c(10, 20), c("a", "a"), c("x", "y"))),
         c("sir", "prior", "level")),
    list(quote(ctabx(matrix(1:4, 2))),
         c("tab", "expd", "res", "stat", "drop", "agg", "correct", "level")),
    list(quote(clprep()), c("preset", "generic", "args", "source")),
    list(quote(enpath(cbind(1:4, c(2, 1, 4, 3)), c(1, 3, 2, 5), nlam = 3)),
         c("lambda", "a0", "beta", "cvm", "cvsd", "lmin", "l1se")),
    list(quote(lgtreg(cbind(c(1, 3, 2, 5, 4, 6)), c(0, 0, 1, 0, 1, 1))),
         c("coef", "se", "dev")),
    list(quote(rfor(cbind(1:6), c(1, 3, 2, 5, 4, 6), ntree = 2)),
         c("pred", "oob")),
    list(quote(dmlx(c(1, 3, 2, 5, 4, 6), c(2, 1, 4, 3, 6, 5),
                    cbind(c(1, 2, 2, 3, 5, 4)), "plr", rep(1:2, 3),
                    lrn = "glm")),
         c("coef", "se", "reps", "psia", "psib", "pred"))
  )
  fns <- vapply(calls, function(x) as.character(x[[1L]][[1L]]), "")
  expect_setequal(unique(fns), setdiff(getNamespaceExports("crimlab"), spx))
  for (k in seq_along(calls)) {
    expect_named(eval(calls[[k]][[1L]]), calls[[k]][[2L]], info = fns[k])
  }
})
