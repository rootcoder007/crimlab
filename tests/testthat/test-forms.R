# Arguments that hold one value per group (per unit, per locus, per arm):
# the forms they accept and the rule that matches them by name (src/bind.cpp,
# "Group arguments"). Names are the identity of a group: reordering the
# names of one argument never changes a result, and names that differ are
# an error. Every expectation compares two calls of the package; none locks
# an observed number.

# r with its rows in the order of the groups g (column 1), row names reset
byg <- function(r, g) {
  r <- r[match(g, r[[1L]]), , drop = FALSE]
  rownames(r) <- NULL
  r
}

# The rule on f, for the named group arguments a (all named alike, in one
# order), the other arguments x and out() to reach the per-group table.
rule <- function(f, a, x = list(), out = identity) {
  call <- function(a) out(do.call(f, c(a, x)))
  g <- names(a[[1L]])
  o <- c(seq_along(g)[-1L], 1L)
  r <- call(a)
  # every argument permuted alike, names kept: the same rows, reordered
  expect_equal(call(lapply(a, `[`, o)), byg(r, g[o]), tolerance = 1e-14,
               info = f)
  for (k in seq_along(a)[-1L]) {
    b <- a
    b[[k]] <- a[[k]][o]
    # one argument permuted with its names: unchanged
    expect_identical(call(b), r, info = f)
    # one argument permuted without names: positions, the values as given
    b[[k]] <- unname(a[[k]][o])
    p <- a
    p[[k]] <- stats::setNames(unname(a[[k]][o]), g)
    expect_identical(call(b), call(p), info = f)
    # names that differ are an error that lists them
    b[[k]] <- stats::setNames(a[[k]], c("zz", g[-1L]))
    expect_error(call(b), paste0("missing: ", g[1L], "; extra: zz"),
                 fixed = TRUE, info = f)
  }
  # a name given twice is an error
  b <- a
  names(b[[1L]])[2L] <- g[1L]
  expect_error(call(b), "names a group twice", info = f)
  # the other forms: a list, a one-row data frame, a 1-d table, integers
  expect_identical(call(lapply(a, as.list)), r, info = f)
  expect_identical(call(lapply(a, function(v) as.data.frame(as.list(v)))), r,
                   info = f)
  expect_identical(call(lapply(a, as.table)), r, info = f)
  w <- vapply(a, function(v) all(v == round(v)), NA)
  b <- a
  b[w] <- lapply(a[w], function(v) stats::setNames(as.integer(v), names(v)))
  expect_identical(call(b), r, info = f)
  # unnamed: positions, groups "1", "2", ... (ref by name becomes a position)
  if (is.character(x[["ref"]])) x[["ref"]] <- match(x[["ref"]], g)
  u <- call(lapply(a, unname))
  expect_identical(u[-1L], r[-1L], info = f)
  expect_identical(u[[1L]][seq_along(g)], as.character(seq_along(g)),
                   info = f)
  # some named and some not: positions, named by the named one
  b <- a
  b[[1L]] <- unname(a[[1L]])
  expect_identical(call(b), r, info = f)
  r
}

test_that("vodtst matches day, night and shift by name", {
  a <- list(day = c(white = 600, black = 400, asian = 100),
            night = c(white = 700, black = 300, asian = 150),
            shift = c(white = 1, black = 1.3, asian = 0.9))
  r <- rule("vodtst", a, list(ref = "black"))
  # the reported misread: night named in another order
  expect_identical(
    vodtst(c(white = 600, black = 400), c(black = 300, white = 700)),
    vodtst(c(white = 600, black = 400), c(white = 700, black = 300))
  )
  expect_identical(vodtst(list(white = 600, black = 400),
                          list(white = 700, black = 300)),
                   vodtst(c(white = 600, black = 400),
                          c(white = 700, black = 300)))
  # a one-row or one-column matrix gives its names
  expect_identical(vodtst(t(a$day), cbind(a$night)),
                   vodtst(a$day, a$night))
  expect_identical(vodtst(matrix(c(1, 2), 1), c(1, 1)),
                   vodtst(c(1, 2), c(1, 1)))
  expect_identical(r[["ref"]], c(FALSE, TRUE, FALSE))
})

test_that("ref is found by name after matching, or by position as given", {
  d <- c(a = 50, b = 30, c = 20)
  n <- c(a = 40, b = 40, c = 20)
  o <- c(3L, 1L, 2L)
  for (f in c("vodtst", "dspbm", "expbnd")) {
    x <- switch(f, vodtst = list(d, n), dspbm = list(d, n, n),
                expbnd = list(d, n, 1.5))
    r <- do.call(f, c(x, ref = "b"))
    p <- do.call(f, c(lapply(x, function(v) if (length(v) > 1L) v[o] else v),
                      ref = "b"))
    expect_equal(p, byg(r, names(d)[o]), tolerance = 1e-14, info = f)
    # a number is a position in the first argument as given
    q <- do.call(f, c(lapply(x, function(v) if (length(v) > 1L) v[o] else v),
                      ref = 3))
    expect_identical(q[["group"]][q[["ref"]]], "b", info = f)
    expect_equal(q, p, tolerance = 1e-14, info = f)
    # a factor is read by its label, never its code
    expect_identical(do.call(f, c(x, ref = factor("b", levels = c("c", "b")))),
                     r, info = f)
    expect_identical(r[["group"]][r[["ref"]]], "b", info = f)
  }
  expect_error(vodtst(d, n, ref = factor(NA)), "`ref`")
})

test_that("groups from a factor match groups from character by name", {
  g <- c("w", "b", "w", "a", "b", "w")
  lv <- c("w", "b", "a")
  d1 <- table(factor(g, levels = lv))
  d2 <- table(g)
  n <- c(w = 5, b = 3, a = 2)
  expect_equal(vodtst(d1, n, ref = "w"), byg(vodtst(d2, n, ref = "w"), lv),
               tolerance = 1e-14)
  expect_identical(vodtst(d2, n[names(d2)]), vodtst(d2, n))
})

test_that("dspbm, expbnd and clrsmp match their groups by name", {
  a <- list(pop = c(a = 9000, b = 1000, c = 500),
            contact = c(a = 900, b = 300, c = 40),
            force = c(a = 9, b = 6, c = 2))
  rule("dspbm", a, list(ref = "a"))
  rule("expbnd", list(y = c(a = 100, b = 300, c = 7),
                      m = c(a = 1000, b = 1000, c = 70)),
       list(gamma = 1.5, ref = "a"))
  rule("clrsmp", list(clra = c(minor = 90, serious = 20, other = 5),
                      reca = c(minor = 100, serious = 50, other = 10),
                      clrb = c(minor = 19, serious = 60, other = 4),
                      recb = c(minor = 20, serious = 140, other = 9)),
       out = function(r) r[["parts"]])
})

test_that("incap matches lam, q, s and share by name", {
  a <- list(lam = c(low = 2, mid = 4, high = 10),
            q = c(low = 0.1, mid = 0.2, high = 0.3),
            s = c(low = 1, mid = 2, high = 3))
  rule("incap", a)
  # the aggregate row of the shares, reordered by name
  r <- incap(c(low = 2, high = 10), 0.1, 1, share = c(low = 0.8, high = 0.2))
  expect_identical(
    incap(c(low = 2, high = 10), 0.1, 1, share = c(high = 0.2, low = 0.8)), r
  )
  expect_equal(r[["lam"]][3], 0.8 * 2 + 0.2 * 10, tolerance = 1e-15)
})

test_that("gtmp matches q and hom to p by locus", {
  p <- c(D3 = 0.1, vWA = 0.2, FGA = 0.05)
  q <- c(D3 = 0.3, vWA = 0.4, FGA = 0.2)
  rule("gtmp", list(p = p, q = q), list(hom = FALSE),
       out = function(r) r[["locus"]])
  h <- c(D3 = FALSE, vWA = TRUE, FGA = FALSE)
  o <- c(2L, 3L, 1L)
  expect_identical(gtmp(p, q, h[o]), gtmp(p, q, h))
  expect_identical(gtmp(p, q, h, 0.02)[["locus"]][["hom"]], unname(h))
  expect_error(gtmp(p, q, c(D3 = FALSE, vWA = TRUE, X = FALSE)), "`hom`")
  # q as a list of numbers or NA, matched by name
  expect_identical(gtmp(p, as.list(q[o])), gtmp(p, q))
})

test_that("hzdsel matches its arms by name", {
  r <- hzdsel(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.7))
  sh <- c(x = 0.5, y = 0.8)
  sl <- c(x = 0.9, y = 0.7)
  expect_identical(hzdsel(0.5, 0.5, 0.1, sh, sl[2:1]), r)
  expect_identical(hzdsel(0.5, 0.5, 0.1, as.list(sh),
                          as.data.frame(as.list(sl))), r)
  expect_error(hzdsel(0.5, 0.5, 0.1, sh, c(x = 0.9, z = 0.7)),
               "missing: y; extra: z", fixed = TRUE)
  xx <- stats::setNames(c(0.5, 0.8), c("x", "x"))
  expect_error(hzdsel(0.5, 0.5, 0.1, xx, sl),
               "`surh` names a group twice")
  # arms named control and treated are read by name in either order
  expect_identical(hzdsel(0.5, 0.5, 0.1, c(treated = 0.8, control = 0.5),
                          c(treated = 0.7, control = 0.9)), r)
  # other names keep the documented order, control first
  expect_identical(hzdsel(0.5, 0.5, 0.1, c(x = 0.5, y = 0.8),
                          c(y = 0.7, x = 0.9)), r)
  expect_identical(hzdsel(0.5, 0.5, 0.1, c(treated = 0.5, control = 0.8),
                          c(0.9, 0.7)), r)
})

test_that("dyadrt matches count rows, columns and pop by name", {
  pop <- c(a = 800, b = 150, c = 50)
  m <- matrix(c(9, 2, 1, 3, 8, 1, 2, 2, 6), 3, dimnames = list(names(pop),
                                                               names(pop)))
  r <- dyadrt(m, pop)
  key <- function(x) x[order(x$off, x$vic), ]
  o <- c(3L, 1L, 2L)
  expect_identical(dyadrt(m[o, o], pop), r)
  expect_identical(dyadrt(m[o, ], pop), r)
  expect_equal(key(dyadrt(m, pop[o])), key(r), tolerance = 1e-14,
               ignore_attr = TRUE)
  # without pop names, the columns follow the rows by name
  expect_identical(dyadrt(m[, o], unname(pop)), dyadrt(m, unname(pop)))
  expect_identical(dyadrt(m, unname(pop))[["off"]], r[["off"]])
  # names only on the columns name the groups; none: "1", "2", ...
  expect_identical(dyadrt(unname(m), pop)[-(1:2)], r[-(1:2)])
  cm <- m
  rownames(cm) <- NULL
  expect_identical(dyadrt(cm, unname(pop))[1:2], r[1:2])
  expect_identical(dyadrt(unname(m), unname(pop))[["off"]][1:3],
                   c("1", "2", "3"))
  # a 2-d table from factors
  off <- factor(c("a", "b", "a", "c"), levels = c("c", "b", "a"))
  vic <- factor(c("b", "a", "a", "c"), levels = c("c", "b", "a"))
  t2 <- table(off, vic)
  expect_identical(dyadrt(t2, pop),
                   dyadrt(unclass(t2)[names(pop), names(pop)] + 0, pop))
  expect_error(dyadrt(m, c(a = 1, b = 1, z = 1)), "missing: z; extra: c",
               fixed = TRUE)
  expect_error(dyadrt(m, c(x = 1, b = 1, z = 1)),
               "missing: x, z; extra: a, c", fixed = TRUE)
  mm <- m
  colnames(mm) <- c("a", "b", "z")
  expect_error(dyadrt(mm, unname(pop)), "the column names of `count`")
  rownames(mm) <- c("a", "a", "b")
  expect_error(dyadrt(mm, pop), "the row names of `count` names a group twice")
  expect_error(dyadrt(m, list(a = 1, b = 1, c = "x")), "`pop`")
})

test_that("per-unit arguments are matched by unique names", {
  q <- c(a = 0.1, b = 0.5, c = 0.9)
  w <- c(a = 3, b = 1, c = 2)
  o <- c(3L, 1L, 2L)
  expect_identical(outtst(q, 0.3, w[o]), outtst(q, 0.3, w))
  expect_equal(outtst(q[o], 0.3, w[o]), outtst(q, 0.3, w), tolerance = 1e-14)
  expect_identical(outtst(q, 0.3, unname(w[o])),
                   outtst(unname(q), 0.3, unname(w[o])))
  expect_error(outtst(q, 0.3, c(a = 1, b = 1, z = 1)), "`w` must name the")
  # repeated names (labels, not identities) keep positions
  r <- stats::setNames(c(0.1, 0.5, 0.9), c("x", "x", "y"))
  expect_identical(outtst(r, 0.3, w[o]), outtst(unname(r), 0.3, unname(w[o])))
  x1 <- c(a = 1, b = 9, c = 2)
  x2 <- c(a = 2, b = 3, c = 1)
  expect_identical(rtmhs(x1, x2[o], 5), rtmhs(x1, x2, 5))
  expect_identical(rtmhs(x1, x2, 5, w[o]), rtmhs(x1, x2, 5, w))
  lr <- c(a = 10, b = 1e4)
  pr <- c(a = 0.01, b = 0.5)
  expect_identical(lrpost(lr, pr[2:1]), lrpost(lr, pr))
  y <- c(a = 1, b = 0, c = 1)
  rel <- c(a = TRUE, b = TRUE, c = FALSE)
  ru <- c(a = TRUE, b = TRUE, c = TRUE)
  expect_identical(sellab(y, rel[o], ru[o], w[o]), sellab(y, rel, ru, w))
  z <- c(a = "x", b = "y", c = "y")
  expect_identical(sntbnd(y, z[o], "y"), sntbnd(y, z, "y"))
  expect_identical(sntbnd(y, as.list(z)[o], "y"), sntbnd(y, z, "y"))
  zf <- factor(z, levels = c("y", "x"))
  expect_identical(sntbnd(y, zf[o], "y", w[o]), sntbnd(y, z, "y", w))
  # factor and character z give the same answer
  expect_identical(sntbnd(y, zf, "y"), sntbnd(y, z, "y"))
  expect_identical(sntbnd(y, zf, factor("y")), sntbnd(y, z, "y"))
})

test_that("sntbnd takes the treatment by name, never by sorting", {
  y <- c(1, 0, 1, 0)
  expect_error(sntbnd(y, c("x", "y", "x", "y")), "`treat` must name the")
  expect_identical(sntbnd(y, c(TRUE, FALSE, TRUE, FALSE))[["term"]],
                   c("FALSE", "TRUE", "effect"))
  expect_identical(sntbnd(y, c(1, 0, 1, 0))[["term"]], c("0", "1", "effect"))
  expect_identical(sntbnd(y, c(0L, 1L, 0L, 1L)),
                   sntbnd(y, c(0L, 1L, 0L, 1L), 1))
  expect_error(sntbnd(y, c(1, 2, 1, 2)), "`treat`")
  expect_error(sntbnd(y, factor(c(0, 1, 0, 1))), "`treat`")
})

test_that("genotypes are matched by locus name", {
  g1 <- rbind(L1 = c(12, 14), L2 = c(8, 9), L3 = c(5, 5))
  g2 <- rbind(L1 = c(12, 15), L2 = c(8, 9), L3 = c(5, 6))
  f <- rbind(L1 = c(0.1, 0.2), L2 = c(0.3, 0.3), L3 = c(0.25, 0.25))
  o <- c(3L, 1L, 2L)
  r <- kinidx(g1, g2, f, "sibling")
  expect_identical(kinidx(g1, g2[o, ], f, "sibling"), r)
  expect_identical(kinidx(g1, g2, f[o, ], "sibling"), r)
  p <- kinidx(g1[o, ], g2[o, ], f[o, ], "sibling")
  expect_equal(p$locus, byg(r$locus, rownames(g1)[o]), tolerance = 1e-14)
  # character alleles move with their rows
  ch <- function(g) {
    storage.mode(g) <- "character"
    g
  }
  expect_identical(kinidx(ch(g1), ch(g2)[o, ], f, "sibling"), r)
  # a table without row names: positions, loci named by the named one
  expect_identical(kinidx(unname(g1), g2, unname(f), "sibling")$locus$locus,
                   rownames(g2))
  g3 <- g2
  rownames(g3)[3] <- "L9"
  expect_error(kinidx(g1, g3, f), "missing: L3; extra: L9", fixed = TRUE)
  rownames(g3)[3] <- "L1"
  expect_error(kinidx(g1, g3, f), "`g2` names a group twice")
  m <- rbind(L1 = c(12, 12), L2 = c(8, 8), L3 = c(5, 6))
  r <- patidx(g1, g2, f, m)
  expect_identical(patidx(g1, g2[o, ], f[o, ], m[o, ]), r)
  # two alleles for one locus keep working
  expect_identical(kinidx(c(12, 14), c(12, 14), c(0.1, 0.2))$locus$locus, "1")
})

test_that("vodtst takes one table of day and night columns, by name", {
  d <- c(white = 600, black = 400)
  n <- c(white = 700, black = 300)
  r <- vodtst(d, n)
  df <- data.frame(day = d, night = n, row.names = names(d))
  expect_identical(vodtst(df), r)
  expect_identical(vodtst(df[, c("night", "day")]), r)
  expect_identical(vodtst(cbind(night = n, day = d)), r)
  expect_identical(vodtst(cbind(day = d, night = n), shift = c(1, 2)),
                   vodtst(d, n, shift = c(white = 1, black = 2)))
  # a 2-d table from xtabs(), light in the columns
  x <- data.frame(g = c("white", "black", "white", "black"),
                  light = c("night", "night", "day", "day"),
                  k = c(700, 300, 600, 400))
  expect_identical(vodtst(stats::xtabs(k ~ g + light, x), ref = "white"),
                   byg(vodtst(d, n), c("black", "white")))
  # automatic row names: groups "1", "2"
  expect_identical(vodtst(data.frame(day = unname(d), night = unname(n))),
                   vodtst(unname(d), unname(n)))
  expect_error(vodtst(cbind(d, n)), "two columns named day and night")
  expect_error(vodtst(cbind(dark = n, day = d)), "named day and night")
  expect_error(vodtst(d), "`night` is missing")
  expect_error(vodtst(matrix(c("1", "2", "3", "4"), 2,
                             dimnames = list(NULL, c("day", "night")))),
               "`day` must be finite numbers, one per group")
  expect_error(vodtst(data.frame(day = c("a", "b"), night = c(1, 2))),
               "`day` must be finite numbers, one per group")
})

test_that("a group argument in another form is refused, listing the forms", {
  msg <- "a numeric vector, a 1-d table, a list of single numbers or a one-row"
  expect_error(vodtst(list(a = 1, b = c(1, 2)), c(1, 1)), msg)
  expect_error(vodtst(data.frame(a = c(1, 2), b = c(1, 2)), c(1, 1)), msg)
  expect_error(vodtst(c("1", "2"), c(1, 1)), msg)
  expect_error(dspbm(list(a = "1"), 1, 1), msg)
  # NA in a list is a missing value, refused as one
  expect_error(vodtst(list(a = 1, b = NA), c(1, 1)), "no missing values")
})

test_that("labelled vectors are refused: codes are not values", {
  lab <- function(x) {
    structure(x, labels = c(one = 1, two = 2),
              class = c("haven_labelled", "vctrs_vctr", class(x)))
  }
  m <- "has value labels"
  expect_error(vodtst(lab(c(1, 2)), c(1, 1)), m)
  expect_error(vodtst(list(a = lab(1), b = 2), c(1, 1)), m)
  expect_error(vodtst(structure(cbind(c(1, 2)), labels = c(one = 1)),
                      c(1, 1)), m)
  expect_error(vodtst(structure(cbind(day = c(1, 2), night = c(1, 1)),
                                labels = c(one = 1))), m)
  expect_error(vodtst(c(1, 2), c(1, 1), ref = lab(1)), m)
  expect_error(orrisk(lab(2), 0.5, 0.5), m)
  expect_error(orrisk(2, lab(0.5), 0.5), m)
  expect_error(misrd(lab(1), 0, 0.5, 1, 1), m)
  expect_error(incapr(c(3, 2, 1), lab(0), 1), m)
  expect_error(sellab(lab(c(1, 0)), c(TRUE, TRUE), c(TRUE, TRUE)), m)
  expect_error(sellab(c(1, 0), structure(c(TRUE, TRUE), labels = c(a = TRUE)),
                      c(TRUE, TRUE)), m)
  expect_error(sntbnd(c(1, 0), lab(c(1, 2)), 1), m)
  expect_error(sntbnd(c(1, 0), c("a", "b"), lab(1)), m)
  expect_error(gtmp(0.1, lab(0.2)), m)
  expect_error(kinidx(lab(c(1, 2)), c(1, 2), c(0.1, 0.2)), m)
  expect_error(dyadrt(structure(matrix(1, 2, 2), labels = c(one = 1)),
                      c(1, 1)), m)
})
