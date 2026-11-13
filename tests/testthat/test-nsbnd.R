# nsbnd(), impfpr() and rnkres(): against the inversion of the recording
# model searched over a grid of the noise boxes, the confusion tables the
# identity is about, brute force over pairs, the earlier package's recorded
# values, and refusals.

ref <- readRDS(test_path("ref-fair.rds"))
# the true rates a recorded one allows over a fine grid of the two boxes
# (with the points a = pobs and b = 1 - pobs where they fall inside): the
# noise rates under which some true rate in [0, 1] gives pobs
grid <- function(q, a, b, k = 201) {
  if (length(a) == 1L) a <- c(0, a)
  if (length(b) == 1L) b <- c(0, b)
  sa <- c(seq(a[1], a[2], length.out = k), q[q >= a[1] & q <= a[2]])
  r <- 1 - q
  sb <- c(seq(b[1], b[2], length.out = k), r[r >= b[1] & r <= b[2]])
  g <- expand.grid(a = sa, b = sb)
  p <- (q - g$a) / (1 - g$a - g$b)
  p <- p[p >= -1e-12 & p <= 1 + 1e-12]
  if (!length(p)) return(c(NaN, NaN))
  range(p)
}

test_that("the bounds are the ends of the identified set", {
  q <- c(0.05, 0.2, 0.35, 0.55, 0.9, 0.97)
  for (bx in list(list(c(0, 0.1), c(0, 0.2)), list(c(0.02, 0.1), c(0.05, 0.3)),
                  list(c(0, 0), c(0, 0.4)), list(c(0.1, 0.25), c(0, 0)))) {
    r <- nsbnd(q, bx[[1]], bx[[2]])$bnd
    for (i in seq_along(q)) {
      e <- grid(q[i], bx[[1]], bx[[2]])
      # the grid holds both ends exactly, as its corners: the set is the
      # interval between them
      expect_equal(c(r$lo[i], r$hi[i]), e, tolerance = 1e-12)
    }
    expect_equal(r$width, r$hi - r$lo, tolerance = 0)
  }
  # the true rate stays within the bounds for any admissible noise
  set.seed(1)
  for (k in 1:200) {
    p <- runif(1)
    a <- runif(1, 0, 0.2)
    b <- runif(1, 0, 0.3)
    r <- nsbnd(p * (1 - b) + (1 - p) * a, 0.2, 0.3)$bnd
    expect_lte(r$lo, p + 1e-12)
    expect_gte(r$hi + 1e-12, p)
  }
})

test_that("known noise gives the true rate; no admissible noise, none", {
  r <- nsbnd(c(0.38, 0.5), c(0.1, 0.1), c(0.2, 0.2))$bnd
  expect_equal(r$lo, ref$tr, tolerance = 1e-15)
  expect_equal(r$hi, ref$tr, tolerance = 1e-15)
  expect_identical(r$width, c(0, 0))
  # below the smallest false-positive rate, or above 1 minus the smallest
  # false-negative one, no true rate is recorded as pobs
  r <- nsbnd(c(0.01, 0.5, 0.99), c(0.05, 0.1), c(0.05, 0.1))
  expect_identical(is.nan(r$bnd$lo), c(TRUE, FALSE, TRUE))
  expect_identical(r$pairs$ord, c(NA, NA, NA_character_))
  # at the ends of [0, 1]
  r <- nsbnd(c(0, 1), 0.1, 0.2)$bnd
  expect_identical(c(r$lo, r$hi), c(0, 1, 0, 1))
})

test_that("orders and breakdowns agree with the earlier package's", {
  r <- nsbnd(c(0.35, 0.55, 0.2), 0.1, 0.2)$bnd
  expect_equal(r$lo, ref$bnd$lower, tolerance = 1e-15)
  expect_equal(r$hi, ref$bnd$upper, tolerance = 1e-15)
  ord <- c("a < b" = "lower", "b < a" = "higher", undecided = "undecided")
  args <- list(c(0.35, 0.55, 0.05, 0.20), c(0.35, 0.55, 0.05, 0.40),
               c(0.6, 0.2, 0.1, 0.1), c(0.30, 0.31, 0, 0.05))
  for (k in seq_along(args)) {
    v <- args[[k]]
    p <- nsbnd(v[1:2], v[3], v[4])$pairs
    expect_identical(p$ord, unname(ord[ref$cmp[[k]]$order]))
    expect_equal(p$bd, ref$cmp[[k]]$breakdown_beta_max, tolerance = 1e-15)
  }
})

test_that("an order is identified exactly when it survives every noise", {
  set.seed(2)
  for (k in 1:60) {
    q <- sort(runif(2, 0.05, 0.9))
    a <- runif(1, 0, 0.15)
    b <- runif(1, 0, 0.4)
    p <- nsbnd(q, a, b)$pairs
    # brute force: the identified sets of the two true rates
    s1 <- grid(q[1], a, b, 81)
    s2 <- grid(q[2], a, b, 81)
    want <- if (s1[2] < s2[1]) "lower" else "undecided"
    expect_identical(p$ord, want)
    # the breakdown: decided just below it, not just above
    if (!is.na(p$bd) && p$bd < 1 - a - 1e-6) {
      expect_identical(nsbnd(q, a, p$bd * (1 - 1e-6))$pairs$ord, "lower")
      up <- min(p$bd * (1 + 1e-6), 1 - a - 1e-9)
      expect_identical(nsbnd(q, a, up)$pairs$ord, "undecided")
    }
  }
  # pairs come in the order 1-2, 1-3, 2-3, with the rates' names
  p <- nsbnd(c(x = 0.1, y = 0.5, z = 0.3), 0.05, 0.05)$pairs
  expect_identical(paste(p$i, p$j), c("x y", "x z", "y z"))
  expect_identical(p$ord, c("lower", "lower", "higher"))
  # equal rates: undecided, no breakdown
  p <- nsbnd(c(0.3, 0.3), 0, 0)$pairs
  expect_identical(p$ord, "undecided")
  expect_true(is.na(p$bd))
  # the cap: the earlier package stops 1e-12 short of 1 - a; nsbnd()
  # reports the end itself, the supremum (no box at it is admissible)
  expect_identical(nsbnd(c(0.01, 0.9), 0.1, 0)$pairs$bd, 0.9)
  expect_identical(nrow(nsbnd(0.3, 0.1, 0.1)$pairs), 0L)
})

test_that("impfpr() is the false positive rate of every confusion table", {
  set.seed(3)
  for (k in 1:100) {
    tab <- runif(4, 1, 100)
    tp <- tab[1]
    fp <- tab[2]
    fn <- tab[3]
    tn <- tab[4]
    p <- (tp + fn) / sum(tab)
    v <- impfpr(p, tp / (tp + fp), fn / (tp + fn))
    expect_equal(v, fp / (fp + tn), tolerance = 1e-13)
  }
  expect_equal(impfpr(0.3, 0.6, 0.25), ref$fpr, tolerance = 1e-15)
  # equal ppv and fnr, unequal base rates: the false positive rates differ
  v <- impfpr(c(a = 0.2, b = 0.4), 0.7, 0.3)
  expect_named(v, c("a", "b"))
  expect_gt(v[["b"]], v[["a"]])
  expect_identical(impfpr(0.5, 1, 0), 0)
  # named inputs are matched by name
  f <- c(a = 0.1, b = 0.2)
  expect_identical(impfpr(c(a = 0.2, b = 0.4), c(b = 0.5, a = 0.7), f),
                   impfpr(c(a = 0.2, b = 0.4), c(a = 0.7, b = 0.5), f))
})

test_that("rnkres() identifies the pairs whose intervals do not meet", {
  e <- c(0.2, 0.35, 0.8, 0.5)
  w <- c(0.1, 0.1, 0.05, 0.2)
  r <- rnkres(e, w)
  id <- abs(outer(e, e, "-")) > outer(w, w, "+")
  diag(id) <- NA
  expect_identical(unname(r$id), id)
  expect_identical(unname(r$id), unname(ref$rank$identified_pairs))
  expect_identical(r$npair, 6L)
  expect_identical(r$nunid, sum(!id[upper.tri(id)]))
  expect_equal(r$share, ref$rank$share_unidentified, tolerance = 1e-15)
  expect_identical(r$res, ref$rank$resolution)
  expect_identical(dimnames(r$id), list(c("1", "2", "3", "4"),
                                        c("1", "2", "3", "4")))
  r <- rnkres(c(a = 0, b = 1), 0.5)
  expect_false(r$id["a", "b"])
  expect_identical(r$share, 1)
  expect_true(rnkres(c(0, 1), 0)$id[1, 2])
})

test_that("nsbnd(), impfpr() and rnkres() refuse what they cannot use", {
  errs("nsbnd", list(
    list("rates from 0 to 1", 1.2, 0.1, 0.1),
    list("must be at least 0", -0.1, 0.1, 0.1),
    list("one number", 0.3, c(0, 0.1, 0.2), 0.1),
    list("lower end first", 0.3, c(0.2, 0.1), 0.1),
    list("add to less than 1", 0.3, 0.5, 0.5)
  ))
  errs("impfpr", list(
    list("base rates in", 1, 0.5, 0.5),
    list("`ppv` must be in", 0.5, 0, 0.5),
    list("`fnr` must be in", 0.5, 0.5, 1),
    list("one length", c(0.1, 0.2), c(0.5, 0.5, 0.5), 0.1),
    list("finite numbers", NA, 0.5, 0.5)
  ))
  errs("rnkres", list(
    list("at least two", 1, 0.1),
    list("one half-width", 1:3, c(0.1, 0.2)),
    list("must be at least 0", 1:3, -1)
  ))
})
