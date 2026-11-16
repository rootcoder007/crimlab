test_that("dagcol() lists the colliders, shielded or not", {
  r <- dagcol(c("a", "b", "a", "d"), c("c", "c", "b", "c"))
  expect_s3_class(r, "dagcol")
  expect_identical(r$nodes, c("a", "c", "b", "d"))
  expect_identical(r$nedge, 4L)
  expect_identical(r$col$a, c("a", "a", "b"))
  expect_identical(r$col$c, c("c", "c", "c"))
  expect_identical(r$col$b, c("b", "d", "d"))
  expect_identical(r$col$shielded, c(TRUE, FALSE, FALSE))
  expect_null(r$query)
  expect_output(print(r), "dagcol: 3 colliders among 4 nodes, 2 unshielded")
  q <- dagcol(c("a", "b"), c("c", "c"), triple = c("b", "c", "a"))
  expect_identical(q$query, c(collider = TRUE, shielded = FALSE))
  expect_output(print(q), "triple: collider TRUE, shielded FALSE")
  ch <- dagcol(c("a", "c"), c("c", "b"), triple = c("a", "c", "b"))
  expect_identical(nrow(ch$col), 0L)
  expect_identical(ch$query, c(collider = FALSE, shielded = FALSE))
  expect_identical(dagcol(c("a", "a", "b"), c("c", "c", "c"))$nedge, 2L)
  expect_identical(dagcol(c(1, 2), c(3, 3))$col$c, "3")
})

test_that("dagcol() agrees with the colliders found by brute force", {
  set.seed(6)
  am <- matrix(stats::rbinom(64, 1, 0.35), 8) * upper.tri(diag(8))
  w <- which(am == 1, arr.ind = TRUE)
  g <- dagcol(LETTERS[w[, 1]], LETTERS[w[, 2]])
  nd <- g$nodes
  id <- match(nd, LETTERS)
  a <- character(0)
  cc <- character(0)
  b <- character(0)
  sh <- logical(0)
  for (k in seq_along(nd)) {
    ps <- nd[am[id, id[k]] == 1]
    if (length(ps) < 2L) next
    pr <- utils::combn(length(ps), 2)
    a <- c(a, ps[pr[1, ]])
    b <- c(b, ps[pr[2, ]])
    cc <- c(cc, rep(nd[k], ncol(pr)))
    i1 <- match(ps[pr[1, ]], LETTERS)
    j1 <- match(ps[pr[2, ]], LETTERS)
    sh <- c(sh, am[cbind(i1, j1)] + am[cbind(j1, i1)] > 0)
  }
  expect_gt(length(a), 0L)
  expect_identical(g$col$a, a)
  expect_identical(g$col$c, cc)
  expect_identical(g$col$b, b)
  expect_identical(g$col$shielded, sh)
  expect_identical(g$nedge, nrow(w))
})

test_that("dagcol() refuses what is not a DAG", {
  msg <- paste("dagcol: `from` and `to` must be labels of the same length,",
               "one edge each, at least one edge.")
  expect_error(dagcol(c("a", "b"), "c"), msg, fixed = TRUE)
  expect_error(dagcol(character(0), character(0)), msg, fixed = TRUE)
  expect_error(dagcol(list("a"), "b"), msg, fixed = TRUE)
  expect_error(dagcol(c("a", NA), c("b", "c")),
               "dagcol: `from` and `to` must have no missing values.",
               fixed = TRUE)
  cyc <- "dagcol: the edges have a cycle; a DAG has none."
  expect_error(dagcol(c("a", "b"), c("b", "a")), cyc, fixed = TRUE)
  expect_error(dagcol(c("a", "b", "c"), c("b", "c", "a")), cyc, fixed = TRUE)
  expect_error(dagcol("a", "a"), cyc, fixed = TRUE)
  tri <- "dagcol: `triple` must be three labels of nodes: a, c, b."
  expect_error(dagcol("a", "b", triple = c("a", "b")), tri, fixed = TRUE)
  expect_error(dagcol("a", "b", triple = c("a", "b", "z")), tri, fixed = TRUE)
  expect_error(dagcol("a", "b", triple = c("a", NA, "b")), tri, fixed = TRUE)
})
