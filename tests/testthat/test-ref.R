# The stored reference values (ref/*.txt) and their tolerance table
# (ref/tol.csv): every value has the row dev/ref/cmp.R checks it by when
# the ref workflow regenerates it on another platform, and every row is
# used. Needs none of the packages that computed the values.

test_that("every stored reference value has a tolerance", {
  tab <- utils::read.csv(test_path("ref", "tol.csv"), stringsAsFactors = FALSE)
  expect_named(tab, c("fixture", "path", "kind", "tol", "test"))
  expect_true(all(tab$kind %in% c(
    "equal", "each", "maxrel", "maxabs", "exact", "report"
  )))
  expect_true(all(is.finite(tab$tol) & tab$tol >= 0))
  expect_false(anyDuplicated(tab[c("fixture", "path")]) > 0)
  fx <- sub("\\.txt$", "", list.files(test_path("ref"), "\\.txt$"))
  expect_setequal(unique(tab$fixture), fx)
  for (f in fx) {
    x <- dget(test_path("ref", paste0(f, ".txt")))
    expect_type(x$pkgs, "character")
    l <- names(refleaf(x[setdiff(names(x), "pkgs")]))
    p <- tab$path[tab$fixture == f]
    i <- refrow(p, l)
    expect_false(anyNA(i), info = paste(f, l[is.na(i)][1]))
    expect_setequal(unique(i), seq_along(p))
  }
})

test_that("leaves and rows are found by path", {
  x <- list(a = 1, b = list(c = "x", 2:3), d = data.frame(e = 1))
  expect_named(refleaf(x), c("a", "b/c", "b/2", "d/e"))
  expect_identical(refleaf(x)[["b/2"]], 2:3)
  p <- c("b", "b/c", "*/e", "z")
  expect_identical(refrow(p, c("a", "b/c", "b/2", "d/e")),
                   c(NA, 2L, 1L, 3L))
})
