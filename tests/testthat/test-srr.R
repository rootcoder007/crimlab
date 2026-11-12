# G2.2: each kind of vector argument refuses a matrix; one column is a vector.
test_that("a vector argument refuses a matrix", {
  m <- matrix(c(6, 4, 5, 5), 2)
  v <- "must be a vector, not a matrix"
  errs("vodtst", list(list(v, m, c(5, 5, 4, 6))))
  errs("sellab", list(list(v, c(1, 0), m > 4, c(TRUE, TRUE)),
                      list(v, matrix(c(1, 0, 1, 0), 2), rep(TRUE, 4),
                           rep(TRUE, 4))))
  errs("gtmp", list(list(v, c(0.1, 0.2, 0.1, 0.2), m / 20),
                    list(v, c(0.1, 0.2, 0.1, 0.2), hom = matrix(TRUE, 2, 2))))
  errs("mixlr", list(list(v, c(a = 0.1, b = 0.2), matrix("a", 2, 2), 1,
                          "a", 2)))
  errs("sntbnd", list(list(v, c(1, 0, 1, 0), matrix(c("a", "b"), 2, 2))))
  expect_identical(vodtst(matrix(c(6, 4)), c(5, 5)), vodtst(c(6, 4), c(5, 5)))
})

# The class, the number of rows (or elements), then each column or element
# name with the first letter of its type.
shape <- function(x) {
  paste(class(x)[1L], NROW(x),
        paste(names(x), substr(vapply(x, typeof, ""), 1L, 1L), sep = ":",
              collapse = " "))
}

#' @srrstats {EA6.0a} The class of each result is checked.
#' @srrstats {EA6.0b} The number of rows of each table is checked.
#' @srrstats {EA6.0c} The column and element names are checked.
#' @srrstats {EA6.0d} The type of every column and element is checked.
#' @srrstats {EA4.0} Numeric input gives double results; labels come back
#'   as character, signs as integer.
#' @srrstats {G5.3} Results hold no NA, NaN or Inf where none is documented.
test_that("results have the documented class, shape and types", {
  d <- c(a = 600, b = 400)
  n <- c(a = 700, b = 300)
  f <- "data.frame 2"
  x <- list(
    v = list(vodtst(d, n), f,
             "group:c ref:l shday:d shnight:d ratio:d sign:i"),
    o = list(orrisk(3, 0.4, 0.5), "numeric 4",
             "exposed:d unexposed:d rr:d or:d"),
    b = list(dspbm(d, n, c(9, 6)), f, "group:c ref:l cd:d fd:d rd:d"),
    y = list(dyadrt(outer(d, d) / 1e3, d), "data.frame 4",
             "off:c vic:c count:d rtoff:d rtpair:d nulloff:d ratio:d"),
    g = list(gtmp(c(0.1, 0.2), c(0.3, NA), c(FALSE, TRUE))$locus, f,
             "locus:c hom:l mp:d"),
    k = list(kinidx(c(1, 2), c(1, 2), c(0.1, 0.2))$locus, "data.frame 1",
             "locus:c idx:d"),
    c = list(clrsmp(c(90, 20), c(100, 50), c(19, 60), c(20, 140))$parts, f,
             "part:c ratea:d rateb:d sharea:d shareb:d rateeff:d mixeff:d"),
    e = list(expbnd(d, n, 1.5), f,
             "group:c ref:l rate:d lo:d hi:d ratio:d rlo:d rhi:d sign:i"),
    i = list(incap(c(2, 10), 0.1, 1), f,
             "group:c lam:d q:d s:d rate:d prev:d marg:d"),
    s = list(sntbnd(c(1, 0, 1, 0), c("x", "x", "y", "y"), "y")[1:2, ], f,
             "term:c share:d lower:d upper:d observed:d"),
    p = list(proof("orrisk"), "data.frame 3", "fn:c theorem:c file:c")
  )
  for (k in names(x)) {
    v <- x[[k]][[1L]]
    expect_identical(shape(v), paste(x[[k]][[2L]], x[[k]][[3L]]), info = k)
    num <- unlist(v[vapply(v, is.numeric, NA)])
    expect_true(all(is.finite(num)), info = k)
  }
})
