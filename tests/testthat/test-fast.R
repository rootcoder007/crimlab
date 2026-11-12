# The fast paths of the binding: names matched by identity when every name
# carries one encoding mark (by text otherwise), and group labels written
# directly. Each must give what the general path gives, refusals included.

test_that("names of one mark match by identity, mixed marks by text", {
  # all ASCII: one mark
  a <- c(x = 1, y = 2, z = 3)
  b <- c(z = 30, x = 10, y = 20)
  expect_identical(vodtst(a, b), vodtst(a, c(x = 10, y = 20, z = 30)))
  # the very same names: nothing to match
  expect_identical(vodtst(a, a), vodtst(a, unname(a)))
  expect_identical(outtst(c(p = 0.2, q = 0.8), 0.5, c(p = 1, q = 3)),
                   outtst(c(0.2, 0.8), 0.5, c(1, 3)))
  # all marked UTF-8 on one side, latin1 on the other: the second is looked
  # up by its text
  u <- c("\u00e9t\u00e9", "\u00fcber")
  l <- iconv(u, "UTF-8", "latin1")
  d <- stats::setNames(c(1, 2), u)
  e <- stats::setNames(c(4, 3), rev(l))
  expect_identical(vodtst(d, e), vodtst(d, stats::setNames(c(3, 4), u)))
  expect_identical(vodtst(e, d)[["group"]], rev(l))
  # mixed marks within one argument: matched by text throughout
  m <- c("x", u[1])
  n <- c(l[1], "x")
  expect_identical(
    vodtst(stats::setNames(c(1, 2), m), stats::setNames(c(5, 6), n)),
    vodtst(stats::setNames(c(1, 2), m), stats::setNames(c(6, 5), m))
  )
})

test_that("the fast paths refuse what the general path refuses", {
  u <- "\u00e9t\u00e9"
  l <- iconv(u, "UTF-8", "latin1")
  # the same text under two marks is one name, given twice
  expect_error(vodtst(stats::setNames(c(1, 1), c(u, l)), c(1, 1)),
               "names a group twice")
  expect_error(vodtst(stats::setNames(c(1, 1), c("x", "x")), c(1, 1)),
               "names a group twice")
  expect_error(vodtst(c(x = 1, y = 1), c(x = 1, z = 1)),
               "missing: y; extra: z")
  expect_error(
    vodtst(stats::setNames(c(1, 1), c(u, "y")),
           stats::setNames(c(1, 1), c(l, "z"))),
    "missing: y; extra: z"
  )
  # a "bytes" string equals only a "bytes" string
  b <- "\xe9"
  Encoding(b) <- "bytes"
  d <- stats::setNames(c(1, 2), c(b, "x"))
  expect_identical(vodtst(d, stats::setNames(c(4, 3), c("x", b))),
                   vodtst(d, c(3, 4)))
  expect_error(vodtst(stats::setNames(c(1, 1), c(b, "x")),
                      stats::setNames(c(1, 1), c(enc2utf8("\u00e9"), "x"))))
})

test_that("unnamed groups are labelled 1, 2, ... at every size", {
  for (n in c(1, 9, 10, 64, 65, 1000)) {
    x <- rep(1, n)
    expect_identical(vodtst(x, x)[["group"]], as.character(seq_len(n)))
    expect_identical(incap(x + 1, 0.5, 1, share = x / n)[["group"]],
                     c(as.character(seq_len(n)), "all"))
  }
  n <- 1e5
  expect_identical(incap(rep(2, n), 0.5, 1, share = rep(1 / n, n))[["group"]],
                   c(as.character(seq_len(n)), "all"))
  expect_identical(
    incap(c(a = 2, b = 3), 0.5, 1, share = c(0.5, 0.5))[["group"]],
    c("a", "b", "all")
  )
  expect_identical(incap(c(2, 3), 0.5, 1, share = c(b = 0.5, c = 0.5)),
                   incap(c(b = 2, c = 3), 0.5, 1, share = c(0.5, 0.5)))
})

test_that("numbers as incident labels give what their text gives", {
  off <- c("a", "b", "c", "a", "d", "b", "c")
  inc <- c(3L, 3L, 1L, 1L, 7L, 7L, 7L)
  ref <- coofg(as.character(inc), off)
  expect_identical(coofg(inc, off), ref)
  expect_identical(coofg(as.double(inc), off), ref)
  # -0 is the label "0", as 0 is
  expect_identical(coofg(c(0, -0, 1, 1, 2, 2, 2), off),
                   coofg(c("0", "0", "1", "1", "2", "2", "2"), off))
  # fractions and large numbers are read through their text
  x <- c(0.5, 0.5, 1e15, 1e15, 2.25, 2.25, 2.25)
  expect_identical(coofg(x, off), coofg(as.character(x), off))
  expect_error(coofg(c(1, NA, 2, 2, 3, 3, 3), off), "missing values")
  expect_error(coofg(c(1L, NA, 2L, 2L, 3L, 3L, 3L), off), "missing values")
  expect_error(coofg(matrix(1:8, 2), c(off, "e")), "a vector")
  expect_identical(coofg(factor(c("x", "x", "y", "y", "z", "z", "z")), off),
                   coofg(c("x", "x", "y", "y", "z", "z", "z"), off))
})

test_that("edge labels are coded alike by address and by text", {
  f <- c("a", "b", "c", "a")
  t <- c("b", "c", "a", "d")
  r <- ntdeg(f, t)
  expect_identical(r[["node"]], c("a", "b", "c", "d"))
  expect_identical(r, ntdeg(f, t, nodes = c("a", "b", "c", "d")))
  # labels under two marks: the same nodes, read by their text
  u <- "été"
  l <- iconv(u, "UTF-8", "latin1")
  expect_identical(ntdeg(c(u, "x"), c("x", l))[["deg"]], c(1L, 1L))
  expect_identical(ntdeg(1:3, c(2L, 3L, 1L)), ntdeg(c("1", "2", "3"),
                                                    c("2", "3", "1")))
})
