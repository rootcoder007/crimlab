# relchk(): every reading of a table against a base-R enumeration written
# here from the textbook formulas, verdicts, rounding, zero cells, the input
# forms and the refusals.

tab <- matrix(c(20, 10, 80, 90), 2,
              dimnames = list(c("A", "B"), c("yes", "no")))

# every reading of x, by the formulas on the four cells, in no set order
allr <- function(x) {
  out <- list()
  for (grp in c("rows", "cols")) {
    y <- if (grp == "rows") x else t(x)
    for (ty in c("or", "rr", "rd")) {
      for (h in rownames(y)) {
        for (g in setdiff(rownames(y), h)) {
          for (e in colnames(y)) {
            a <- y[g, e]
            b <- sum(y[g, ]) - a
            c <- y[h, e]
            d <- sum(y[h, ]) - c
            v <- switch(ty, or = (a * d) / (b * c),
                        rr = (a / (a + b)) / (c / (c + d)),
                        rd = a / (a + b) - c / (c + d))
            out[[length(out) + 1L]] <- data.frame(
              grp = grp, type = ty, group = g, ref = h, event = e, value = v
            )
          }
        }
      }
    }
  }
  do.call(rbind, out)
}

# readings keyed and sorted, so that two enumerations compare
keyed <- function(r) {
  r <- r[c("grp", "type", "group", "ref", "event", "value")]
  r <- r[order(paste(r$grp, r$type, r$group, r$ref, r$event)), ]
  rownames(r) <- NULL
  r
}

test_that("every reading of a 2 x 2 table is the textbook value", {
  r <- relchk(tab, 2.25)
  expect_identical(nrow(r$reads), 24L)
  expect_equal(keyed(r$reads), keyed(allr(tab)), tolerance = 1e-14)
  # by hand
  v <- function(g, h, e, ty, grp = "rows") {
    s <- r$reads
    s$value[s$grp == grp & s$type == ty & s$group == g & s$ref == h &
              s$event == e]
  }
  expect_equal(v("A", "B", "yes", "or"), 2.25, tolerance = 1e-15)
  expect_equal(v("B", "A", "yes", "or"), 4 / 9, tolerance = 1e-15)
  expect_equal(v("A", "B", "no", "or"), 4 / 9, tolerance = 1e-15)
  expect_equal(v("B", "A", "no", "or"), 2.25, tolerance = 1e-15)
  expect_equal(v("A", "B", "yes", "rr"), 2, tolerance = 1e-15)
  expect_equal(v("B", "A", "no", "rr"), 1.125, tolerance = 1e-15)
  expect_equal(v("A", "B", "yes", "rd"), 0.1, tolerance = 1e-15)
  expect_equal(v("yes", "no", "A", "rr", "cols"), 17 / 12, tolerance = 1e-15)
  expect_identical(r$reads$match, abs(r$reads$value - 2.25) <= 0.005)
})

test_that("every relabelling of a 3 x 3 table gives the same readings", {
  x <- matrix(c(12, 7, 30, 5, 22, 9, 40, 3, 16), 3,
              dimnames = list(c("p", "q", "s"), c("lo", "mid", "hi")))
  base <- relchk(x, 1.2, digits = 1)
  expect_equal(keyed(base$reads), keyed(allr(x)), tolerance = 1e-14)
  hits <- function(r) {
    m <- r$reads[r$reads$match, ]
    sort(paste(m$grp, m$type, m$group, m$ref, m$event))
  }
  perms <- list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2),
                c(3, 2, 1))
  for (pr in perms) {
    for (pc in perms) {
      r <- relchk(x[pr, pc], 1.2, digits = 1)
      expect_equal(keyed(r$reads), keyed(base$reads), tolerance = 1e-14)
      expect_identical(hits(r), hits(base))
    }
  }
})

test_that("the odds ratio inverts with the reference or the event", {
  x <- matrix(c(7, 13, 29, 11), 2, dimnames = list(c("g", "h"), c("e", "f")))
  r <- relchk(x, 1)$reads
  o <- r[r$grp == "rows" & r$type == "or", ]
  ab <- o$value[o$group == "g" & o$event == "e"]
  expect_equal(o$value[o$group == "h" & o$event == "e"], 1 / ab,
               tolerance = 1e-14)
  expect_equal(o$value[o$group == "g" & o$event == "f"], 1 / ab,
               tolerance = 1e-14)
  expect_equal(o$value[o$group == "h" & o$event == "f"], ab,
               tolerance = 1e-14)
  expect_equal(ab, (7 * 11) / (29 * 13), tolerance = 1e-15)
})

test_that("the verdict says whether the stated labelling reproduces it", {
  r <- relchk(tab, 2.25, ref = "B", event = "yes")
  expect_identical(r$verdict, "holds")
  expect_true(r$ok)
  expect_identical(r$claim, c(group = NA, ref = "B", event = "yes"))
  expect_identical(sum(r$reads$claim), 1L)
  # 0.44 is the same contrast with the reference swapped
  r <- relchk(tab, 0.44, ref = "B", event = "yes")
  expect_identical(r$verdict, "other")
  expect_false(r$ok)
  m <- r$reads[r$reads$match & r$reads$grp == "rows", ]
  expect_setequal(paste(m$group, m$ref, m$event), c("B A yes", "A B no"))
  # 2.0 at one decimal is the risk ratio, not the odds ratio
  r <- relchk(tab, 2.0, type = "or", digits = 1)
  expect_identical(r$verdict, "other")
  expect_identical(unique(r$reads$type[r$reads$match]), "rr")
  expect_identical(relchk(tab, 2.0, type = "rr", digits = 1)$verdict, "holds")
  expect_identical(relchk(tab, 7)$verdict, "none")
  # by position, a factor by its label, and group as well as ref
  expect_identical(relchk(tab, 2.25, ref = 2, event = 1)$verdict, "holds")
  expect_identical(relchk(tab, 2.25, ref = factor("B"))$verdict, "holds")
  expect_identical(relchk(tab, 2.25, group = "B", ref = "A")$verdict, "holds")
  expect_identical(relchk(tab, 2.25, group = "B", ref = "A",
                          event = "yes")$verdict, "other")
  expect_identical(relchk(tab, 0.1, type = "rd", group = "A")$verdict,
                   "holds")
})

test_that("labels that differ by case are different groups", {
  x <- matrix(c(20, 10, 80, 90), 2, dimnames = list(c("a", "A"), c("y", "n")))
  expect_identical(relchk(x, 2.25, ref = "A")$verdict, "holds")
  expect_identical(relchk(x, 2.25, ref = "a")$verdict, "holds")
  expect_identical(relchk(x, 2.25, ref = "a", event = "y")$verdict, "other")
  expect_error(relchk(x, 2.25, ref = "b"), "`ref` names no row: b.")
})

test_that("a rounded value matches at either end of its last decimal", {
  # odds ratio 469 / 200 = 2.345 for row 1 against row 2, column 1
  x <- matrix(c(469, 1, 200, 1), 2)
  m <- function(v, ...) {
    r <- relchk(x, v, ...)$reads
    r$match[r$grp == "rows" & r$type == "or" & r$group == "1" &
              r$event == "1"]
  }
  expect_true(m(2.35))
  expect_true(m(2.34))
  expect_false(m(2.33))
  expect_false(m(2.36))
  expect_true(m(2.3, digits = 1))
  expect_false(m(2.3, digits = 2))
  expect_true(m(2.345, digits = 3))
  expect_false(m(2.345 + 1e-12, tol = 0))
  # just above and below the half
  y <- matrix(c(23451, 1, 10000, 1), 2)
  expect_true(relchk(y, 2.35, ref = 2, event = 1)$ok)
  expect_false(relchk(y, 2.34, ref = 2, event = 1)$ok)
  y <- matrix(c(23449, 1, 10000, 1), 2)
  expect_true(relchk(y, 2.34, ref = 2, event = 1)$ok)
  expect_false(relchk(y, 2.35, ref = 2, event = 1)$ok)
  # the decimals taken by default, and the tolerance they give
  expect_identical(relchk(x, 2.3)$digits, 1L)
  expect_equal(relchk(x, 2.3)$tol, 0.05, tolerance = 1e-15)
  expect_identical(relchk(x, 2)$digits, 0L)
  expect_identical(relchk(x, 1 / 3)$digits, NA_integer_)
  expect_identical(relchk(x, 1 / 3)$tol, 0)
  expect_identical(relchk(x, 2.3, tol = 0.01)$digits, NA_integer_)
  expect_identical(relchk(x, 2.3, tol = 0.01)$tol, 0.01)
})

test_that("zero cells give Inf or NaN and never match", {
  x <- matrix(c(0, 0, 5, 5, 5, 0), 3)  # rows 1, 2: 0 and 5; row 3: 5 and 0
  r <- relchk(x, 0, tol = 1e300)$reads
  expect_equal(keyed(r), keyed(allr(`dimnames<-`(x, list(c("1", "2", "3"),
                                                         c("1", "2"))))),
               tolerance = 1e-14)
  s <- r[r$grp == "rows" & r$group == "1" & r$ref == "2" & r$event == "1", ]
  expect_identical(s$value, c(NaN, NaN, 0))
  expect_true(all(is.infinite(r$value[r$grp == "rows" & r$type == "or" &
                                        r$group == "3" & r$event == "1"])))
  expect_false(any(r$match[!is.finite(r$value)]))
  expect_true(all(r$match[is.finite(r$value)]))
})

test_that("the table forms give one answer", {
  r <- relchk(tab, 2.25)
  df <- data.frame(yes = c(20, 10), no = c(80, 90), row.names = c("A", "B"))
  expect_identical(relchk(df, 2.25), r)
  expect_identical(relchk(as.table(tab), 2.25), r)
  l <- list(yes = c(A = 20, B = 10), no = c(80, 90))
  expect_identical(relchk(l, 2.25), r)
  storage.mode(tab) <- "integer"
  expect_identical(relchk(tab, 2.25), r)
  # counts from records, through xtabs()
  d <- data.frame(g = rep(c("A", "B"), each = 100),
                  y = c(rep(c("yes", "no"), c(20, 80)),
                        rep(c("yes", "no"), c(10, 90))))
  x <- xtabs(~ g + y, d)
  expect_identical(relchk(x, 2.25, ref = "B", event = "yes")$verdict, "holds")
  # unnamed rows and columns are "1", "2", ...
  u <- relchk(unname(tab), 2.25)$reads
  expect_setequal(u$group, c("1", "2"))
  # weighted counts are read as they are
  w <- relchk(tab / 3, 2.25)
  expect_equal(w$reads$value, r$reads$value, tolerance = 1e-14)
})

test_that("relchk() refuses what it cannot read", {
  big <- matrix(1, 1000, 1000)
  errs("relchk", list(
    list("`x` must be counts", "a", 1),
    list("`x` must be counts", array(1, c(2, 2, 2)), 1),
    list("`x` must be counts", matrix("1", 2, 2), 1),
    list("`x` must be counts", data.frame(a = factor(1:2), b = 1:2), 1),
    list("`x` must be counts", list(a = 1:2, b = 1:3), 1),
    list("`x` must be counts", list(a = matrix(1, 2, 2), b = 1:2), 1),
    list("has value labels", structure(tab, labels = c(a = 1)), 1),
    list("has value labels",
         list(a = structure(1:2, labels = c(a = 1)), b = 1:2), 1),
    list("at least two rows and two columns", matrix(1:3, 1), 1),
    list("at least two rows and two columns", list(), 1),
    list("too large", big, 1),
    list("missing count", matrix(c(1, NA, 1, 1), 2), 1),
    list("finite counts", matrix(c(1, Inf, 1, 1), 2), 1),
    list("negative count", matrix(c(1, -1, 1, 1), 2), 1),
    list("no count in row B", tab * c(1, 0), 1),
    list("no count in column no", t(t(tab) * c(1, 0)), 1),
    list("names a row twice: A", `rownames<-`(tab, c("A", "A")), 1),
    list("names a column twice: y", `colnames<-`(tab, c("y", "y")), 1),
    list("missing row label", `rownames<-`(tab, c("A", NA)), 1),
    list("`value` must be finite", tab, NA_real_),
    list("`value` must have length 1", tab, c(1, 2)),
    list("`type` must be", tab, 1, "OR"),
    list("`type` must be", tab, 1, c("or", "rr")),
    list("not both", tab, 1, digits = 1, tol = 0.1),
    list("`digits` must be at most 15", tab, 1, digits = 16),
    list("`digits` must be one whole number", tab, 1, digits = -1),
    list("`tol` must be at least 0", tab, 1, tol = -1),
    list("two different rows", tab, 1, group = "A", ref = 1),
    list("`event` names no column: maybe", tab, 1, event = "maybe"),
    list("`event` must be one column name or a position from 1 to 2", tab,
         1, event = 3),
    list("`group` must be one row name", tab, 1, group = NA)
  ))
})

test_that("print.relchk() states the claim, the verdict and the matches", {
  expect_output(print(relchk(tab, 2.25, ref = "B", event = "yes")),
                paste0("to 2 decimals\nclaim: or, ref B, event yes\n",
                       "verdict: the claim holds"))
  expect_output(print(relchk(tab, 0.4, ref = "B")), "to 1 decimal\n")
  expect_output(print(relchk(tab, 0.44, tol = 0.01)),
                "within 0.01\nclaim: or\nverdict: the claim holds")
  expect_output(print(relchk(tab, 2.0, type = "or", digits = 1)),
                "only readings the claim does not allow match")
  out <- capture.output(print(relchk(tab, 7)))
  expect_identical(out[3], "verdict: no reading matches")
  expect_length(out, 3L)
  expect_output(print(relchk(matrix(c(0, 0, 5, 5, 5, 0), 3), 1)),
                "20 of 54 readings undefined")
  r <- relchk(tab, 7)
  expect_output(expect_invisible(print(r)))
})
