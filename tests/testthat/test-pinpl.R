# Point in polygon and the crosswalk. The kernels are checked by hand in
# src/t_geou.cpp; here the binding, against the rules written out in R.

sq <- function(x0, y0, id, s = 1, ring = 0) {
  data.frame(x = x0 + c(0, s, s, 0), y = y0 + c(0, 0, s, s), id = id,
             ring = ring)
}

test_that("points fall in their regions, holes and parts included", {
  poly <- rbind(sq(0, 0, "A", 4), sq(1, 1, "A", 1, ring = 1),
                sq(5, 0, "B"), sq(7, 0, "B", ring = 1), sq(3, 3, "C", 2))
  p <- cbind(x = c(0.5, 1.5, 5.5, 7.5, 9, 3.5, 4), y = c(0.5, 1.5, 0.5, 0.5, 9,
                                                         3.5, 2))
  u <- pinpl(p, poly = poly, pop = c(100, 50, 25))
  expect_s3_class(u, "pinpl")
  expect_named(u, c("reg", "nreg", "regs", "cov"))
  expect_identical(u$reg, c("A", NA, "B", "B", NA, "A", "A"))
  expect_identical(u$nreg, c(1L, 0L, 1L, 1L, 0L, 2L, 1L))
  expect_identical(u$regs$units, c(3, 2, 0))
  expect_identical(u$regs$hasunit, c(TRUE, TRUE, FALSE))
  expect_equal(u$regs$popshr, c(100, 50, 25) / 175, tolerance = 1e-15)
  expect_equal(unname(u$cov), c(3, 2, 5, 175, 150, 25, 150 / 175),
               tolerance = 1e-15)
  expect_output(print(u), "not a rate denominator")
  # against a ray-casting rule written out, on random points and a star
  set.seed(111)
  a <- seq(0, 2 * pi, length.out = 13)[-13]
  rr <- rep(c(1, 0.4), 6)
  star <- data.frame(lon = rr * cos(a), lat = rr * sin(a), id = 1)
  q <- cbind(x = runif(300, -1, 1), y = runif(300, -1, 1))
  inside <- function(px, py, vx, vy) {
    j <- c(length(vx), seq_len(length(vx) - 1))
    sum((vy > py) != (vy[j] > py) &
          px < vx + (py - vy) * (vx[j] - vx) / (vy[j] - vy)) %% 2 == 1
  }
  w <- mapply(inside, q[, 1], q[, 2], MoreArgs = list(star$lon, star$lat))
  s <- pinpl(q, poly = star)
  expect_identical(!is.na(s$reg), w)
  expect_true(all(is.na(s$regs$pop)))
  expect_true(is.na(s$cov[["covshr"]]))
  expect_output(print(s), "pinpl: ")
})

test_that("the crosswalk splits by overlap and aggregates exactly", {
  a <- arwt(c(a = 120, b = 40, c = 60, d = 100),
            from = c("a", "b", "c", "d", "d"),
            to = c("N1", "N1", "N2", "N2", "N3"), area = c(1, 1, 1, 7, 3),
            top = c(1, 2))
  expect_s3_class(a, "arwt")
  expect_identical(a$out$to, c("N1", "N2", "N3"))
  expect_equal(a$out$value, c(160, 130, 30), tolerance = 1e-15)
  expect_identical(a$map$to, c("N1", "N1", "N2", "N2"))
  expect_equal(a$map$share, c(1, 1, 1, 0.7), tolerance = 1e-15)
  expect_false(a$exact)
  expect_identical(a$lost, 0)
  v <- c(30, 130, 160)
  g <- (3 + 1 - 2 * sum(cumsum(v)) / sum(v)) / 3
  expect_equal(unname(a$conc), c(g, 160 / 320, 290 / 320), tolerance = 1e-15)
  expect_named(a$conc, c("gini", "top1", "top2"))
  e <- arwt(c(5, 7, 9), from = c("x", "y"), to = c("P", "P"),
            id = c("x", "y", "z"))
  expect_true(e$exact)
  expect_identical(e$lost, 9)
  expect_identical(e$out$value, 12)
  expect_true(is.na(e$map$to[3]))
  expect_output(print(e), "exact")
  expect_output(print(a), "split by overlap")
  z <- arwt(c(a = 0, b = 0), from = c("a", "b"), to = c("P", "Q"))
  expect_identical(z$conc[["gini"]], 0)
  expect_true(is.na(z$conc[["top10"]]))
})

test_that("pinpl() and arwt() refuse bad input in words", {
  poly <- sq(0, 0, "A")
  p <- cbind(x = 0.5, y = 0.5)
  errs("pinpl", list(
    list("data frame of x, y and id", p, poly = 1:3),
    list("data frame of x, y and id", p, poly = poly[, c("x", "y")]),
    list("at least 3 vertices", p, poly = poly[1:2, ]),
    list("`pop`", p, poly = poly, pop = c(1, 2))
  ))
  errs("arwt", list(
    list("`id` must name each source", c(1, 2), "a", "P"),
    list("each source unit once", c(1, 2), "a", "P", id = c("a", "a")),
    list("not in `id`", c(a = 1), "b", "P"),
    list("`from` must hold at least one pair", c(a = 1), character(0),
         character(0)),
    list("`x` must be finite numbers", numeric(0), "a", "P",
         id = character(0)),
    list("`area`", c(a = 1), "a", "P", area = -1),
    list("`top`", c(a = 1), "a", "P", top = 0)
  ))

})
