# dwnexp() against the three rules written out here in base R: bearings by
# atan2() on the plane and the initial great-circle bearing on the sphere,
# the circular mean of the wind, the angle about the circle and plainly,
# the haversine distance, and a polygon meeting a wedge by its vertices,
# the apex, its sides and its arc; hand cases; and the refusals.

# the bearing from (x0, y0) to (x1, y1), degrees clockwise from north
rbrg <- function(x0, y0, x1, y1, ll = FALSE) {
  if (ll) {
    r <- pi / 180
    b <- atan2(sin((x1 - x0) * r) * cos(y1 * r),
               cos(y0 * r) * sin(y1 * r) -
                 sin(y0 * r) * cos(y1 * r) * cos((x1 - x0) * r))
  } else {
    b <- atan2(x1 - x0, y1 - y0)
  }
  (b * 180 / pi) %% 360
}

rdst <- function(x0, y0, x1, y1, ll = FALSE) {
  if (!ll) return(sqrt((x1 - x0)^2 + (y1 - y0)^2))
  r <- pi / 180
  h <- sin((y1 - y0) * r / 2)^2 +
    cos(y0 * r) * cos(y1 * r) * sin((x1 - x0) * r / 2)^2
  2 * 6371.0088 * asin(sqrt(pmin(1, h)))
}

rdif <- function(a, b) {
  d <- abs(a - b) %% 360
  ifelse(d > 180, 360 - d, d)
}

test_that("the sector rule counts the sources a place lies downwind of", {
  set.seed(31)
  s <- data.frame(x = runif(40, -100, 100), y = runif(40, -100, 100),
                  day = sample(1:4, 40, TRUE), dir = runif(40, 0, 360))
  p <- data.frame(x = runif(25, -150, 150), y = runif(25, -150, 150))
  for (ll in c(FALSE, TRUE)) {
    sx <- if (ll) s$x / 2 else s$x
    sy <- if (ll) s$y / 3 else s$y
    px <- if (ll) p$x / 2 else p$x
    py <- if (ll) p$y / 3 else p$y
    r <- dwnexp(sx, sy, s$day, s$dir, px, py, half = 25, rad = 120,
                lonlat = ll)
    to <- (s$dir + 180) %% 360
    want <- unlist(lapply(1:4, function(d) {
      k <- which(s$day == d)
      vapply(seq_len(25), function(j) {
        b <- rbrg(sx[k], sy[k], px[j], py[j], ll)
        dd <- rdst(sx[k], sy[k], px[j], py[j], ll)
        sum(rdif(b, to[k]) <= 25 & dd <= 120)
      }, 0)
    }))
    expect_identical(r$n, as.numeric(want), info = ll)
    expect_identical(r$any, as.integer(want > 0))
    expect_identical(r$per, rep(1:4, each = 25))
    expect_identical(r$id, rep(1:25, 4))
  }
})

test_that("the sector rule: hand cases and the count by month", {
  # the wind from the south blows to the north: (0, 100) in, (50, 100) at
  # 26.57 degrees in, (60, 100) at 30.96 out, the source's point in
  r <- dwnexp(0, 0, NULL, 180, c(0, 50, 60, 0), c(100, 100, 100, 0))
  expect_identical(r$any, c(1L, 1L, 0L, 1L))
  # from = FALSE takes 180 as the direction blown to
  q <- dwnexp(0, 0, NULL, 180, c(0, 0), c(100, -100), from = FALSE)
  expect_identical(q$any, c(0L, 1L))
  # two plumes reach the first place on day 1, one on day 2: the
  # county-month count of downwind days is 2 (Bao and Da 2026, eq 1)
  s <- dwnexp(c(0, 1, 0), c(0, 0, 0), c(1, 1, 2), c(180, 180, 180), 0, 50)
  expect_identical(s$n, c(2, 1))
  expect_identical(sum(s$any), 2L)
})

# the edge e = (x0, y0, x1, y1) crosses the circle of radius rad about
# (ax, ay) at a point within hf of th
rarc <- function(e, ax, ay, th, hf, rad) {
  dx <- e[3] - e[1]
  dy <- e[4] - e[2]
  fx <- e[1] - ax
  fy <- e[2] - ay
  a <- dx^2 + dy^2
  b <- 2 * (fx * dx + fy * dy)
  disc <- b^2 - 4 * a * (fx^2 + fy^2 - rad^2)
  if (disc < 0) return(FALSE)
  t <- (-b + c(-1, 1) * sqrt(disc)) / (2 * a)
  t <- t[t >= 0 & t <= 1]
  any(rdif(rbrg(ax, ay, e[1] + t * dx, e[2] + t * dy), th) <= hf)
}

# a point in the wedge from (ax, ay) about th, half-angle hf, radius rad
rinw <- function(px, py, ax, ay, th, hf, rad) {
  d <- sqrt((px - ax)^2 + (py - ay)^2)
  d == 0 | (d <= rad & rdif(rbrg(ax, ay, px, py), th) <= hf)
}

# whether segments p and q (x0, y0, x1, y1) cross
rsg <- function(p, q) {
  o <- function(a, b, c) {
    sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
  }
  o(p[1:2], p[3:4], q[1:2]) != o(p[1:2], p[3:4], q[3:4]) &&
    o(q[1:2], q[3:4], p[1:2]) != o(q[1:2], q[3:4], p[3:4])
}

# a polygon (x, y) meets the wedge, written out: a vertex in it, the apex
# in the polygon (even-odd), or an edge across a side or the arc
rwedge <- function(x, y, ax, ay, th, hf, rad) {
  if (any(rinw(x, y, ax, ay, th, hf, rad))) return(TRUE)
  n <- length(x)
  j <- c(n, seq_len(n - 1))
  cross <- (y > ay) != (y[j] > ay) &
    ax < (x[j] - x) * (ay - y) / (y[j] - y) + x
  if (sum(cross) %% 2 == 1) return(TRUE)
  far <- max(sqrt((x - ax)^2 + (y - ay)^2))
  len <- if (is.finite(rad)) rad else 2 * far + 1
  sd <- lapply(c(th - hf, th + hf), function(a) {
    c(ax, ay, ax + len * sin(a * pi / 180), ay + len * cos(a * pi / 180))
  })
  hit <- vapply(seq_len(n), function(i) {
    e <- c(x[j[i]], y[j[i]], x[i], y[i])
    rsg(e, sd[[1]]) || rsg(e, sd[[2]]) ||
      (is.finite(rad) && rarc(e, ax, ay, th, hf, rad))
  }, TRUE)
  any(hit)
}

test_that("the wedge rule is the polygon test written out", {
  set.seed(32)
  cx <- runif(60, -150, 150)
  cy <- runif(60, -150, 150)
  w <- runif(60, 5, 40)
  sq <- data.frame(id = rep(1:60, each = 4),
                   x = rep(cx, each = 4) + c(-1, 1, 1, -1) * rep(w, each = 4),
                   y = rep(cy, each = 4) + c(-1, -1, 1, 1) * rep(w, each = 4))
  for (rad in c(100, Inf)) {
    r <- dwnexp(c(0, 20), c(0, -10), c(1, 2), c(200, 120), sq$x, sq$y,
                sq$id, rule = "wedge", rad = rad)
    to <- c(20, 300)
    want <- unlist(lapply(1:2, function(d) {
      vapply(1:60, function(i) {
        k <- sq$id == i
        rwedge(sq$x[k], sq$y[k], c(0, 20)[d], c(0, -10)[d], to[d], 30, rad)
      }, TRUE)
    }))
    expect_identical(r$any, as.integer(want), info = rad)
  }
  # the hand case of the help page: a square north of the source
  h <- dwnexp(0, 0, NULL, 180, c(-10, 10, 10, -10), c(40, 40, 60, 60),
              rep(1, 4), rule = "wedge", rad = 100)
  expect_identical(h$any, 1L)
})

test_that("the angle rule is Ayesh's instrument written out", {
  set.seed(33)
  s <- data.frame(x = runif(12, 70, 74), y = runif(12, 30, 32),
                  yr = rep(2015:2016, each = 6), dir = runif(12, 200, 280))
  p <- data.frame(x = c(s$x[1:3], runif(5, 69, 76)),
                  y = c(s$y[1:3], runif(5, 29, 33)))
  for (wrap in c(TRUE, FALSE)) {
    r <- dwnexp(s$x, s$y, s$yr, s$dir, p$x, p$y, rule = "angle",
                lonlat = TRUE, wrap = wrap)
    want <- do.call(rbind, lapply(2015:2016, function(yy) {
      k <- s$yr == yy
      to <- (s$dir[k] + 180) * pi / 180
      md <- (atan2(sum(sin(to)), sum(cos(to))) * 180 / pi) %% 360
      b <- rbrg(mean(s$x[k]), mean(s$y[k]), p$x, p$y, TRUE)
      a <- if (wrap) rdif(b, md) else abs(b - md)
      d <- vapply(seq_len(8), function(j) {
        dd <- rdst(s$x[k], s$y[k], p$x[j], p$y[j], TRUE)
        mean(dd[dd > 0])
      }, 0)
      cbind(a, d)
    }))
    expect_equal(r$angle, unname(want[, 1]), tolerance = 1e-12)
    expect_equal(r$dist, unname(want[, 2]), tolerance = 1e-12)
    expect_identical(r$z, r$angle * r$dist)
  }
  q <- dwnexp(s$x, s$y, s$yr, s$dir, p$x, p$y, rule = "angle", dwt = "div")
  expect_identical(q$z, q$angle / q$dist)
  q <- dwnexp(s$x, s$y, s$yr, s$dir, p$x, p$y, rule = "angle", dwt = "none")
  expect_identical(q$z, q$angle)
})

test_that("dwnexp() refusals name the function", {
  e <- function(m, ...) expect_error(dwnexp(...), m, fixed = TRUE)
  e("dwnexp: `rule` must be one of", 0, 0, NULL, 0, 1, 1, rule = "cone")
  e("dwnexp: `sx`, `sy` and `sdir` must have one length.", c(0, 1), 0,
    NULL, 0, 1, 1)
  e("dwnexp: `sper` must give the period of each source, with no NA.",
    0, 0, NA, 0, 1, 1)
  e("dwnexp: `rx` and `ry` must have one length.", 0, 0, NULL, 0, 1,
    c(1, 2))
  e("dwnexp: `rid` must give the place of each row of `rx`, with no NA.",
    0, 0, NULL, 0, 1, 1, rid = c(1, 2))
  e("dwnexp: `half` must be between 0 and 180.", 0, 0, NULL, 0, 1, 1,
    half = 200)
  e("dwnexp: `rad` must be a positive number or Inf.", 0, 0, NULL, 0, 1,
    1, rad = 0)
  e("dwnexp: `rid` must name each place once.", 0, 0, NULL, 0, c(1, 2),
    c(1, 2), rid = c(1, 1))
  e("dwnexp: `rule = \"wedge\"` needs `rid`, the polygons.", 0, 0, NULL,
    0, c(1, 2, 3), c(1, 2, 3), rule = "wedge")
  e("dwnexp: each polygon in `rid` needs three vertices or more.", 0, 0,
    NULL, 0, c(1, 2), c(1, 2), rid = c(1, 1), rule = "wedge")
  e("dwnexp: `rule = \"wedge\"` needs planar coordinates (`lonlat = FALSE`).",
    0, 0, NULL, 0, c(1, 2, 3), c(1, 2, 3), rid = c(1, 1, 1),
    rule = "wedge", lonlat = TRUE)
  e("dwnexp: `from` must be TRUE or FALSE.", 0, 0, NULL, 0, 1, 1,
    from = NA)
  e("dwnexp: `sdir` must be finite, with no missing values.", 0, 0, NULL,
    NA_real_, 1, 1)
})

test_that("tailwd() is the speed times the cosine off the course", {
  set.seed(34)
  d <- runif(30, 0, 360)
  s <- runif(30, 0, 8)
  expect_equal(tailwd(d, s), s * cos((d + 180 - 250) * pi / 180),
               tolerance = 1e-14)
  # hand values: from 70 at 4 m/s a 4 m/s tailwind; from 340 or 160 none;
  # from 250 a 4 m/s headwind
  expect_equal(tailwd(c(70, 340, 160, 250), rep(4, 4)), c(4, 0, 0, -4),
               tolerance = 1e-14)
  # by group: the maximum (Hener 2022) or the mean, NA readings left out
  g <- rep(c("b", "a", "c"), each = 10)
  d[3] <- NA
  m <- tailwd(d, s, by = g)
  t1 <- s * cos((d + 180 - 250) * pi / 180)
  expect_identical(m$by, c("a", "b", "c"))
  expect_equal(m$tw, unname(c(max(t1[11:20]), max(t1[1:10], na.rm = TRUE),
                              max(t1[21:30]))), tolerance = 1e-14)
  a <- tailwd(d, s, course = 70, by = g, agg = "mean")
  t2 <- s * cos((d + 180 - 70) * pi / 180)
  expect_equal(a$tw, unname(c(mean(t2[11:20]), mean(t2[1:10], na.rm = TRUE),
                              mean(t2[21:30]))), tolerance = 1e-14)
  expect_identical(tailwd(NA_real_, 1, by = 1)$tw, NA_real_)
  expect_error(tailwd(1, c(1, 2)),
               "tailwd: `dir` and `spd` must be numeric vectors of one length.",
               fixed = TRUE)
  expect_error(tailwd(1, 1, by = c(1, 2)),
               "tailwd: `by` must give the group of each reading, with no NA.",
               fixed = TRUE)
  expect_error(tailwd(1, 1, agg = "min"), "tailwd: `agg` must be one of",
               fixed = TRUE)
})

test_that("tailwd() and dwnexp() refuse speeds below 0 and latitudes past 90", {
  expect_error(tailwd(10, -2),
               "tailwd: `spd` must hold wind speeds, 0 or more (or NA).",
               fixed = TRUE)
  # the tailwind of a speed 0 is 0, of NA is NA
  expect_identical(tailwd(c(10, 10), c(0, NA)), c(0, NA))
  m <- paste("dwnexp: with `lonlat = TRUE`, `sy` and `ry` must be latitudes,",
             "-90 to 90.")
  expect_error(dwnexp(0, 100, NULL, 180, c(0, 10), c(80, 10), lonlat = TRUE),
               m, fixed = TRUE)
  expect_error(dwnexp(0, 10, NULL, 180, c(0, 10), c(80, -91), lonlat = TRUE),
               m, fixed = TRUE)
})

test_that("labels that are not vectors and directions past 1e9 are refused", {
  errs("dwnexp", list(
    list("`sper` must give the period of each source", 0, 0, function() 1,
         90, 1, 1),
    list("`sdir` must be between -1e\\+09 and 1e\\+09", 0, 0, NULL, 1e308,
         1, 1)
  ))
  errs("tailwd", list(
    list("`by` must give the group of each reading", c(10, 20), c(1, 2),
         by = list(1, 2))
  ))
})
