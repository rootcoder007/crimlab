# The spatial functions. Their kernels are checked in C++ (src/t_geo.cpp,
# t_geo2.cpp, t_cgr.cpp, t_nnd.cpp, t_rpk.cpp, t_hsp.cpp, t_stm.cpp,
# t_jtc.cpp); here the binding: the input forms, the conversions, the
# messages, and agreement with other R implementations where one exists.
# Each comparison states its tolerance; references come from base R, and
# from geosphere, sf and spatstat as values stored in ref/sp.txt by
# dev/ref/sp.R (versions named there).

rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

test_that("coordinates are read by name, never by position when named", {
  x <- c(3, 1, 4, 1, 5)
  y <- c(9, 2, 6, 5, 3)
  r <- cgstat(x, y)
  expect_identical(cgstat(cbind(x, y)), r)
  expect_identical(cgstat(cbind(y = y, x = x)), r)
  expect_identical(cgstat(data.frame(id = 1:5, y = y, x = x)), r)
  expect_identical(cgstat(unname(cbind(x, y))), r)
  expect_identical(cgstat(cbind(as.integer(x), as.integer(y))), r)
  # a classed column (units, say) is read as its numbers
  u <- structure(x, class = "units")
  expect_identical(cgstat(data.frame(x = I(u), y = y)), r)
  # one-column matrices are vectors
  expect_identical(cgstat(matrix(x), matrix(y)), r)
})

test_that("spatial input is refused when it cannot be read safely", {
  m <- cbind(a = 1:3, b = 1:3)
  ll <- data.frame(lon = c(-79, -78), lat = c(43, 44))
  errs("cgstat", list(
    list("`y` is missing", c(1, 2)),
    list("columns named x and y", m),
    list("columns named x and y", matrix(1:9, 3)),
    list("set `lonlat = TRUE`", ll),
    list("`lonlat` must be TRUE or FALSE", ll, lonlat = NA),
    list("`lon` must be longitudes", c(-181, 0), c(0, 0), lonlat = TRUE),
    list("`lat` must be latitudes", c(1, 0), c(0, 91), lonlat = TRUE),
    list("antipode of their centre", c(0, 0, 180), c(0, 0, 0),
         lonlat = TRUE),
    list("`x` must be finite", c(1, NA), c(1, 2)),
    list("`y` must be finite", c(1, 2), c(1, Inf)),
    list("`y` must have length 2", c(1, 2), c(1, 2, 3)),
    list("`x` must be finite", data.frame(x = I(list(1, 2)), y = 1:2)),
    list("`x` must be finite", character(0), character(0)),
    list("`w` must have a positive sum", 1:2, 1:2, w = c(0, 0)),
    list("`w` must be at least 0", 1:2, 1:2, w = c(-1, 2)),
    list("`k` must be 0, 1 or 2", 1:3, 1:3, k = 3),
    list("must exceed `k`", 1:2, 1:2, k = 2)
  ))
  errs("cvhull", list(list("magnitude", c(1e-200, 1, 2), c(0, 1, 2))))
  errs("hvdist", list(
    list("columns named lon and lat", cbind(x = 1, y = 2), c(0, 0)),
    list("same number of points", rbind(c(0, 0), c(1, 1)),
         rbind(c(0, 0), c(1, 1), c(2, 2))),
    list("`r` must be one number", c(0, 0), c(1, 1), r = 0)
  ))
  errs("eqproj", list(
    list("both `lon0` and `lat0`", 1:2, 1:2, lon0 = 1),
    list("`lat0` must be one number from -90 to 90", 1:2, 1:2, lon0 = 1,
         lat0 = 100),
    list("antipode of the centre", c(0, 1), c(0, 0), lon0 = -180, lat0 = 0)
  ))
  errs("eqinv", list(
    list("farther than 2 r", 2e7, 0, lon0 = 0, lat0 = 0),
    list("this function takes planar coordinates", ll, lon0 = 0, lat0 = 0)
  ))
  errs("webpt", list(
    list("out of range", c(-1, 1), c(0, 0), start = c(200, 0),
         lonlat = TRUE),
    list("antipode of their centre", c(-1, 1), c(0, 0), start = c(180, 0),
         lonlat = TRUE),
    list("`start` must have length 2", 1:3, 1:3, start = 1)
  ))
})

#' @srrstats {SP6.2} Distances from pole to pole, the projection's round
#'   trip at both poles (src/t_geo2.cpp too) and a centre of points
#'   around the north pole are tested.
#' @noRd
test_that("great-circle distances agree with geosphere and closed forms", {
  r <- 6371008.8
  expect_equal(hvdist(c(0, 0), c(1, 0)), r * pi / 180, tolerance = 1e-15)
  # SP6.2: the poles
  expect_equal(hvdist(c(10, 90), c(-150, -90)), r * pi, tolerance = 1e-15)
  expect_lt(abs(hvdist(c(10, 90), c(-150, 90))), 1e-9)
  set.seed(3)
  a <- cbind(lon = runif(100, -180, 180), lat = runif(100, -90, 90))
  b <- cbind(lon = runif(100, -180, 180), lat = runif(100, -90, 90))
  # geosphere: distVincentySphere(a, b, r = r) and (a, b[1, ], r = r)
  ref <- refget("sp")$gc
  expect_equal(sum(a, b), ref$inp, tolerance = 1e-12, info = "inputs changed")
  expect_lt(rel(hvdist(a, b), ref$g), 1e-12)
  expect_lt(rel(hvdist(as.data.frame(a), b[1, ]), ref$g1), 1e-12)
})

#' @srrstats {SP6.0} The projection and its inverse return the original
#'   longitudes and latitudes, poles included.
#' @srrstats {SP2.2b} eqproj() agrees with the same projection done by sf
#'   (PROJ; its coordinates stored with the tests), so results move
#'   between the two.
#' @noRd
test_that("the projection inverts and agrees with PROJ through sf", {
  set.seed(4)
  lon <- runif(200, -100, -60)
  lat <- runif(200, 20, 60)
  p <- eqproj(lon, lat)
  a <- attr(p, "centre")
  b <- eqinv(p, lon0 = a[1], lat0 = a[2])
  expect_lt(max(abs(b$lon - lon)), 1e-10)
  expect_lt(max(abs(b$lat - lat)), 1e-10)
  q <- eqproj(cbind(lon = c(0, 120, -60), lat = c(90, -90, 0)), lon0 = 0,
              lat0 = 0)
  expect_equal(eqinv(q$x, q$y, lon0 = 0, lat0 = 0)$lat, c(90, -90, 0),
               tolerance = 1e-12)
  # sf: st_transform() of the points to "+proj=laea" about the centre a
  # (lat_0 and lon_0 written with %.17g), +R=6371008.8
  ref <- refget("sp")$proj
  expect_equal(sum(lon, lat), ref$inp, tolerance = 1e-12,
               info = "inputs changed")
  expect_lt(max(abs(ref$x - p$x)), 1e-6)
  expect_lt(max(abs(ref$y - p$y)), 1e-6)
})

test_that("centrography agrees with base R", {
  set.seed(5)
  x <- rnorm(100, 50, 10)
  y <- runif(100, 1, 9)
  w <- rpois(100, 3) + 1
  r <- cgstat(x, y)
  expect_equal(unname(r[c("w", "mx", "my", "sx", "sy")]),
               c(100, mean(x), mean(y), sd(x), sd(y)), tolerance = 1e-13)
  expect_equal(r[["sdd"]], sqrt(var(x) + var(y)), tolerance = 1e-13)
  expect_equal(r[["gx"]], exp(mean(log(x))), tolerance = 1e-13)
  expect_equal(r[["hy"]], 1 / mean(1 / y), tolerance = 1e-13)
  # integer weights are replicated points
  xr <- rep(x, w)
  yr <- rep(y, w)
  expect_equal(cgstat(x, y, w = w), cgstat(xr, yr), tolerance = 1e-13)
  expect_equal(cgstat(x, y, k = 0)[["sx"]], sqrt(mean((x - mean(x))^2)),
               tolerance = 1e-13)
  expect_identical(cgmedn(x, y), c(x = median(x), y = median(y)))
  e <- sdell(x, y, c = 2)
  ev <- eigen(cov(cbind(x, y)), symmetric = TRUE)
  expect_equal(unname(e[c("a", "b")]), 2 * sqrt(ev$values),
               tolerance = 1e-12)
  v <- ev$vectors[, 1] * sign(ev$vectors[2, 1])
  expect_equal(e[["angle"]], atan2(v[2], v[1]), tolerance = 1e-12)
  expect_equal(e[["area"]], pi * e[["a"]] * e[["b"]], tolerance = 1e-15)
  h <- cvhull(x, y)
  expect_setequal(h$idx, grDevices::chull(x, y))
  i <- h$idx
  j <- c(i[-1L], i[1L])
  expect_equal(h$area, sum(x[i] * y[j] - x[j] * y[i]) / 2, tolerance = 1e-13)
  expect_equal(h$perim, sum(sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2)),
               tolerance = 1e-13)
})

test_that("the centre of minimum distance is the known minimiser", {
  # the equilateral triangle: its centre; a square: its centre
  r <- webpt(c(0, 1, 0.5), c(0, 0, sqrt(3) / 2))
  expect_equal(unname(r[c("x", "y")]), c(0.5, sqrt(3) / 6), tolerance = 1e-9)
  r <- webpt(cbind(x = c(0, 2, 2, 0), y = c(0, 0, 2, 2)))
  expect_equal(unname(r[c("x", "y", "obj")]), c(1, 1, 4 * sqrt(2)),
               tolerance = 1e-9)
  # a point holding more than half the weight is the minimiser, certified
  r <- webpt(c(0, 5, 9), c(0, 1, 7), w = c(3, 1, 1), start = c(4, 4))
  expect_identical(unname(r[c("x", "y", "gap", "code")]), c(0, 0, 0, 2))
  # the objective is not larger than at any data point or the mean
  set.seed(6)
  x <- rexp(60)
  y <- rexp(60)
  r <- webpt(x, y)
  f <- function(a, b) sum(sqrt((x - a)^2 + (y - b)^2))
  expect_lte(r[["obj"]], min(mapply(f, c(x, mean(x)), c(y, mean(y)))))
  expect_equal(r[["obj"]], f(r[["x"]], r[["y"]]), tolerance = 1e-13)
  expect_lt(r[["gap"]], 1e-6)
  expect_identical(webpt(x, y, maxit = 0)[["code"]], 0)
})

#' @srrstats {SP6.1} Functions with `lonlat` are tested on planar and on
#'   longitude-latitude input.
#' @srrstats {SP6.1b} With `lonlat = TRUE` the result equals projecting
#'   first with eqproj(), computing in metres and mapping centres back with
#'   eqinv().
#' @srrstats {SP6.1a} Treating degrees as planar coordinates gives a wrong
#'   standard distance; the projected one matches great-circle distances.
#' @noRd
test_that("longitude and latitude are projected, and only on request", {
  set.seed(7)
  lon <- runif(80, -79.6, -79.2)
  lat <- runif(80, 43.6, 43.8)
  p <- eqproj(lon, lat)
  a <- attr(p, "centre")
  for (f in list(cgstat, sdell, webpt)) {
    g <- f(lon, lat, lonlat = TRUE)
    h <- f(p$x, p$y)
    i <- 1:2 + identical(f, cgstat)
    b <- eqinv(h[[i[1]]], h[[i[2]]], lon0 = a[1], lat0 = a[2])
    expect_equal(unname(g[i]), c(b$lon, b$lat), tolerance = 1e-13)
  }
  g <- cgstat(cbind(lon = lon, lat = lat), lonlat = TRUE)
  expect_true(all(is.nan(g[c("gx", "gy", "hx", "hy")])))
  m <- cgmedn(lon, lat, lonlat = TRUE)
  expect_equal(m[[2L]], median(lat), tolerance = 1e-4)
  expect_identical(cvhull(lon, lat, lonlat = TRUE)$idx, cvhull(p)$idx)
  # the standard distance in metres matches great-circle distances to the
  # mean centre; in degrees it is off by five orders of magnitude
  d <- hvdist(cbind(lon, lat), g[c("mx", "my")])
  expect_equal(g[["sdd"]], sqrt(sum(d^2) / 79), tolerance = 1e-4)
  expect_gt(g[["sdd"]] / cgstat(lon, lat)[["sdd"]], 5e4)
  # SP6.2: points around the north pole have their centre near it (the
  # planar mean in the projection about latitude 89, not the spherical
  # centre, so within 1e-4 degrees, not exactly)
  q <- cgstat(c(0, 90, 180, -90), rep(89, 4), lonlat = TRUE)
  expect_lt(90 - q[["my"]], 1e-4)
  r <- webpt(c(-1, 1, 0), c(0, 0, 1), start = c(0, 0.5), lonlat = TRUE)
  expect_gt(r[["code"]], 0)
})

test_that("nearest neighbours and K agree with spatstat", {
  set.seed(8)
  x <- runif(150)
  y <- 2 * runif(150)
  w <- c(0, 1, 0, 2)
  r <- nnidx(cbind(x, y), win = w, kmax = 3, dist = TRUE)
  expect_named(r, c("order", "ce", "dist"))
  expect_identical(dim(r$dist), c(150L, 3L))
  expect_equal(r$order$mean, colMeans(r$dist), tolerance = 1e-14)
  expect_null(nnidx(x, y, win = w)$dist)
  rr <- seq(0, 0.3, by = 0.01)
  k <- ripk(x, y, win = w, r = rr)
  expect_equal(k$l, sqrt(k$k / pi), tolerance = 1e-15)
  # spatstat: X is ppp(x, y, c(0, 1), c(0, 2)); nndist(X, k = 1:3),
  # clarkevans(X, correction = c("none", "Donnelly")) and Kest(X, r = rr)
  # with each correction
  ref <- refget("sp")$nn
  expect_equal(sum(x, y), ref$inp, tolerance = 1e-12, info = "inputs changed")
  expect_equal(r$dist, ref$dist, tolerance = 1e-14)
  expect_equal(r$order$ratio[1], ref$naive, tolerance = 1e-13)
  expect_equal(r$ce[["ratio"]], ref$donnelly, tolerance = 1e-13)
  expect_equal(k$k, ref$iso, tolerance = 1e-12)
  expect_equal(ripk(x, y, w, rr, "border")$k, ref$border, tolerance = 1e-12)
  expect_equal(ripk(x, y, w, rr, "none")$k, ref$none, tolerance = 1e-12)
  # the modified border estimator asked for alone (see dev/doss/dst.R)
  expect_equal(ripk(x, y, w, rr, "modif")$k, ref$modif, tolerance = 1e-12)
})

test_that("simulations follow the seed and detect clustering", {
  set.seed(9)
  x <- runif(60)
  y <- runif(60)
  w <- c(0, 1, 0, 1)
  rr <- seq(0.02, 0.2, 0.02)
  set.seed(1)
  a <- ripenv(x, y, w, rr, nsim = 19)
  set.seed(1)
  expect_identical(ripenv(x, y, w, rr, nsim = 19), a)
  expect_identical(ripenv(x, y, w, rr, nsim = 19, seed = 4),
                   ripenv(x, y, w, rr, nsim = 19, seed = 4))
  expect_equal(a$env$obs, ripk(x, y, w, rr)$l - rr, tolerance = 1e-14)
  expect_true(all(a$env$lo <= a$env$mean & a$env$mean <= a$env$hi))
  expect_equal(a$test[["t"]], max(abs(a$env$obs)), tolerance = 1e-15)
  # tight clusters: the global test and the index reject randomness
  cx <- rep(c(0.2, 0.8, 0.5), each = 20) + rnorm(60, sd = 0.01)
  cy <- rep(c(0.3, 0.7, 0.5), each = 20) + rnorm(60, sd = 0.01)
  expect_identical(ripenv(cx, cy, w, rr, nsim = 19, seed = 2)$test[["p"]],
                   0.05)
  n <- nnidx(cx, cy, w, nsim = 19, seed = 3)$order
  expect_identical(n$pclu, 0.05)
  expect_identical(n$preg, 1)
  expect_lt(n$ratio, 0.2)
})

test_that("distance functions refuse what they cannot use", {
  w <- c(0, 1, 0, 1)
  errs("nnidx", list(
    list("at least 2 points", 0.5, 0.5, w),
    list("inside `win`", c(0.5, 2), c(0.5, 0.5), w),
    list("with xmin < xmax", c(0.5, 0.6), c(0.5, 0.5), c(1, 0, 0, 1)),
    list("`win` must have length 4", c(0.5, 0.6), c(0.5, 0.5), c(0, 1)),
    list("`kmax` must be below", c(0.5, 0.6), c(0.5, 0.5), w, kmax = 2),
    list("`nsim` must be one whole number", c(0.5, 0.6), c(0.5, 0.5), w,
         nsim = -1),
    list("`dist` must be TRUE or FALSE", c(0.5, 0.6), c(0.5, 0.5), w,
         dist = 1),
    list("`seed` must be one whole number", c(0.5, 0.6), c(0.5, 0.5), w,
         nsim = 1, seed = -1),
    list("takes planar coordinates", data.frame(lon = 1:2, lat = 1:2),
         win = w)
  ))
  errs("ripk", list(
    list("strictly increasing", c(0.5, 0.6), c(0.5, 0.5), w, c(0.1, 0.1)),
    list("`r` must be at least 0", c(0.5, 0.6), c(0.5, 0.5), w, -1),
    list("`corr` must be one of", c(0.5, 0.6), c(0.5, 0.5), w, 0.1,
         corr = "ripley")
  ))
  errs("ripenv", list(
    list("`nsim` must be one whole number from 1", c(0.5, 0.6), c(0.5, 0.5),
         w, 0.1, nsim = 0)
  ))
})

#' @srrstats {SP2.3} Points written to a GeoPackage and read back by sf
#'   are taken as input (the objects sf made are stored with the tests).
#' @srrstats {SP2.2b} Results go back to sf: projected points become sf
#'   points with sf::st_as_sf() (its coordinates stored).
#' @noRd
test_that("sf points go in, with their CRS checked against lonlat", {
  set.seed(10)
  d <- data.frame(lon = runif(20, -80, -79), lat = runif(20, 43, 44))
  # sf objects: s is st_as_sf(d, coords = c("lon", "lat"), crs = 4326), g
  # s written with st_write() to a GeoPackage and read back with
  # st_read(), p st_transform(s, 32617), pg st_geometry(p), xy
  # st_coordinates(p), ex the X of st_as_sf(eqproj(d), coords = c("x",
  # "y")), u two POINT Z, l a LINESTRING
  ref <- refget("sp")$obj
  expect_equal(sum(d), ref$inp, tolerance = 1e-12, info = "inputs changed")
  s <- ref$s
  expect_identical(cgstat(s, lonlat = TRUE), cgstat(d, lonlat = TRUE))
  expect_identical(hvdist(s, c(-79.5, 43.5)), hvdist(d, c(-79.5, 43.5)))
  g <- ref$g
  expect_equal(cgstat(g, lonlat = TRUE), cgstat(d, lonlat = TRUE),
               tolerance = 1e-15)
  p <- ref$p
  xy <- ref$xy
  expect_identical(cgstat(p), cgstat(xy))
  expect_identical(cgstat(ref$pg), cgstat(xy))
  expect_equal(ref$ex, eqproj(d)$x, tolerance = 1e-8)
  u <- ref$u
  l <- ref$l
  errs("cgstat", list(
    list("geographic CRS: set `lonlat = TRUE`", s),
    list("projected CRS", p, lonlat = TRUE),
    list("POINT geometries", l),
    list("two-dimensional points", u)
  ))
  errs("nnidx", list(list("geographic CRS: this function takes planar", s,
                          win = c(-80, -79, 43, 44))))
  errs("hvdist", list(list("projected CRS", p, c(0, 0))))
})

test_that("every spatial function answers one valid call in its shape", {
  w <- c(0, 1, 0, 1)
  x <- c(0.1, 0.4, 0.8, 0.3)
  y <- c(0.2, 0.9, 0.5, 0.6)
  calls <- list(
    list(quote(hvdist(c(0, 0), c(1, 1))), NULL, "numeric 1"),
    list(quote(eqproj(c(0, 1), c(0, 1))), c("x", "y"), "data.frame 2"),
    list(quote(eqinv(c(0, 1), c(0, 1), lon0 = 0, lat0 = 0)),
         c("lon", "lat"), "data.frame 2"),
    list(quote(cgstat(x, y)), c("w", "mx", "my", "sx", "sy", "sdd", "gx",
                               "gy", "hx", "hy"), "numeric 10"),
    list(quote(cgmedn(x, y)), c("x", "y"), "numeric 2"),
    list(quote(webpt(x, y)), c("x", "y", "obj", "gap", "iter", "code"),
         "numeric 6"),
    list(quote(sdell(x, y)), c("mx", "my", "a", "b", "angle", "area"),
         "numeric 6"),
    list(quote(cvhull(x, y)), c("idx", "area", "perim"), "list 3"),
    list(quote(nnidx(x, y, w)), c("order", "ce", "dist"), "list 3"),
    list(quote(ripk(x, y, w, c(0.1, 0.2))), c("r", "k", "l"),
         "data.frame 2"),
    list(quote(ripenv(x, y, w, c(0.1, 0.2), nsim = 3, seed = 1)),
         c("env", "test"), "list 2"),
    list(quote(fzmode(x, y, 0.5)), NULL, "integer 4"),
    list(quote(kmpart(x, y, 2, seed = 1)), c("cl", "cen", "fit"), "list 3"),
    list(quote(nnhclu(x, y, fixd = 0.5, minpts = 2)), c("cl", "clusters"),
         "list 2"),
    list(quote(stacel(x, y, w, r = 0.5, minpts = 2)), c("cl", "clusters"),
         "list 2"),
    list(quote(mantst(x, y, t = 1:4, nsim = 3, seed = 1)),
         c("r", "pairs", "pup", "plo", "p2"), "numeric 5"),
    list(quote(stmavg(x, y, t = 1:4, span = 3)), c("mx", "my", "n"),
         "data.frame 4"),
    list(quote(jtdecay(1:2, "negexp", c(1, 1))), NULL, "numeric 2"),
    list(quote(jtsurf(x, y, w, 2, 2, "negexp", c(1, 1))),
         c("raw", "prob", "peak"), "list 3"),
    list(quote(jtcal(x, 1:3, 1)), NULL, "numeric 3"),
    list(quote(jtacc(matrix(1:4, 2), w, c(0.5, 0.5))),
         c("p", "hitbest", "hitworst", "err", "peak", "cell"), "numeric 6"),
    list(quote(bjtc(1:2, diag(2), 1)), c("pj", "po", "pod", "prod", "risk"),
         "data.frame 2")
  )
  fns <- vapply(calls, function(a) as.character(a[[1L]][[1L]]), "")
  expect_setequal(fns, spf)
  for (k in seq_along(calls)) {
    v <- eval(calls[[k]][[1L]])
    expect_named(v, calls[[k]][[2L]])
    expect_identical(paste(class(v)[1L], NROW(v)), calls[[k]][[3L]],
                     info = fns[k])
    num <- unlist(v[vapply(v, is.numeric, NA)])
    expect_true(all(is.finite(num[!is.na(num)])), info = fns[k])
  }
  expect_type(cvhull(x, y)$idx, "integer")
  expect_type(kmpart(x, y, 2, seed = 1)$cl, "integer")
  expect_type(nnhclu(x, y, fixd = 0.5, minpts = 2)$cl, "integer")
})

#' @srrstats {SP6.5} kmpart() is compared with stats::kmeans(), the
#'   non-spatial K-means of base R, from the same starting centres.
#' @noRd
test_that("hot-spot functions agree with brute force and base R", {
  set.seed(11)
  x <- c(rnorm(40, 0, 0.5), rnorm(40, 4, 0.5), runif(20, -2, 6))
  y <- c(rnorm(40, 0, 0.5), rnorm(40, 4, 0.5), runif(20, -2, 6))
  d <- as.matrix(dist(cbind(x, y)))
  expect_identical(fzmode(x, y, r = 0.7), as.integer(rowSums(d <= 0.7)))
  expect_identical(fzmode(c(1, 1, 2), c(1, 1, 2), r = 0), c(2L, 2L, 1L))
  s <- cbind(x = x[c(1, 50, 90)], y = y[c(1, 50, 90)])
  r <- kmpart(x, y, k = 3, init = s)
  b <- stats::kmeans(cbind(x, y), s, iter.max = 100, algorithm = "Lloyd")
  expect_identical(r$cl, unname(b$cluster))
  expect_equal(as.matrix(r$cen[, c("x", "y")]), unname(b$centers),
               tolerance = 1e-13, ignore_attr = TRUE)
  expect_equal(r$fit[["ss"]], b$tot.withinss, tolerance = 1e-12)
  expect_equal(r$cen$ss, b$withinss, tolerance = 1e-12)
  expect_identical(r$cen$n, as.numeric(b$size))
  w <- kmpart(x, y, k = 3, init = s, upd = "weber")
  # each Weber centre is its group's centre of minimum distance
  for (j in 1:3) {
    i <- w$cl == j
    expect_equal(w$cen$sd[j], webpt(x[i], y[i])[["obj"]], tolerance = 1e-8)
  }
  expect_identical(kmpart(x, y, 2, seed = 3), kmpart(x, y, 2, seed = 3))
  g <- kmpart(x, y, k = 2, init = "cells", cell = c(1, 3))
  expect_setequal(round(g$cen$x), c(0, 4))
  # two tight groups, a fixed threshold: two clusters at order 1
  h <- nnhclu(c(rep(0, 12), rep(10, 12)) + rep(seq(0, 0.11, 0.01), 2),
              rep(0, 24), fixd = 0.5, minpts = 5, maxord = 2)
  expect_identical(h$clusters$points, c(12, 12))
  expect_identical(h$cl[, 1], rep(1:2, each = 12))
  expect_equal(h$clusters$mx, c(0.055, 10.055), tolerance = 1e-14)
  expect_true(all(is.na(h$cl[, 2])))
  # the random-expectation threshold, p = 0.5: 0.5 sqrt(A / n)
  q <- nnhclu(x, y, area = 64, minpts = 5, maxord = 1)
  expect_equal(q$clusters$thresh[1], 0.5 * sqrt(64 / 100), tolerance = 1e-14)
  e <- stacel(x, y, win = c(-2, 6, -2, 6), r = 0.75, top = 10, minpts = 5)
  expect_identical(nrow(e$clusters), length(unique(na.omit(e$cl))))
  for (j in seq_len(nrow(e$clusters))) {
    i <- which(e$cl == j)
    expect_equal(e$clusters$mx[j], mean(x[i]), tolerance = 1e-13)
    expect_identical(e$clusters$points[j], length(i) + 0)
  }
  # lon/lat: the same as on the projected points, centres back in degrees
  lon <- -79 + x / 100
  lat <- 43 + y / 100
  p <- eqproj(lon, lat)
  expect_identical(fzmode(lon, lat, r = 500, lonlat = TRUE),
                   fzmode(p, r = 500))
  k1 <- kmpart(lon, lat, 2, seed = 1, lonlat = TRUE)
  k2 <- kmpart(p, k = 2, seed = 1)
  expect_identical(k1$cl, k2$cl)
  a <- attr(p, "centre")
  expect_equal(k1$cen$x, eqinv(k2$cen$x, k2$cen$y, lon0 = a[1],
                               lat0 = a[2])$lon, tolerance = 1e-13)
  k3 <- kmpart(lon, lat, 2, init = cbind(lon = lon[c(1, 50)],
                                         lat = lat[c(1, 50)]), lonlat = TRUE)
  expect_identical(k3$cl, kmpart(p, k = 2, init = as.matrix(p[c(1, 50), ]))$cl)
  n1 <- nnhclu(lon, lat, area = 64e6, minpts = 5, lonlat = TRUE)
  n2 <- nnhclu(p, area = 64e6, minpts = 5)
  expect_identical(n1$cl, n2$cl)
})

test_that("hot-spot functions refuse what they cannot use", {
  x <- c(0, 1, 2, 3)
  errs("fzmode", list(list("`r` must be one number from 0", x, x, r = -1)))
  errs("kmpart", list(
    list("`k` must be at most", x, x, k = 5),
    list("`init` must be one of", x, x, k = 2, init = "random"),
    list("`k` starting centres", x, x, k = 2, init = cbind(x = 1, y = 1)),
    list("`cell` must be finite", x, x, k = 2, init = "cells"),
    list("positive cell side", x, x, k = 2, init = "cells", cell = c(0, 1)),
    list("`upd` must be one of", x, x, k = 2, upd = "median"),
    list("fewer distinct points", c(1, 1, 2), c(1, 1, 2), k = 3),
    list("fewer than `k` cells", x, x, k = 2, init = "cells",
         cell = c(10, 0))
  ))
  errs("nnhclu", list(
    list("at least 2 points", 1, 1, fixd = 1),
    list("give one of `area`", x, x),
    list("give one of `area`", x, x, area = 1, fixd = 1),
    list("`p` must be one number strictly between 0 and 1", x, x, fixd = 1,
         p = 1),
    list("`maxord` must be at most 20", x, x, fixd = 1, maxord = 21),
    list("`k` must be 0, 1 or 2", x, x, fixd = 1, k = 3)
  ))
  errs("stacel", list(
    list("at most 10000", x, x, win = c(0, 3, 0, 3), nx = 20000, r = 1),
    list("`k` must be 0, 1 or 2", x, x, win = c(0, 3, 0, 3), r = 1, k = 3)
  ))
})

test_that("space-time functions agree with their definitions", {
  set.seed(12)
  x <- runif(40)
  y <- runif(40)
  t <- 10 * x + rnorm(40)
  r <- mantst(x, y, t, nsim = 99, seed = 1)
  ds <- as.vector(dist(cbind(x, y)))
  dt <- as.vector(dist(t))
  expect_equal(r[["r"]], cor(ds, dt), tolerance = 1e-12)
  expect_identical(r[["pairs"]], 780)
  expect_identical(r[["pup"]], 0.01)
  expect_equal(mantst(x, y, t, cs = 0.1, ct = 1, nsim = 0)[["r"]],
               cor(1 / (ds + 0.1), 1 / (dt + 1)), tolerance = 1e-12)
  expect_true(is.nan(mantst(x, y, t, nsim = 0)[["p2"]]))
  m <- stmavg(x, y, t, span = 3)
  o <- order(t)
  i <- o[c(4, 5, 6)]
  expect_equal(m$mx[o[5]], mean(x[i]), tolerance = 1e-15)
  expect_identical(m$n[o[1]], 2)
  expect_identical(m$n[o[2]], 3)
  # lon/lat: the means of the projected points, back in degrees
  lon <- -79 + x / 10
  lat <- 43 + y / 10
  p <- eqproj(lon, lat)
  a <- attr(p, "centre")
  q <- stmavg(lon, lat, t, span = 5, lonlat = TRUE)
  s <- stmavg(p, t = t, span = 5)
  expect_equal(q$mx, eqinv(s$mx, s$my, lon0 = a[1], lat0 = a[2])$lon,
               tolerance = 1e-13)
  expect_equal(mantst(lon, lat, t, nsim = 0, lonlat = TRUE)[["r"]],
               mantst(p, t = t, nsim = 0)[["r"]], tolerance = 1e-14)
})

test_that("journey-to-crime functions agree with their definitions", {
  d <- c(0, 0.3, 1, 2.5, 7)
  expect_equal(jtdecay(d, "negexp", c(2, 0.7)), 2 * exp(-0.7 * d),
               tolerance = 1e-14)
  expect_equal(jtdecay(d, "normal", c(2, 1, 0.5)), 2 * dnorm(d, 1, 0.5),
               tolerance = 1e-14)
  expect_equal(jtdecay(d, "lognormal", c(2, 0.1, 0.5)),
               2 * dlnorm(d, 0.1, 0.5), tolerance = 1e-14)
  expect_equal(jtdecay(d, "linear", c(1, -0.4)), pmax(0, 1 - 0.4 * d),
               tolerance = 1e-15)
  expect_equal(jtdecay(d, "trunc", c(1, 0.5, 0.8)),
               ifelse(d <= 0.5, d / 0.5, exp(-0.8 * (d - 0.5))),
               tolerance = 1e-14)
  tab <- c(0.2, 1, 0.5, 0.1)
  expect_equal(jtdecay(d, "table", 1, tab),
               approx(0:3, tab, d, yleft = 0, yright = 0)$y, tolerance = 1e-15)
  # the surface by its definition over the cell centres
  ix <- c(2, 3, 4.5)
  iy <- c(2, 4, 3)
  s <- jtsurf(ix, iy, win = c(0, 6, 0, 4), nx = 12, ny = 8, type = "negexp",
              par = c(1, 1))
  cx <- (seq_len(12) - 0.5) * 0.5
  cy <- (seq_len(8) - 0.5) * 0.5
  g <- outer(cx, cy, function(a, b) {
    exp(-sqrt((a - ix[1])^2 + (b - iy[1])^2)) +
      exp(-sqrt((a - ix[2])^2 + (b - iy[2])^2)) +
      exp(-sqrt((a - ix[3])^2 + (b - iy[3])^2))
  })
  expect_equal(s$raw, g, tolerance = 1e-14)
  expect_equal(s$prob, g / sum(g), tolerance = 1e-14)
  k <- which.max(g)
  expect_identical(unname(s$peak), c(k, cx[(k - 1) %% 12 + 1],
                                     cy[(k - 1) %/% 12 + 1]))
  l <- jtsurf(ix, iy, win = c(0, 6, 0, 4), nx = 12, ny = 8, type = "negexp",
              par = c(1, 1), mode = "lik", metric = "manhattan")
  ll <- outer(cx, cy, function(a, b) {
    -(abs(a - ix[1]) + abs(b - iy[1])) - (abs(a - ix[2]) + abs(b - iy[2])) -
      (abs(a - ix[3]) + abs(b - iy[3]))
  })
  expect_equal(l$raw, ll, tolerance = 1e-14)
  expect_equal(l$prob, exp(ll - max(ll)) / sum(exp(ll - max(ll))),
               tolerance = 1e-13)
  # calibration: the reflected normal kernel estimate by hand
  tr <- c(0.5, 1, 1.2, 3)
  tt <- seq(0, 5, 0.5)
  w <- c(1, 2, 1, 1)
  h <- 0.4
  f <- vapply(tt, function(u) {
    sum(w * (dnorm(u, tr, h) + dnorm(u, -tr, h))) / sum(w)
  }, 0)
  expect_equal(jtcal(tr, tt, h, w = w), f, tolerance = 1e-14)
  expect_equal(jtcal(tr, -1, h, refl = FALSE), mean(dnorm(-1, tr, h)),
               tolerance = 1e-14)
  # accuracy on a grid with ties
  pr <- matrix(c(1, 3, 3, 0, 2, 3), 3, 2)
  a <- jtacc(pr, win = c(0, 3, 0, 2), base = c(1.5, 1.5))
  expect_identical(unname(a[c("p", "hitbest", "hitworst", "peak", "cell")]),
                   c(2, 4 / 6, 4 / 6, 2, 5))
  b <- jtacc(pr, win = c(0, 3, 0, 2), base = c(2.5, 0.5))
  expect_identical(unname(b[c("hitbest", "hitworst")]), c(1 / 6, 3 / 6))
  expect_equal(b[["err"]], sqrt((2.5 - 1.5)^2), tolerance = 1e-15)
  # the Bayesian blocks by hand
  od <- matrix(c(5, 1, 0, 2, 4, 1, 0, 1, 6), 3)
  pj <- c(0.5, 0.3, 0.2)
  r <- bjtc(pj, od, dest = c(1, 2, 2))
  po <- rowSums(od) / sum(od)
  pod <- od[, 1] + 2 * od[, 2]
  pod <- pod / sum(pod)
  pr <- pj * pod / sum(pj * pod)
  rk <- pr / po
  expect_equal(r$po, po, tolerance = 1e-15)
  expect_equal(r$pod, pod, tolerance = 1e-15)
  expect_equal(r$prod, pr, tolerance = 1e-15)
  expect_equal(r$risk, rk / sum(rk), tolerance = 1e-15)
})

test_that("space-time and journey-to-crime functions refuse bad input", {
  x <- c(0, 1, 2, 3)
  errs("mantst", list(
    list("from 3 to 46340", 1:2, 1:2, t = 1:2),
    list("`cs` must be one number from 0", x, x, t = x, cs = -1),
    list("all equal", x, x, t = rep(1, 4)),
    list("`t` must have length 4", x, x, t = 1:3)
  ))
  errs("stmavg", list(list("`span` must be odd", x, x, t = x, span = 2)))
  errs("jtdecay", list(
    list("`type` must be one of", 1, "exp", 1),
    list("`par` must have length 2", 1, "negexp", 1),
    list("at least 2 values", 1, "table", 1, tab = 1),
    list("outside the domain", 1, "negexp", c(1, -1))
  ))
  w <- c(0, 4, 0, 4)
  errs("jtsurf", list(
    list("`metric` must be one of", x, x, w, 2, 2, "negexp", c(1, 1),
         metric = "l1"),
    list("`mode` must be one of", x, x, w, 2, 2, "negexp", c(1, 1),
         mode = "max"),
    list("0 everywhere", x, x, w, 2, 2, "linear", c(-1, 0)),
    list("at most 1e9 cells", x, x, w, 1e5, 1e5, "negexp", c(1, 1))
  ))
  errs("jtcal", list(
    list("`kern` must be one of", 1, 1, 1, kern = "gauss"),
    list("`refl` must be TRUE or FALSE", 1, 1, 1, refl = NA),
    list("`t` must be at least 0", 1, -1, 1)
  ))
  m <- matrix(1, 2, 2)
  errs("jtacc", list(
    list("must be a matrix of cell values", 1:4, w, c(1, 1)),
    list("`base` must lie inside", m, w, c(5, 1)),
    list("`prob` must be finite numbers from 0", m - 2, w, c(1, 1)),
    list("`prob` must be finite numbers from 0", matrix("a", 2, 2), w,
         c(1, 1))
  ))
  errs("bjtc", list(
    list("one row per origin zone", 1:2, matrix(1, 3, 2), 1),
    list("whole numbers from 1 to 2", 1:2, m, 3),
    list("`dest` must be at least 1", 1:2, m, 0),
    list("a block sums to 0", 1:2, matrix(c(1, 1, 0, 0), 2), 2)
  ))
})

test_that("named per-point vectors are matched by name, as elsewhere", {
  x <- c(a = 1, b = 2, c = 4, d = 7)
  y <- c(a = 1, b = 3, c = 5, d = 2)
  w <- c(a = 1, b = 2, c = 1, d = 3)
  o <- c(3, 1, 4, 2)
  r <- cgstat(x, y, w = w)
  expect_identical(cgstat(x, y[o], w = w[o]), r)
  expect_identical(sdell(x, y[o], w = w[o]), sdell(x, y, w = w))
  expect_identical(webpt(x, y[o], w = w[o]), webpt(x, y, w = w))
  expect_identical(mantst(x, y[o], t = w[o], nsim = 0),
                   mantst(x, y, t = w, nsim = 0))
  expect_identical(stmavg(x, y[o], t = w[o], span = 3),
                   stmavg(x, y, t = w, span = 3))
  expect_identical(jtcal(x, 1:3, 1, w = w[o]), jtcal(x, 1:3, 1, w = w))
  expect_identical(dirmn(x, y[o], origin = c(0, 0), w = w[o]),
                   dirmn(x, y, origin = c(0, 0), w = w))
  # unnamed or repeated names: positions
  expect_identical(cgstat(x, unname(y)), cgstat(unname(x), unname(y)))
  rep <- stats::setNames(c(1, 2, 4, 7), c("a", "a", "b", "b"))
  expect_identical(cgstat(rep, y), cgstat(unname(rep), unname(y)))
  yi <- c(a = 1L, b = 3L, c = 5L, d = 2L)
  expect_identical(cgstat(x, yi[o]), cgstat(x, yi))
  errs("cgstat", list(
    list("must name the same units as `x`", x, c(a = 1, b = 3, c = 5, e = 2)),
    list("`y` must have length 4", x, c(a = 1, b = 3, c = 5))
  ))
})
