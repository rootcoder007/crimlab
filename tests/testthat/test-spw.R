# Spatial weights, Geary's C, the global and local G, local Moran and the
# directional mean. The kernels are checked against dense re-derivations
# in src/t_swt.cpp; here the binding, and agreement with spdep (weights,
# statistics, moments; 1e-10) and circular (directions; 1e-12).

#' @srrstats {SP6.3} Neighbours by k nearest, by distance band and by
#'   contiguity are each compared with spdep's knearneigh(), dnearneigh()
#'   and a contiguity list, and the statistics are compared under each.
#' @srrstats {SP6.4} Binary, row-standardised and inverse-distance weights
#'   are compared with spdep's styles "B" and "W" and with weights written
#'   out, and Geary's C is shown to change with them.
#' @noRd
test_that("weights agree with spdep's neighbours", {
  set.seed(31)
  p <- cbind(x = runif(40), y = runif(40))
  w <- spwt(p, k = 3, style = "raw")
  expect_identical(w$n, 40L)
  expect_identical(tabulate(w$i, 40), rep(3L, 40))
  d <- as.matrix(dist(p))
  expect_equal(w$d, d[cbind(w$i, w$j)], tolerance = 1e-15)
  b <- spwt(p, type = "band", d = 0.2)
  expect_equal(as.vector(tapply(b$w, factor(b$i, 1:40), sum)),
               ifelse(tabulate(b$i, 40) > 0, 1, NA), tolerance = 1e-15)
  m <- spwt(p, type = "band", d = 0.2, metric = "manhattan", style = "raw")
  l1 <- as.matrix(dist(p, method = "manhattan"))
  expect_identical(length(m$i), as.integer(sum(l1 <= 0.2) - 40))
  iw <- spwt(p, k = 3, idw = 2, style = "raw")
  expect_equal(iw$w, iw$d^-2, tolerance = 1e-14)
  cn <- spwt(nb = list(2L, c(1L, 3L), 0L), type = "contig")
  expect_identical(cn$i, c(1L, 2L, 2L))
  expect_identical(cn$w, c(1, 0.5, 0.5))
  expect_true(all(is.na(cn$d)))
  cp <- spwt(p[1:3, ], nb = list(2, c(1, 3), 2), type = "contig", idw = 1)
  expect_equal(cp$d, d[cbind(c(1, 2, 2, 3), c(2, 1, 3, 2))], tolerance = 1e-15)
  ll <- spwt(-79 + p[, 1] / 10, 43 + p[, 2] / 10, type = "band", d = 2000,
             lonlat = TRUE)
  q <- eqproj(-79 + p[, 1] / 10, 43 + p[, 2] / 10)
  expect_identical(ll$j, spwt(q, type = "band", d = 2000)$j)
  skip_if_not_installed("spdep")
  kn <- spdep::knn2nb(spdep::knearneigh(p, k = 3))
  expect_identical(lapply(split(w$j, w$i), sort), lapply(unclass(kn), sort),
                   ignore_attr = TRUE)
  dn <- spdep::dnearneigh(p, 0, 0.2)
  bb <- split(b$j, factor(b$i, 1:40))
  for (i in 1:40) expect_setequal(bb[[i]], setdiff(dn[[i]], 0L))
})

test_that("Geary's C, global and local G, local Moran agree with spdep", {
  set.seed(32)
  p <- cbind(x = runif(50), y = runif(50))
  v <- p[, 1] + p[, 2] + rnorm(50, sd = 0.3)
  pos <- exp(2 * p[, 1]) + rexp(50)
  w <- spwt(p, k = 4)
  g <- gearc(v, w)
  expect_lt(g[["c"]], 1)
  expect_gt(gearc(v, spwt(p, k = 4, idw = 1))[["c"]], 0)
  skip_if_not_installed("spdep")
  lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(p, k = 4)))
  expect_identical(gearc(v, lw), g)
  r <- spdep::geary.test(v, lw)
  e <- r$estimate
  expect_equal(g[["c"]], unname(e[1]), tolerance = 1e-12)
  expect_equal(g[["vr"]], unname(e[3]), tolerance = 1e-10)
  expect_equal(g[["zr"]], unname(r$statistic), tolerance = 1e-10)
  expect_equal(g[["pr"]], r$p.value, tolerance = 1e-10)
  n <- spdep::geary.test(v, lw, randomisation = FALSE)
  expect_equal(g[["vn"]], unname(n$estimate[3]), tolerance = 1e-10)
  # a band wide enough that no unit is isolated: spdep's globalG.test
  # counts only units with neighbours in n(n - 1), crimlab all units (the
  # definition's sum over every pair i != j), so they differ when a unit
  # has no neighbour
  dn <- spdep::dnearneigh(p, 0, 0.3)
  expect_true(all(spdep::card(dn) > 0))
  lb <- spdep::nb2listw(dn, style = "B")
  wb <- spwt(p, type = "band", d = 0.3, style = "raw")
  o <- gorg(pos, wb)
  s <- spdep::globalG.test(pos, lb, zero.policy = TRUE)
  expect_equal(o[["g"]], unname(s$estimate[1]), tolerance = 1e-12)
  expect_equal(o[["e"]], unname(s$estimate[2]), tolerance = 1e-12)
  expect_equal(o[["v"]], unname(s$estimate[3]), tolerance = 1e-10)
  expect_equal(o[["z"]], as.vector(s$statistic), tolerance = 1e-10)
  lm <- lmorn(v, w)
  sl <- spdep::localmoran(v, lw, conditional = FALSE)
  expect_equal(lm$ii, unname(sl[, 1]), tolerance = 1e-12)
  expect_equal(lm$e, unname(sl[, 2]), tolerance = 1e-12)
  expect_equal(lm$v, unname(sl[, 3]), tolerance = 1e-10)
  expect_equal(lm$p, unname(sl[, 5]), tolerance = 1e-10)
  gs <- locgi(pos, wb, star = TRUE)
  ls <- spdep::localG(pos, spdep::nb2listw(spdep::include.self(
    dn), style = "B"))
  expect_equal(gs, as.vector(ls), tolerance = 1e-10)
  g0 <- locgi(pos, wb, star = FALSE)
  l0 <- spdep::localG(pos, lb)
  ok <- is.finite(g0)
  expect_equal(g0[ok], as.vector(l0)[ok], tolerance = 1e-10)
})

test_that("permutation tests follow the seed and find a cluster", {
  set.seed(33)
  p <- cbind(x = runif(60), y = runif(60))
  v <- ifelse(p[, 1] > 0.7 & p[, 2] > 0.7, 5, 1) + rnorm(60, sd = 0.3)
  w <- spwt(p, k = 5)
  set.seed(2)
  a <- lmorn(v, w, nsim = 99)
  set.seed(2)
  expect_identical(lmorn(v, w, nsim = 99), a)
  hot <- p[, 1] > 0.75 & p[, 2] > 0.75
  expect_true(all(a$quad[hot] == 1))
  expect_lt(median(a$pperm[hot]), 0.05)
  expect_identical(gearc(v, w, nsim = 99, seed = 1)[["pperm"]], 0.01)
  expect_identical(gorg(v, w, nsim = 99, seed = 1)[["pperm"]], 0.01)
})

test_that("directional means agree with circular", {
  a <- c(0.1, 0.4, 6.1, 1.2, 5.9)
  r <- dirmn(angle = a)
  m <- atan2(mean(sin(a)), mean(cos(a)))
  expect_equal(r[["theta"]], m %% (2 * pi), tolerance = 1e-14)
  expect_equal(r[["r"]], sqrt(mean(cos(a))^2 + mean(sin(a))^2),
               tolerance = 1e-14)
  expect_equal(r[["bearing"]], (90 - r[["theta"]] * 180 / pi) %% 360,
               tolerance = 1e-12)
  q <- dirmn(c(1, 0, -1), c(0, 1, 0), origin = c(0, 0), w = c(1, 2, 1))
  expect_equal(q[["theta"]], pi / 2, tolerance = 1e-15)
  expect_identical(q[["bearing"]], 0)
  l <- dirmn(c(-79.1, -79), c(43, 43.1), origin = c(-79.05, 43.05),
             lonlat = TRUE)
  expect_gt(l[["r"]], 0)
  skip_if_not_installed("circular")
  cc <- circular::circular(a)
  expect_equal(r[["theta"]], as.numeric(circular::mean.circular(cc)) %%
                 (2 * pi), tolerance = 1e-12)
  expect_equal(r[["r"]], circular::rho.circular(cc), tolerance = 1e-12)
  expect_equal(r[["var"]], circular::var.circular(cc), tolerance = 1e-12)
  expect_equal(r[["sd"]], circular::sd.circular(cc), tolerance = 1e-12)
})

test_that("weights and statistics refuse what they cannot use", {
  p <- cbind(x = c(0, 1, 2, 3), y = c(0, 1, 0, 1))
  errs("spwt", list(
    list("`type` must be one of", p, type = "queen"),
    list("`idw` must be one number from 0", p, idw = -1),
    list("`style` must be one of", p, style = "W"),
    list("a list of neighbour numbers", type = "contig", nb = 1:3),
    list("a list of neighbour numbers", type = "contig", nb = list("a")),
    list("one element per point", p, type = "contig", nb = list(2, 1)),
    list("need the points", type = "contig", nb = list(2, 1), idw = 1),
    list("unit numbers from 1 to 2", type = "contig", nb = list(3, 1)),
    list("must give the points", type = "knn"),
    list("at least 2 points", cbind(x = 1, y = 1)),
    list("`k` must be below", p, k = 4),
    list("`d` must be one number strictly between", p, type = "band", d = 0),
    list("two neighbours coincide", cbind(x = c(0, 0, 1), y = c(0, 0, 1)),
         k = 1, idw = 1),
    list("cannot be its own neighbour", type = "contig", nb = list(1, 1)),
    list("lists a pair of units twice", type = "contig",
         nb = list(c(2, 2), 1))
  ))
  w <- spwt(p, k = 2)
  lw <- structure(list(neighbours = list(2L, 1L, 4L, 3L),
                       weights = list(1, 1, 1, 1)), class = "listw")
  errs("gearc", list(
    list("at least 4 values", 1:3, w),
    list("the list spwt", 1:4, list(1)),
    list("is for 4 units, not 5", 1:5, w),
    list("unit numbers from 1 to 4", 1:4, list(n = 4, i = 5, j = 1, w = 1)),
    list("one row of neighbours", 1:4, structure(list(neighbours = list(),
         weights = list()), class = "listw")),
    list("integer neighbour lists", 1:4, structure(list(
      neighbours = list(2, 1, 4, 3), weights = list(1, 1, 1, 1)),
      class = "listw")),
    list("one weight per neighbour", 1:4, structure(list(
      neighbours = list(2L, 1L, 4L, 3L), weights = list(1:2, 1, 1, 1)),
      class = "listw")),
    list("a neighbour outside 1 to 4", 1:4, structure(list(
      neighbours = list(9L, 1L, 4L, 3L), weights = list(1, 1, 1, 1)),
      class = "listw")),
    list("C is undefined", rep(1, 4), lw)
  ))
  expect_identical(gearc(1:4, lw)[["s0"]], 4)
  lz <- structure(list(neighbours = list(2L, 1L, 4L, 3L, 0L),
                       weights = list(1, 1, 1, 1, NULL)), class = "listw")
  expect_identical(gearc(1:5, lz)[["s0"]], 4)
  errs("gorg", list(
    list("at least 4 values", 1:3, w),
    list("`x` must be at least 0", c(-1, 1, 1, 1), w),
    list("G is undefined", c(1, 0, 0, 0), w)
  ))
  errs("lmorn", list(
    list("at least 3 values", 1:2, w),
    list("local Moran is undefined", rep(2, 4), w)
  ))
  errs("locgi", list(
    list("at least 3 values", 1:2, w),
    list("`star` must be TRUE or FALSE", 1:4, w, star = NA)
  ))
  errs("dirmn", list(
    list("not both", c(1, 2), c(1, 2), angle = 1),
    list("not both"),
    list("lies on `origin`", c(0, 1), c(0, 1)),
    list("`w` must have a positive sum", angle = 1:2, w = c(0, 0))
  ))
})

test_that("every weights function answers one valid call in its shape", {
  p <- cbind(x = c(0, 1, 2, 3, 1), y = c(0, 1, 0, 1, 3))
  v <- c(1, 4, 2, 5, 3)
  w <- spwt(p, k = 2)
  calls <- list(
    list(quote(spwt(p, k = 2)), c("n", "i", "j", "d", "w"), "list 5"),
    list(quote(gearc(v, w)), c("c", "vn", "vr", "zn", "zr", "pn", "pr",
                               "pperm", "s0"), "numeric 9"),
    list(quote(gorg(v, w)), c("g", "e", "v", "z", "p", "pperm"),
         "numeric 6"),
    list(quote(lmorn(v, w)), c("ii", "e", "v", "z", "p", "pperm", "quad"),
         "data.frame 5"),
    list(quote(locgi(v, w)), NULL, "numeric 5"),
    list(quote(dirmn(p, origin = c(1, 1.5))),
         c("theta", "bearing", "r", "var", "sd", "w", "n"), "numeric 7")
  )
  fns <- vapply(calls, function(a) as.character(a[[1L]][[1L]]), "")
  expect_setequal(fns, spw)
  for (k in seq_along(calls)) {
    r <- eval(calls[[k]][[1L]])
    expect_named(r, calls[[k]][[2L]])
    expect_identical(paste(class(r)[1L], NROW(r)), calls[[k]][[3L]],
                     info = fns[k])
  }
})
