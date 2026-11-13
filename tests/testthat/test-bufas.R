# The buffer design. The kernel is checked by hand in src/t_buf.cpp; here
# the binding, against the rule written out in R.

test_that("buffer groups follow the rule", {
  set.seed(112)
  p <- cbind(x = runif(300, 0, 5000), y = runif(300, 0, 5000))
  s <- cbind(x = c(1000, 3500), y = c(1200, 3000))
  cp <- cbind(x = c(2500, 4500), y = c(500, 4000))
  b <- bufas(p, site = s, comp = cp, r = 330)
  expect_s3_class(b, "bufas")
  near <- function(q) {
    apply(p, 1, function(v) sqrt(min((q[, 1] - v[1])^2 + (q[, 2] - v[2])^2)))
  }
  ds <- near(s)
  dc <- near(cp)
  g <- ifelse(ds <= 330, "treated",
              ifelse(dc <= 330 & ds > 990, "comparison", "excluded"))
  expect_identical(as.character(b$pts$grp), g)
  expect_identical(levels(b$pts$grp), c("treated", "comparison", "excluded"))
  expect_equal(b$pts$dsite, ds, tolerance = 1e-14)
  expect_equal(b$pts$dcomp, dc, tolerance = 1e-14)
  expect_identical(unname(b$n), as.numeric(table(factor(g, levels(b$pts$grp)))))
  expect_identical(b$par[["buf"]], 990)
  expect_output(print(b), "bufas: radius 330, buffer 990")
  ll <- eqinv(p, lon0 = -79.4, lat0 = 43.7)
  sl <- eqinv(s, lon0 = -79.4, lat0 = 43.7)
  cl <- eqinv(cp, lon0 = -79.4, lat0 = 43.7)
  bl <- bufas(ll, site = sl, comp = cl, r = 330, buf = 990, lonlat = TRUE)
  # the same about the points' own centre, projected first
  rad <- pi / 180
  c0 <- atan2(mean(sin(ll[, 1] * rad)), mean(cos(ll[, 1] * rad))) / rad
  c1 <- mean(ll[, 2])
  pr <- function(v) eqproj(v, lon0 = c0, lat0 = c1)
  expect_identical(bl$pts$grp, bufas(pr(ll), site = pr(sl), comp = pr(cl),
                                     r = 330)$pts$grp)
})

test_that("bufas() refuses bad input in words", {
  p <- cbind(x = 0.5, y = 0.5)
  errs("bufas", list(
    list("`r`", p, site = p, comp = p, r = 0),
    list("`buf`", p, site = p, comp = p, r = 1, buf = -1)
  ))})
