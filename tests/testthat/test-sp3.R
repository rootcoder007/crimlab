# Head-bang smoothing, the correlated random walk, risk-adjusted clustering,
# the Manhattan options and the zone inputs. Kernels are checked in
# src/t_sp2.cpp; here each is compared with its definition written out in
# R (independently of the C++), to 1e-12 or exactly.

l1 <- function(p) as.matrix(dist(p, method = "manhattan"))

test_that("head-bang smoothing follows its rule, written out in R", {
  set.seed(41)
  n <- 40
  p <- cbind(x = runif(n), y = runif(n))
  v <- rexp(n)
  w <- runif(n, 0.5, 2)
  hb <- function(p, v, w = NULL, nn = 6, ang = 135, ntrip = 0) {
    d <- as.matrix(dist(p))
    out <- v
    wmed <- function(a, u) {
      o <- order(a)
      a[o][which(cumsum(u[o]) >= sum(u) / 2)[1]]
    }
    for (i in seq_len(nrow(p))) {
      nb <- order(d[i, ], seq_len(nrow(p)))
      nb <- nb[nb != i][seq_len(nn)]
      tr <- NULL
      for (a in 1:(nn - 1)) for (b in (a + 1):nn) {
        u <- p[nb[a], ] - p[i, ]
        z <- p[nb[b], ] - p[i, ]
        cs <- sum(u * z) / sqrt(sum(u^2) * sum(z^2))
        if (cs <= cos(ang * pi / 180) + 1e-12) {
          tr <- rbind(tr, c(nb[a], nb[b], sqrt(sum(u^2)) + sqrt(sum(z^2))))
        }
      }
      if (is.null(tr)) next
      tr <- tr[order(tr[, 3]), , drop = FALSE]
      if (ntrip > 0) tr <- tr[seq_len(min(ntrip, nrow(tr))), , drop = FALSE]
      lo <- pmin(v[tr[, 1]], v[tr[, 2]])
      hi <- pmax(v[tr[, 1]], v[tr[, 2]])
      if (is.null(w)) {
        L <- median(lo)
        H <- median(hi)
      } else {
        wl <- ifelse(v[tr[, 1]] <= v[tr[, 2]], w[tr[, 1]], w[tr[, 2]])
        wh <- ifelse(v[tr[, 1]] <= v[tr[, 2]], w[tr[, 2]], w[tr[, 1]])
        L <- wmed(lo, wl)
        H <- wmed(hi, wh)
      }
      out[i] <- median(c(L, v[i], H))
    }
    out
  }
  expect_equal(hdbang(p, v = v, nn = 6), hb(p, v), tolerance = 1e-15)
  expect_equal(hdbang(p, v = v, nn = 6, ntrip = 3), hb(p, v, ntrip = 3),
               tolerance = 1e-15)
  expect_equal(hdbang(p, v = v, w = w, nn = 6), hb(p, v, w), tolerance = 1e-15)
  expect_equal(hdbang(p, v = v, nn = 6, iter = 2),
               hb(p, hb(p, v)), tolerance = 1e-15)
  # a spike on a lattice is removed; the rest stays
  g <- expand.grid(x = 1:6, y = 1:6)
  s <- rep(1, 36)
  s[15] <- 50
  expect_identical(hdbang(g, v = s), rep(1, 36))
  errs("hdbang", list(
    list("at least 3 points", 1:2, 1:2, v = 1:2),
    list("`nn` must be below", 1:4, 1:4, v = 1:4, nn = 4),
    list("from 90 to 180", 1:4, c(1, 3, 2, 4), v = 1:4, nn = 2, ang = 45)
  ))
})

test_that("the correlated random walk is fitted and simulated as defined", {
  set.seed(42)
  a <- cumsum(rnorm(30, 0, 0.6))
  st <- rexp(30) + 0.1
  p <- cbind(x = cumsum(st * cos(a)), y = cumsum(st * sin(a)))
  r <- crwalk(p, h = 4, nsim = 500, seed = 3)
  l <- sqrt(diff(p[, 1])^2 + diff(p[, 2])^2)
  hd <- atan2(diff(p[, 2]), diff(p[, 1]))
  tu <- (diff(hd) + pi) %% (2 * pi) - pi
  expect_equal(unname(r$fit[c("ml", "ml2", "c", "s")]),
               c(mean(l), mean(l^2), mean(cos(tu)), mean(sin(tu))),
               tolerance = 1e-12)
  k <- 1:4
  cc <- mean(cos(tu))
  ex <- k * mean(l^2) + 2 * mean(l)^2 * cc / (1 - cc) *
    (k - (1 - cc^k) / (1 - cc))
  expect_equal(r$msd$exp, ex, tolerance = 1e-12)
  ob <- vapply(k, function(j) {
    mean((p[-(1:j), 1] - p[1:(30 - j), 1])^2 +
           (p[-(1:j), 2] - p[1:(30 - j), 2])^2)
  }, 0)
  expect_equal(r$msd$obs, ob, tolerance = 1e-12)
  expect_identical(crwalk(p, h = 4, nsim = 500, seed = 3), r)
  expect_identical(dim(r$sim$x), c(500L, 4L))
  expect_equal(r$msd$x, colMeans(r$sim$x), tolerance = 1e-12)
  # longitude and latitude: the same walk in metres, centres in degrees
  q <- eqproj(-79 + p[, 1] / 1000, 43 + p[, 2] / 1000)
  rl <- crwalk(-79 + p[, 1] / 1000, 43 + p[, 2] / 1000, h = 2, nsim = 10,
               seed = 1, lonlat = TRUE)
  rq <- crwalk(q, h = 2, nsim = 10, seed = 1)
  expect_equal(rl$fit, rq$fit, tolerance = 1e-12)
  errs("crwalk", list(
    list("at least 3 locations", 1:2, 1:2),
    list("drop repeats first", c(1, 1, 2), c(1, 1, 2))
  ))
})

test_that("risk-adjusted thresholds follow the baseline", {
  set.seed(43)
  town <- cbind(x = rnorm(300), y = rnorm(300))
  inc <- rbind(town[1:50, ], cbind(x = rnorm(12, 3, 0.05),
                                   y = rnorm(12, 3, 0.05)))
  bw <- rpois(300, 3) + 1
  r <- rnnhcl(inc, base = town, bw = bw, h = 0.4, p = 0.3, minpts = 5)
  n <- nrow(inc)
  dens <- vapply(seq_len(n), function(i) {
    sum(bw * exp(-((inc[i, 1] - town[, 1])^2 + (inc[i, 2] - town[, 2])^2) /
                   (2 * 0.4^2))) / (2 * pi * 0.4^2)
  }, 0)
  lam <- n * dens / sum(bw)
  f <- 0.5 + qnorm(0.3) * sqrt((4 - pi) / (4 * pi * n))
  expect_equal(r$thr, f / sqrt(lam), tolerance = 1e-12)
  # the tight group where few people live is a cluster; the town is not
  expect_false(anyNA(r$cl[51:62]))
  expect_length(unique(r$cl[51:62]), 1L)
  # longitude and latitude: the same as on the projected points and
  # baseline (projected about the incidents' centre)
  lo <- -79 + inc[, 1] / 100
  la <- 43 + inc[, 2] / 100
  bl <- cbind(lon = -79 + town[, 1] / 100, lat = 43 + town[, 2] / 100)
  q <- eqproj(lo, la)
  a <- attr(q, "centre")
  qb <- eqproj(bl, lon0 = a[1], lat0 = a[2])
  ll <- rnnhcl(lo, la, base = bl, bw = bw, h = 400, p = 0.3, minpts = 5,
               lonlat = TRUE)
  pq <- rnnhcl(q, base = qb, bw = bw, h = 400, p = 0.3, minpts = 5)
  expect_identical(ll$cl, pq$cl)
  expect_equal(ll$thr, pq$thr, tolerance = 1e-12)
  errs("rnnhcl", list(
    list("at least 2 points", 1, 1, base = cbind(0, 0), h = 1),
    list("`k` must be 0, 1 or 2", inc, base = town, h = 1, k = 3),
    list("positive sum", inc, base = town[1:2, ], bw = c(0, 0), h = 1)
  ))
})

test_that("Manhattan distance options agree with brute force", {
  set.seed(44)
  p <- cbind(x = runif(120, 0, 2), y = runif(120))
  d <- l1(p)
  diag(d) <- Inf
  r <- nnidx(p, win = c(0, 2, 0, 1), kmax = 2, metric = "manhattan",
             nsim = 19, seed = 1, dist = TRUE)
  srt <- t(apply(d, 1, sort))
  expect_equal(r$order$mean, colMeans(srt[, 1:2]), tolerance = 1e-12)
  expect_equal(r$dist, unname(srt[, 1:2]), tolerance = 1e-15)
  lam <- 120 / 2
  expect_equal(r$order$exp, c(0.5, 0.75) * sqrt(pi / (2 * lam)),
               tolerance = 1e-14)
  expect_true(is.nan(r$ce[["donnelly"]]))
  expect_equal(r$ce[["z"]], (r$order$mean[1] - r$order$exp[1]) /
                 sqrt((4 - pi) / (8 * lam * 120)), tolerance = 1e-12)
  rr <- c(0.05, 0.1, 0.2)
  k <- ripk(p, win = c(0, 2, 0, 1), r = rr, corr = "none",
            metric = "manhattan")
  diag(d) <- Inf
  expect_equal(k$k, 2 * vapply(rr, function(s) sum(d <= s), 0) / (120 * 119),
               tolerance = 1e-12)
  expect_equal(k$l, sqrt(k$k / 2), tolerance = 1e-15)
  b <- ripk(p, win = c(0, 2, 0, 1), r = rr, corr = "border",
            metric = "manhattan")
  bd <- pmin(p[, 1], 2 - p[, 1], p[, 2], 1 - p[, 2])
  bk <- vapply(rr, function(s) {
    i <- bd >= s
    sum(d[i, ] <= s) / sum(i) / (120 / 2)
  }, 0)
  expect_equal(b$k, bk, tolerance = 1e-12)
  t <- p[, 1] + runif(120)
  ds <- as.vector(as.dist(l1(p)))
  expect_equal(mantst(p, t = t, nsim = 0, metric = "manhattan")[["r"]],
               cor(ds, as.vector(dist(t))), tolerance = 1e-12)
  dd <- l1(p)
  expect_identical(fzmode(p, r = 0.1, metric = "manhattan"),
                   as.integer(rowSums(dd <= 0.1)))
  errs("ripk", list(
    list("no Manhattan form", p, win = c(0, 2, 0, 1), r = rr,
         metric = "manhattan")
  ))
  errs("nnidx", list(
    list("`metric` must be one of", p, win = c(0, 2, 0, 1), metric = "l1")
  ))
})

test_that("zones are counts at centroids", {
  set.seed(45)
  z <- cbind(x = runif(30), y = runif(30))
  cnt <- rpois(30, 2)
  cnt[3] <- 0
  zr <- z[rep(seq_len(30), cnt), ]
  first <- match(seq_len(30), rep(seq_len(30), cnt))
  d <- as.matrix(dist(z))
  expect_equal(fzmode(z, r = 0.2, count = cnt), as.vector((d <= 0.2) %*% cnt),
               tolerance = 1e-15)
  expect_identical(fzmode(z, r = 0, count = cnt), cnt + 0)
  s <- z[c(1, 10, 20), ]
  a <- kmpart(z, k = 3, init = s, count = cnt)
  b <- kmpart(zr, k = 3, init = s)
  expect_identical(a$cl, b$cl[first])
  expect_identical(a$cen, b$cen)
  h <- nnhclu(z, area = 1, minpts = 3, count = cnt)
  g <- nnhclu(zr, area = 1, minpts = 3)
  expect_identical(h$cl, g$cl[first, , drop = FALSE])
  expect_identical(h$clusters, g$clusters)
  e <- stacel(z, win = c(0, 1, 0, 1), r = 0.15, minpts = 3, count = cnt)
  f <- stacel(zr, win = c(0, 1, 0, 1), r = 0.15, minpts = 3)
  expect_identical(e$cl, f$cl[first])
  expect_true(is.na(e$cl[3]))
  errs("fzmode", list(
    list("whole numbers, one per zone", z, r = 0.1, count = rep(0.5, 30)),
    list("`count` must have a positive sum", z, r = 0.1, count = rep(0, 30))
  ))
  errs("kmpart", list(
    list("sum to at most 1e8", z[1:2, ], k = 1, count = c(2e8, 1))
  ))
})

test_that("the new functions answer one valid call in their shape", {
  p <- cbind(x = c(0, 1, 2, 3, 1.5), y = c(0, 1, 0.2, 1, 3))
  calls <- list(
    list(quote(hdbang(p, v = 1:5, nn = 3)), NULL, "numeric 5"),
    list(quote(crwalk(p, nsim = 5, seed = 1)), c("fit", "msd", "sim"),
         "list 3"),
    list(quote(rnnhcl(p, base = p, h = 1, minpts = 2)),
         c("cl", "clusters", "thr"), "list 3")
  )
  fns <- vapply(calls, function(a) as.character(a[[1L]][[1L]]), "")
  expect_setequal(fns, sp3x)
  for (k in seq_along(calls)) {
    r <- eval(calls[[k]][[1L]])
    expect_named(r, calls[[k]][[2L]])
    expect_identical(paste(class(r)[1L], NROW(r)), calls[[k]][[3L]],
                     info = fns[k])
  }
})
