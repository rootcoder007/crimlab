# The dynamic and heterogeneous panels: gmmpnl() against difference and
# system GMM written out here from Arellano and Bond (1991), Blundell and
# Bond (1998), Windmeijer (2005) and Roodman (2009), against the
# Anderson-Hsiao estimator where GMM is just identified, and against the
# invariance of Arellano and Bover (1995) to the transformation; wsctst()
# against the two regressions of Drukker (2003) done with lm() and a
# cluster sandwich here; ccepnl() against lm() unit by unit and pooled
# with unit-specific coefficients on the averages, and its variances
# written out from Pesaran (2006); simulated panels with known parameters;
# and the refusals.

# an AR(1) panel with unit effects, started near its stationary mean
dsim <- function(n = 30, tt = 6, a = 0.5, seed = 11) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(tt), u = seq_len(n))
  mu <- rnorm(n)
  d$y <- 0
  for (r in seq_len(nrow(d))) {
    i <- d$u[r]
    d$y[r] <- if (d$t[r] == 1) {
      mu[i] / (1 - a) + rnorm(1)
    } else {
      a * d$y[r - 1] + mu[i] + rnorm(1)
    }
  }
  d
}

# the value of y of unit i in period s, NA when there is none
yat <- function(d, i, s) {
  v <- d$y[d$u == i & d$t == s]
  if (length(v)) v else NA
}

# The rows of unit i for the reference: first differences of y on the
# lagged difference and, in system GMM, levels on the lag and a constant;
# instruments keyed by name, as Roodman (2009, secs 3.2 and 3.4) builds
# them
grows <- function(d, i, sys, a, b, coll) {
  ts <- sort(d$t[d$u == i])
  lev <- ts[!is.na(vapply(ts - 1, function(s) yat(d, i, s), 0))]
  fd <- lev[(lev - 1) %in% lev]
  out <- list()
  for (t in fd) {
    z <- numeric()
    for (l in seq(a, min(b, t - 1))) {
      v <- yat(d, i, t - l)
      if (t - l >= 1 && !is.na(v) && v != 0) {
        z[if (coll) paste0("g", l) else paste0("g", t, "_", l)] <- v
      }
    }
    x <- yat(d, i, t - 1) - yat(d, i, t - 2)
    out[[length(out) + 1]] <- list(kind = "d", t = t, z = z,
                                    x = if (sys) c(x, 0) else x,
                                    y = yat(d, i, t) - yat(d, i, t - 1))
  }
  for (t in if (sys) lev else numeric()) {
    w1 <- yat(d, i, t - a + 1)
    w0 <- yat(d, i, t - a)
    z <- c(cons = 1)
    if (!is.na(w1) && !is.na(w0) && w1 != w0) {
      z[if (coll) "lv" else paste0("lv", t)] <- w1 - w0
    }
    out[[length(out) + 1]] <- list(kind = "l", t = t, z = z,
                                    x = c(yat(d, i, t - 1), 1),
                                    y = yat(d, i, t))
  }
  out
}

# H of a unit's rows (Roodman 2009, eqs 23 and 26)
hblock <- function(r, hm) {
  k <- vapply(r, `[[`, "", "kind")
  t <- vapply(r, `[[`, 0, "t")
  n <- length(r)
  h <- diag(n)
  if (hm == "id") return(h)
  for (p in seq_len(n)) {
    for (q in seq_len(n)) {
      if (k[p] == "d" && k[q] == "d") {
        h[p, q] <- if (t[p] == t[q]) 2 else if (abs(t[p] - t[q]) == 1) -1 else 0
      } else if (hm == "ab" && k[p] != k[q]) {
        td <- if (k[p] == "d") t[p] else t[q]
        tl <- if (k[p] == "d") t[q] else t[p]
        h[p, q] <- if (tl == td) 1 else if (tl == td - 1) -1 else 0
      }
    }
  }
  h
}

# GMM with weight a, its residuals and pieces
gfit <- function(zx, zy, a) {
  mi <- solve(t(zx) %*% a %*% zx)
  list(mi = mi, p = t(zx) %*% a, b = mi %*% t(zx) %*% a %*% zy)
}

gref <- function(d, sys = FALSE, a = 2, b = Inf, coll = FALSE, hm = "ab",
                 nar = 2) {
  us <- unique(d$u)
  rs <- lapply(us, function(i) grows(d, i, sys, a, b, coll))
  keys <- unique(unlist(lapply(rs, function(r) {
    unlist(lapply(r, function(x) names(x$z)))
  })))
  zm <- function(r) {
    m <- matrix(0, length(r), length(keys), dimnames = list(NULL, keys))
    for (j in seq_along(r)) m[j, names(r[[j]]$z)] <- r[[j]]$z
    m
  }
  zs <- lapply(rs, zm)
  xs <- lapply(rs, function(r) do.call(rbind, lapply(r, `[[`, "x")))
  ys <- lapply(rs, function(r) vapply(r, `[[`, 0, "y"))
  hs <- lapply(rs, hblock, hm = hm)
  zhz <- Reduce(`+`, Map(function(z, h) t(z) %*% h %*% z, zs, hs))
  zx <- Reduce(`+`, Map(function(z, x) t(z) %*% x, zs, xs))
  zy <- Reduce(`+`, Map(function(z, y) t(z) %*% y, zs, ys))
  a1 <- solve(zhz)
  f1 <- gfit(zx, zy, a1)
  e1 <- Map(function(x, y) y - x %*% f1$b, xs, ys)
  u1 <- Map(function(z, e) t(z) %*% e, zs, e1)
  s <- Reduce(`+`, lapply(u1, function(u) u %*% t(u)))
  a2 <- solve(s)
  f2 <- gfit(zx, zy, a2)
  e2 <- Map(function(x, y) y - x %*% f2$b, xs, ys)
  v1 <- f1$mi %*% f1$p %*% s %*% t(f1$p) %*% f1$mi
  # Windmeijer (2005): D's column p
  ze2 <- zy - zx %*% f2$b
  g <- a2 %*% ze2
  dd <- vapply(seq_len(ncol(zx)), function(p) {
    w <- Reduce(`+`, Map(function(z, x, e) {
      ax <- t(z) %*% x[, p]
      be <- t(z) %*% e
      ax %*% t(be) + be %*% t(ax)
    }, zs, xs, e1))
    f2$mi %*% f2$p %*% w %*% g
  }, numeric(ncol(zx)))
  dd <- matrix(dd, ncol(zx))
  v2 <- f2$mi
  vc <- v2 + dd %*% v2 + v2 %*% t(dd) + dd %*% v1 %*% t(dd)
  ze1 <- zy - zx %*% f1$b
  dk <- lapply(rs, function(r) vapply(r, `[[`, "", "kind") == "d")
  sig <- sum(unlist(Map(function(e, k) e[k]^2, e1, dk))) /
    sum(unlist(Map(function(h, k) diag(h)[k], hs, dk)))
  hx <- NA
  lvk <- grepl("^lv", keys)
  if (sys && any(lvk)) {
    k <- !lvk
    fu <- gfit(zx[k, , drop = FALSE], zy[k], solve(s[k, k]))
    zu <- zy[k] - zx[k, , drop = FALSE] %*% fu$b
    hx <- drop(t(zu) %*% solve(s[k, k]) %*% zu)
  }
  # Arellano and Bond (1991, eqs 8-9) on the differenced residuals
  arz <- function(st, l) {
    f <- if (st == 1) f1 else f2
    e <- if (st == 1) e1 else e2
    v <- if (st == 1) v1 else vc
    num <- 0
    t1 <- 0
    q <- 0
    zv <- 0
    for (j in seq_along(rs)) {
      r <- rs[[j]]
      tt <- vapply(r, `[[`, 0, "t")
      w <- vapply(seq_along(r), function(h) {
        if (!dk[[j]][h]) return(0)
        o <- which(dk[[j]] & tt == tt[h] - l)
        if (length(o)) e[[j]][o] else 0
      }, 0)
      c0 <- sum(w * e[[j]] * dk[[j]])
      num <- num + c0
      t1 <- t1 + c0^2
      q <- q + colSums(w * dk[[j]] * xs[[j]])
      zv <- zv + t(zs[[j]]) %*% e[[j]] * c0
    }
    vv <- t1 - 2 * drop(t(q) %*% f$mi %*% f$p %*% zv) +
      drop(t(q) %*% v %*% q)
    num / sqrt(vv)
  }
  list(b1 = drop(f1$b), b2 = drop(f2$b), v1 = v1, v2 = v2, vc = vc,
       sargan = drop(t(ze1) %*% a1 %*% ze1) / sig,
       hansen = drop(t(ze2) %*% a2 %*% ze2), hx = hx, m = length(keys),
       ar = function(st) vapply(seq_len(nar), function(l) arz(st, l), 0))
}

tst <- function(r, nm) r$tests$stat[r$tests$test == nm]

test_that("difference GMM is the estimator written out, both steps", {
  d <- dsim()
  g <- gref(d)
  r1 <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE)
  expect_identical(r1$coef$term, "L1.y")
  expect_equal(r1$fit$inst, g$m)
  expect_equal(r1$coef$est, g$b1, tolerance = 1e-9)
  expect_equal(r1$coef$se, sqrt(diag(g$v1)), tolerance = 1e-8)
  expect_equal(tst(r1, "Sargan"), g$sargan, tolerance = 1e-8)
  expect_equal(tst(r1, "Hansen"), g$hansen, tolerance = 1e-8)
  expect_equal(r1$tests$stat[1:2], g$ar(1), tolerance = 1e-8)
  expect_equal(r1$tests$p[1:2], 2 * pnorm(-abs(g$ar(1))), tolerance = 1e-8)
  expect_equal(r1$tests$p[4], pchisq(g$hansen, g$m - 1, lower.tail = FALSE),
               tolerance = 1e-8)
  r2 <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE,
               steps = 2)
  expect_equal(r2$coef$est, g$b2, tolerance = 1e-9)
  expect_equal(unname(r2$vcov), unname(g$vc), tolerance = 1e-8)
  expect_equal(r2$tests$stat[1:2], g$ar(2), tolerance = 1e-8)
  r3 <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE,
               steps = 2, wc = FALSE)
  expect_equal(unname(r3$vcov), unname(g$v2), tolerance = 1e-8)
  expect_identical(r1$fit$n, 120)
  expect_identical(r1$fit$groups, 30)
  expect_identical(r1$fit$omin, 4)
})

test_that("system GMM with a gap, collapsed and limited lags, every H", {
  d <- dsim(seed = 12)
  d <- d[-c(9, 40), ]
  for (hm in c("ab", "dpd", "id")) {
    g <- gref(d, sys = TRUE, b = 3, coll = TRUE, hm = hm)
    r <- gmmpnl(d$y, unit = d$u, time = d$t, tdum = FALSE, lagmax = 3,
                collapse = TRUE, hmat = hm, steps = 2)
    expect_identical(r$coef$term, c("L1.y", "(const)"))
    expect_equal(r$fit$inst, g$m)
    expect_equal(r$coef$est, g$b2, tolerance = 1e-9, info = hm)
    expect_equal(unname(r$vcov), unname(g$vc), tolerance = 1e-7, info = hm)
    expect_equal(tst(r, "Sargan"), g$sargan, tolerance = 1e-8, info = hm)
    expect_equal(tst(r, "Hansen"), g$hansen, tolerance = 1e-8, info = hm)
    expect_equal(tst(r, "Hansen without level GMM"), g$hx, tolerance = 1e-8)
    expect_equal(tst(r, "Difference-in-Hansen"), g$hansen - g$hx,
                 tolerance = 1e-7)
    expect_equal(r$tests$stat[1:2], g$ar(2), tolerance = 1e-7, info = hm)
  }
  r1 <- gmmpnl(d$y, unit = d$u, time = d$t, tdum = FALSE, lagmax = 3,
               collapse = TRUE)
  g1 <- gref(d, sys = TRUE, b = 3, coll = TRUE)
  expect_equal(r1$coef$est, g1$b1, tolerance = 1e-9)
  expect_equal(unname(r1$vcov), unname(g1$v1), tolerance = 1e-8)
  # system GMM counts the rows in levels: 180 less the two gaps and the
  # first period of each unit, and the period after each gap
  expect_identical(r1$fit$n, 180 - 2 - 30 - 2)
  expect_identical(r1$tests$df[r1$tests$test == "Difference-in-Hansen"], 1)
})

test_that("uncollapsed system GMM is the estimator written out", {
  d <- dsim(seed = 13)
  g <- gref(d, sys = TRUE)
  r <- gmmpnl(d$y, unit = d$u, time = d$t, tdum = FALSE)
  expect_equal(r$fit$inst, g$m)
  expect_equal(r$coef$est, g$b1, tolerance = 1e-9)
  expect_equal(r$tests$stat[1:2], g$ar(1), tolerance = 1e-7)
})

test_that("with T = 3 difference GMM is Anderson and Hsiao's estimator", {
  # one equation, in period 3, instrumented by y_1: just identified, so
  # both steps give sum y1 dy3 / sum y1 dy2 and J is 0
  d <- dsim(n = 25, tt = 3, seed = 14)
  w <- reshape(d, idvar = "u", timevar = "t", direction = "wide")
  ah <- sum(w$y.1 * (w$y.3 - w$y.2)) / sum(w$y.1 * (w$y.2 - w$y.1))
  for (s in 1:2) {
    r <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE,
                steps = s)
    expect_equal(r$coef$est, ah, tolerance = 1e-12)
    expect_equal(tst(r, "Hansen"), 0, tolerance = 1e-10)
  }
  expect_identical(r$fit$inst, 1)
  expect_true(all(is.na(r$tests$stat[1:2])))
})

test_that("first differences and deviations agree on a balanced panel", {
  # Arellano and Bover (1995, sec. 3): with all lags as instruments the
  # estimator is invariant to the transformation, and so are the moments
  # behind the tests
  d <- dsim(seed = 15)
  for (s in 1:2) {
    a <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE,
                steps = s)
    b <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE,
                steps = s, trans = "fod")
    expect_identical(b$trans, "fod")
    expect_equal(b$coef$est, a$coef$est, tolerance = 1e-9)
    expect_equal(b$coef$se, a$coef$se, tolerance = 1e-8)
    # Sargan's variance estimate is not invariant: the other tests are
    k <- a$tests$test != "Sargan"
    expect_equal(b$tests$stat[k], a$tests$stat[k], tolerance = 1e-7)
  }
})

test_that("system GMM recovers the persistence of a simulated panel", {
  # 300 units, T = 6, a = 0.5 in the data: the estimate should lie within
  # four standard errors (a normal draw beyond that has chance 6e-5)
  d <- dsim(n = 300, a = 0.5, seed = 16)
  r <- gmmpnl(d$y, unit = d$u, time = d$t, lagmax = 3, collapse = TRUE)
  expect_lt(abs(r$coef$est[1] - 0.5), 4 * r$coef$se[1])
  r2 <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", steps = 2,
               trans = "fod", lagmax = 4)
  expect_lt(abs(r2$coef$est[1] - 0.5), 4 * r2$coef$se[1])
})

test_that("regressors, instruments, dummies and the note", {
  d <- dsim(n = 12, seed = 17)
  set.seed(18)
  d$x <- rnorm(nrow(d))
  d$w <- rnorm(nrow(d))
  d$zz <- rnorm(nrow(d))
  r <- gmmpnl(d$y, d[c("x", "w")], d$u, d$t, gmm = "x", z = d["zz"],
              collapse = TRUE, lagmin = c(2, 1), lagmax = c(3, 2))
  # time dummies: period 1 has no level row, and the last dummy is the
  # constant less the others
  expect_identical(r$coef$term,
                   c("L1.y", "x", "w", "(const)", paste0("t", 2:5)))
  expect_identical(r$dropped, c("t1", "t6"))
  expect_identical(r$dummies, paste0("t", 2:5))
  # y: lags 2-3 and one level column; x: lags 1-2 and one level column;
  # w, zz, the constant and four dummies IV-style
  expect_identical(r$fit$inst, 13)
  expect_identical(r$inst[1:3], c("L2.y", "L3.y", "D.L1.y"))
  expect_identical(r$inst[4:6], c("L1.x", "L2.x", "D.L0.x"))
  expect_true(all(c("w", "zz", "(const)", "t3") %in% r$inst))
  expect_identical(r$note[1], paste("13 instruments for 12 groups: more",
                                 "instruments than groups weakens the",
                                 "Hansen test (Roodman 2009, sec. 2.6)."))
  expect_output(print(r), "note: 13 instruments")
  # iv = character(0) leaves w a regressor without an instrument of its
  # own, and ar = 0 drops the Arellano-Bond tests
  q <- gmmpnl(d$y, d[c("x", "w")], d$u, d$t, gmm = "x", iv = character(0),
              collapse = TRUE, lagmax = 3, ar = 0, type = "diff")
  expect_false("w" %in% q$inst)
  expect_identical(q$tests$test, c("Sargan", "Hansen"))
  u <- gmmpnl(d$y, d["x"], d$u, d$t, iv = "x", ylag = 2, collapse = TRUE,
              type = "diff", tdum = FALSE)
  expect_identical(u$coef$term, c("L1.y", "L2.y", "x"))
  expect_identical(u$inst, c("L2.y", "L3.y", "L4.y", "L5.y", "x"))
  n <- gmmpnl(d$y, d["x"], d$u, d$t, ylag = 0, z = d["zz"], tdum = FALSE,
              type = "diff")
  expect_identical(n$coef$term, "x")
  expect_identical(n$inst, c("x", "zz"))
  expect_identical(n$tests$df[n$tests$test == "Hansen"], 1)
})

test_that("a singular weight matrix is inverted generally and said so", {
  # four groups and fifteen instruments: the two-step moment variance has
  # rank four at most
  d <- dsim(n = 4, tt = 7, seed = 19)
  r <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", tdum = FALSE,
              steps = 2)
  expect_identical(r$fit$inst, 15)
  expect_length(r$note, 2)
  expect_match(r$note[2], "generalized inverse")
})

test_that("gmmpnl() refusals name the function", {
  d <- dsim(n = 6, tt = 4, seed = 20)
  d$x <- d$t
  d$k <- rep(1:6, each = 4)
  y <- d$y
  u <- d$u
  t <- d$t
  expect_error(gmmpnl(y, unit = u, time = t, ylag = 0), "gmmpnl: .*ylag")
  expect_error(gmmpnl(y, unit = u, time = t, ylag = 1001), "at most 1000")
  expect_error(gmmpnl(y, d["x"], u, t, gmm = "q"), "not a column")
  expect_error(gmmpnl(y, d["x"], u, t, gmm = c("x", "x")), "twice")
  expect_error(gmmpnl(y, d["x"], u, t, gmm = NA_character_), "NA")
  expect_error(gmmpnl(y, d["x"], u, t, gmm = 1), "names of columns")
  expect_error(gmmpnl(y, d["x"], u, t, gmm = "x", iv = "x"), "both name")
  expect_error(gmmpnl(y, unit = u, time = t, lagmin = c(2, 3)), "one per")
  expect_error(gmmpnl(y, unit = u, time = t, lagmin = 1.5), "whole numbers")
  expect_error(gmmpnl(y, unit = u, time = t, lagmin = 3, lagmax = 2),
               "at least `lagmin`")
  expect_error(gmmpnl(y, unit = u, time = t, steps = 3), "1 or 2")
  expect_error(gmmpnl(y, unit = u, time = t, ar = 101), "at most 100")
  expect_error(gmmpnl(y, unit = u, time = t, type = "x"), "`type`")
  expect_error(gmmpnl(y, unit = u, time = t, trans = "x"), "`trans`")
  expect_error(gmmpnl(y, unit = u, time = t, hmat = "x"), "`hmat`")
  expect_error(gmmpnl(y, unit = u, time = t, ylag = 3, type = "diff"),
               "no unit has")
  expect_error(gmmpnl(y, d["k"], u, t, type = "diff"), "regressor k")
  expect_error(gmmpnl(y, d["x"], u, t, iv = character(0), collapse = TRUE,
                      lagmax = 2, type = "diff", tdum = FALSE),
               "fewer than the 2")
  # sum y1 (y2 - y1) = 0: the one instrument identifies nothing
  e <- data.frame(u = rep(1:2, each = 3), t = rep(1:3, 2),
                  y = c(1, 2, 0, 1, 0, 5))
  expect_error(gmmpnl(e$y, unit = e$u, time = e$t, type = "diff",
                      tdum = FALSE), "do not identify")
  # 70 periods give some 2400 instruments uncollapsed
  b <- dsim(n = 3, tt = 70, seed = 21)
  expect_error(gmmpnl(b$y, unit = b$u, time = b$t, type = "diff",
                      tdum = FALSE), "too large \\(2346")
  # sizes refused before they are made
  w <- data.frame(u = rep(1:2, each = 2001), t = rep(1:2001, 2),
                  y = sin(1:4002))
  expect_error(gmmpnl(w$y, unit = w$u, time = w$t), "regressors would")
  f <- data.frame(u = rep(1:2, each = 11600), t = rep(1:11600, 2),
                  y = sin(1:23200))
  expect_error(gmmpnl(f$y, unit = f$u, time = f$t, trans = "fod",
                      type = "diff", tdum = FALSE), "too large")
  g <- data.frame(u = rep(1:2, each = 8200), t = rep(1:8200, 2),
                  y = sin(1:16400))
  expect_error(gmmpnl(g$y, unit = g$u, time = g$t, type = "diff",
                      tdum = FALSE), "too large")
  h <- data.frame(u = rep(1:2, each = 6000), t = rep(1:6000, 2),
                  y = sin(1:12000))
  expect_error(gmmpnl(h$y, unit = h$u, time = h$t, tdum = FALSE,
                      collapse = TRUE, lagmax = 2), "too large")
})

test_that("print and plot of gmmpnl() return their input and what they drew", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  d <- dsim(seed = 22)
  r <- gmmpnl(d$y, unit = d$u, time = d$t, collapse = TRUE)
  expect_output(x <- print(r), "gmmpnl: system GMM \\(fd\\), step 1")
  expect_identical(x, r)
  p <- plot(r)
  expect_identical(p$label, "L1.y")
  expect_identical(p$est, r$coef$est[1])
  q <- gmmpnl(d$y, unit = d$u, time = d$t, type = "diff", collapse = TRUE)
  expect_output(print(q), "difference GMM")
  expect_error(plot(structure(list(), class = "gmmpnl")), "not a complete")
})

# CR1 of a regression without weights, written out
cr1 <- function(x, e, cl) {
  b <- solve(crossprod(x))
  s <- rowsum(x * e, cl)
  g <- nrow(s)
  n <- nrow(x)
  b %*% crossprod(s) %*% b * g / (g - 1) * (n - 1) / (n - ncol(x))
}

wsim <- function(seed = 23, r = 0) {
  set.seed(seed)
  d <- expand.grid(t = 1:7, u = 1:40)
  d$x <- rnorm(280)
  d$e <- 0
  for (i in seq_len(280)) {
    d$e[i] <- if (d$t[i] == 1) rnorm(1) else r * d$e[i - 1] + rnorm(1)
  }
  d$y <- 2 * d$x + rnorm(40)[d$u] + d$e
  d[-c(10, 100), ]
}

test_that("wsctst() is Drukker's two regressions", {
  d <- wsim()
  r <- wsctst(d$y, d["x"], d$u, d$t)
  # first differences where the period before is there
  key <- paste(d$u, d$t)
  pre <- match(paste(d$u, d$t - 1), key)
  k <- !is.na(pre)
  f <- data.frame(u = d$u[k], t = d$t[k], dy = d$y[k] - d$y[pre[k]],
                  dx = d$x[k] - d$x[pre[k]])
  m <- lm(dy ~ dx - 1, f)
  expect_equal(r$coef$est, unname(coef(m)), tolerance = 1e-12)
  expect_equal(r$coef$se, sqrt(drop(cr1(cbind(f$dx), residuals(m), f$u))),
               tolerance = 1e-10)
  f$e <- residuals(m)
  fk <- paste(f$u, f$t)
  lag <- match(paste(f$u, f$t - 1), fk)
  j <- !is.na(lag)
  e1 <- f$e[j]
  e0 <- f$e[lag[j]]
  m2 <- lm(e1 ~ e0 - 1)
  se <- sqrt(drop(cr1(cbind(e0), residuals(m2), f$u[j])))
  g <- length(unique(f$u[j]))
  fs <- ((coef(m2)[[1]] + 0.5) / se)^2
  s <- r$stat
  expect_equal(s[["coef"]], coef(m2)[[1]], tolerance = 1e-12)
  expect_equal(s[["se"]], se, tolerance = 1e-10)
  expect_equal(s[["F"]], fs, tolerance = 1e-10)
  expect_equal(s[["p"]], pf(fs, 1, g - 1, lower.tail = FALSE),
               tolerance = 1e-10)
  expect_equal(s[["df2"]], g - 1)
  expect_equal(s[["pairs"]], sum(j))
  expect_equal(s[["ndiff"]], nrow(f))
  # with a constant, and without the small-sample factor
  m3 <- lm(e1 ~ e0)
  c3 <- wsctst(d$y, d["x"], d$u, d$t, cons = TRUE)
  expect_equal(c3$stat[["coef"]], coef(m3)[[2]], tolerance = 1e-12)
  n <- wsctst(d$y, d["x"], d$u, d$t, ssc = "none")
  fac <- g / (g - 1) * (sum(j) - 1) / (sum(j) - 1)
  expect_equal(n$stat[["se"]]^2 * fac, se^2, tolerance = 1e-12)
  expect_identical(n$ssc, "none")
  expect_null(wsctst(d$y, unit = d$u, time = d$t)$coef)
})

test_that("wsctst() has the right size and power on simulated panels", {
  # with no correlation the lag coefficient lies within four standard
  # errors of -0.5 (a normal draw beyond that has chance 6e-5); with AR(1)
  # errors of 0.6 the coefficient is near -0.2 and, with 40 units, the
  # test rejects at 1% (seed fixed)
  w0 <- wsim(24)
  r0 <- wsctst(w0$y, w0["x"], w0$u, w0$t)
  expect_lt(abs(r0$stat[["coef"]] + 0.5), 4 * r0$stat[["se"]])
  d <- wsim(25, r = 0.6)
  r1 <- wsctst(d$y, d["x"], d$u, d$t)
  expect_lt(r1$stat[["p"]], 0.01)
  expect_output(print(r1), "wsctst: Wooldridge test")
})

test_that("wsctst() refusals name the function", {
  d <- data.frame(u = c(1, 1, 1, 2, 2), t = c(1, 2, 3, 1, 3),
                  y = c(1, 3, 2, 4, 5))
  expect_error(wsctst(d$y, unit = d$u, time = d$t), "wsctst: too few")
  expect_error(wsctst(d$y, unit = d$u, time = d$t, ssc = "x"), "`ssc`")
  expect_error(wsctst(d$y, unit = d$u, time = d$t, cons = NA), "`cons`")
  expect_error(wsctst(d$y, matrix(0, 5, 2001), d$u, d$t), "at most 2000")
})

csim <- function(seed = 26, n = 20, tt = 10, fac = TRUE) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(tt), u = seq_len(n))
  f <- if (fac) rnorm(tt)[d$t] else 0
  g <- rnorm(n, 1)[d$u]
  d$x <- g * f + rnorm(n * tt)
  d$b <- (0.7 + rnorm(n, sd = 0.1))[d$u]
  d$y <- rnorm(n)[d$u] + d$b * d$x + g * f + rnorm(n * tt)
  d
}

test_that("ccepnl() is lm() on the averages, unit by unit and pooled", {
  d <- csim()
  d$ybar <- ave(d$y, d$t)
  d$xbar <- ave(d$x, d$t)
  r <- ccepnl(d$y, d["x"], d$u, d$t)
  bi <- vapply(1:20, function(i) {
    coef(lm(y ~ x + ybar + xbar, d[d$u == i, ]))[["x"]]
  }, 0)
  expect_equal(unname(r$bi[, "x"]), bi, tolerance = 1e-10)
  expect_equal(r$mg$est, mean(bi), tolerance = 1e-12)
  expect_equal(r$mg$se, sd(bi) / sqrt(20), tolerance = 1e-10)
  p <- lm(y ~ x + factor(u) + factor(u):ybar + factor(u):xbar, d)
  expect_equal(r$pooled$est, coef(p)[["x"]], tolerance = 1e-10)
  # eqs 67-69 with w_i = 1 / N: var = Psi^-1 R Psi^-1 / N
  a <- vapply(1:20, function(i) {
    s <- d[d$u == i, ]
    sum(residuals(lm(x ~ ybar + xbar, s))^2) / 10
  }, 0)
  ps <- mean(a)
  rr <- sum(a^2 * (bi - mean(bi))^2) / 19
  expect_equal(r$pooled$se, sqrt(rr / ps^2 / 20), tolerance = 1e-10)
  expect_equal(r$pooled$lo, r$pooled$est - qnorm(0.975) * r$pooled$se,
               tolerance = 1e-14)
  expect_identical(r$fit$h, 3)
  expect_equal(unname(r$w), rep(0.05, 20), tolerance = 1e-15)
})

test_that("ccepnl() weights and observed factors enter as Pesaran sets out", {
  d <- csim(seed = 27)
  w <- (1 + d$u %% 3) / 40
  r <- ccepnl(d$y, d["x"], d$u, d$t, w = w)
  wu <- (1 + (1:20) %% 3)
  wu <- wu / sum(wu)
  d$ybar <- ave(d$y * wu[d$u], d$t, FUN = sum)
  d$xbar <- ave(d$x * wu[d$u], d$t, FUN = sum)
  d$tr <- d$t
  bi <- vapply(1:20, function(i) {
    coef(lm(y ~ x + ybar + xbar, d[d$u == i, ]))[["x"]]
  }, 0)
  expect_equal(unname(r$bi[, 1]), bi, tolerance = 1e-10)
  expect_equal(unname(r$w), wu, tolerance = 1e-14)
  a <- vapply(1:20, function(i) {
    sum(residuals(lm(x ~ ybar + xbar, d[d$u == i, ]))^2) / 10
  }, 0)
  ps <- sum(wu * a)
  bp <- sum(wu * vapply(1:20, function(i) {
    s <- d[d$u == i, ]
    sum(residuals(lm(x ~ ybar + xbar, s)) *
          residuals(lm(y ~ ybar + xbar, s))) / 10
  }, 0)) / ps
  expect_equal(r$pooled$est, bp, tolerance = 1e-10)
  wt <- wu / sqrt(mean(wu^2))
  rr <- sum(wt^2 * a^2 * (bi - mean(bi))^2) / 19
  expect_equal(r$pooled$se, sqrt(sum(wu^2) * rr / ps^2), tolerance = 1e-10)
  # a trend as an observed factor, and no constant
  f <- ccepnl(d$y, d["x"], d$u, d$t, fac = d["tr"], cons = FALSE)
  bt <- vapply(1:20, function(i) {
    s <- d[d$u == i, ]
    s$yb <- ave(d$y, d$t)[d$u == i]
    s$xb <- ave(d$x, d$t)[d$u == i]
    coef(lm(y ~ x + tr + yb + xb - 1, s))[["x"]]
  }, 0)
  expect_equal(unname(f$bi[, 1]), bt, tolerance = 1e-10)
})

test_that("CCE without factors comes near the mean slope", {
  # no common factor: the unit slopes are 0.7 plus N(0, 0.1^2) draws, so
  # CCEMG and CCEP lie within four of their standard errors of 0.7
  d <- csim(seed = 28, n = 60, tt = 15, fac = FALSE)
  r <- ccepnl(d$y, d["x"], d$u, d$t)
  expect_lt(abs(r$mg$est - 0.7), 4 * r$mg$se)
  expect_lt(abs(r$pooled$est - 0.7), 4 * r$pooled$se)
  # with a factor in y and x, CCEMG stays near 0.7
  e <- csim(seed = 29, n = 60, tt = 15)
  q <- ccepnl(e$y, e["x"], e$u, e$t)
  expect_lt(abs(q$mg$est - 0.7), 4 * q$mg$se)
})

test_that("ccepnl() refusals, print and plot", {
  d <- csim(n = 5, tt = 6)
  expect_error(ccepnl(d$y[-1], d["x"][-1, , drop = FALSE], d$u[-1], d$t[-1]),
               "ccepnl: the panel must be balanced")
  expect_error(ccepnl(d$y, d["x"], d$u, d$t, fac = seq_len(30)),
               "constant within each period")
  expect_error(ccepnl(d$y, d["x"], d$u, d$t, w = seq_len(30)),
               "constant within each unit")
  expect_error(ccepnl(d$y, d["x"], d$u, d$t, w = rep(0, 30)), "`w`")
  expect_error(ccepnl(d$y, cbind(d$x, d$x^2, d$t), d$u, d$t), "periods must")
  expect_error(ccepnl(d$y, ave(d$x, d$t), d$u, d$t), "collinear")
  expect_error(ccepnl(d$y, matrix(0, 30, 2001), d$u, d$t), "at most 2000")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  r <- ccepnl(d$y, d["x"], d$u, d$t)
  expect_output(x <- print(r), "pooled \\(CCEP\\)")
  expect_identical(x, r)
  p <- plot(r)
  expect_identical(p$group, c("CCEMG", "CCEP"))
  expect_identical(p$est, c(r$mg$est, r$pooled$est))
  expect_error(plot(structure(list(), class = "ccepnl")), "not a complete")
})
