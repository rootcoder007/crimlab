# The dossier functions. Their kernels are checked in C++ (src/t_fair.cpp,
# t_tlv.cpp, t_aor.cpp, t_bnf.cpp, t_sk.cpp, t_rpt.cpp, t_sens.cpp,
# t_trj.cpp, t_mvt.cpp, t_tfr.cpp, t_lws.cpp, t_chm.cpp, t_gst.cpp,
# t_spt.cpp); here the binding, against base R, stats and spdep where they
# compute the same thing, and hand-computed cases. Tolerances are stated
# per comparison.

test_that("fairness metrics agree with tables and prop.test", {
  set.seed(21)
  g <- sample(c("b", "a"), 300, TRUE)
  y <- rbinom(300, 1, ifelse(g == "a", 0.3, 0.5))
  s <- runif(300)
  r <- fair(s, y == 1, factor(g), cut = 0.4, base = "b")
  expect_identical(r$cnt$group, c("a", "b"))
  for (k in c("a", "b")) {
    i <- g == k
    p <- s[i] > 0.4
    tp <- sum(p & y[i] == 1)
    fp <- sum(p & y[i] == 0)
    tn <- sum(!p & y[i] == 0)
    fn <- sum(!p & y[i] == 1)
    j <- match(k, r$cnt$group)
    expect_identical(unlist(r$cnt[j, -1], use.names = FALSE),
                     c(tp, fp, tn, fn, sum(i)) + 0)
    expect_equal(r$met$tpr[j], tp / (tp + fn), tolerance = 1e-15)
    expect_equal(r$met$ppv[j], tp / (tp + fp), tolerance = 1e-15)
    mcc <- (tp * tn - fp * fn) /
      sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    expect_equal(r$met$mcc[j], mcc, tolerance = 1e-14)
    w <- suppressWarnings(prop.test(tp, tp + fn, correct = FALSE)$conf.int)
    expect_equal(c(r$lo$tpr[j], r$hi$tpr[j]), as.vector(w), tolerance = 1e-12)
  }
  expect_equal(r$rat$fpr, r$met$fpr / r$met$fpr[2], tolerance = 1e-15)
  expect_identical(fair(s, y, g)$cnt$group, c("a", "b"))
})

test_that("the two-level model agrees with its formulas", {
  set.seed(22)
  m <- 8
  th <- matrix(rnorm(2 * m, 5, 1), m, 2)
  it <- rep(seq_len(m), c(2, 3, 2, 4, 2, 3, 2, 2))
  x <- th[it, ] + matrix(rnorm(2 * length(it), 0, 0.3), ncol = 2)
  b <- tlvc(x, it)
  n <- nrow(x)
  xb <- rowsum(x, it) / as.vector(table(it))
  r <- x - xb[it, ]
  u <- crossprod(r) / (n - m)
  ni <- as.vector(table(it))
  d <- sweep(xb, 2, colMeans(x))
  sb <- crossprod(d * sqrt(ni))
  cc <- (sb - (m - 1) * u) / (n - sum(ni^2) / n)
  expect_equal(b$mu, colMeans(x), tolerance = 1e-14)
  expect_equal(b$u, u, tolerance = 1e-13, ignore_attr = TRUE)
  expect_equal(b$c, cc, tolerance = 1e-12, ignore_attr = TRUE)
  ldn <- function(v, mu, s) {
    l <- chol(s)
    z <- backsolve(l, v - mu, transpose = TRUE)
    -sum(log(diag(l))) - 0.5 * sum(z^2) - length(v) / 2 * log(2 * pi)
  }
  bg <- list(mu = colMeans(x), u = diag(c(0.09, 0.08)),
             c = matrix(c(1, 0.2, 0.2, 1.1), 2))
  xa <- c(5, 5.2)
  xc <- c(5.1, 5)
  ybar <- (2 * xa + 3 * xc) / 5
  ln <- ldn(xa - xc, 0, bg$u * (1 / 2 + 1 / 3)) +
    ldn(ybar, bg$mu, bg$u / 5 + bg$c) - ldn(xa, bg$mu, bg$u / 2 + bg$c) -
    ldn(xc, bg$mu, bg$u / 3 + bg$c)
  v <- tlvlr(xa, xc, 2, 3, bg)
  expect_equal(v[["ln"]], ln, tolerance = 1e-12)
  expect_equal(v[["log10"]], ln / log(10), tolerance = 1e-12)
})

test_that("aoristic weights, Benford digits and Scott-Knott groups", {
  f <- as.POSIXct(c("2024-03-04 22:00", "2024-03-05 09:30", NA),
                  tz = "UTC")
  t <- as.POSIXct(c("2024-03-05 02:00", NA, "2024-03-05 10:00"), tz = "UTC")
  a <- aorwt(f, t, bw = 3600, nb = 24, origin = 0, each = TRUE)
  expect_equal(a$tot[c(23, 24, 1, 2)], rep(0.25, 4), tolerance = 1e-15)
  expect_identical(a$tot[10], 1)
  expect_identical(unname(a$cnt), c(1, 1, 0, 0, 0))
  expect_identical(dim(a$w), c(3L, 24L))
  # one hour later on a clock an hour ahead from 2024-03-05 00:00
  ch <- as.numeric(as.POSIXct("2024-03-05", tz = "UTC"))
  b <- aorwt(f[1], t[1], bw = 3600, nb = 24, origin = 0, tt = ch,
             toff = 3600)
  expect_equal(b$tot[c(23, 24, 2, 3)], rep(0.25, 4), tolerance = 1e-15)
  d <- bdig(c(0.29, 290, 0, NA, -5, 1e-30), nd = 2)
  expect_identical(d$dig, c(29L, 29L, 0L, 0L, 0L, 0L))
  expect_identical(unname(d$cnt), c(1, 1, 1, 1))
  expect_identical(bdig(c(-5, 7), sgn = "neg")$dig, c(5L, 0L))
  expect_identical(bdig(c(-5, 7), sgn = "both", dp = 0)$dig, c(5L, 7L))
  set.seed(23)
  x <- exp(runif(400, 0, 8))
  dg <- bdig(x)$dig
  bt <- bnft(dg, nsim = 99, seed = 2)
  ben <- log10(1 + 1 / 1:9)
  ob <- tabulate(dg, 9)
  ch <- suppressWarnings(chisq.test(ob, p = ben))
  expect_equal(bt$stat[["chisq"]], unname(ch$statistic), tolerance = 1e-12)
  expect_equal(bt$stat[["p"]], ch$p.value, tolerance = 1e-10)
  expect_equal(bt$stat[["mad"]], mean(abs(ob / 400 - ben)), tolerance = 1e-14)
  expect_equal(bt$obs$exp, 400 * ben, tolerance = 1e-14)
  expect_true(is.nan(bnft(dg)$stat[["pmc"]]))
  set.seed(24)
  v <- c(rnorm(6, 1.5180, 4e-5), rnorm(6, 1.5195, 4e-5))
  expect_identical(skp(v, 4e-5^2, nsim = 199, seed = 1), rep(1:2, each = 6))
  expect_identical(skp(v[1:6], 4e-5^2, nsim = 199, seed = 1), rep(1L, 6))
})

test_that("constant reporters, sensitivity bounds and lowess", {
  r <- rptc(c("b", "a", "a", "b", "b", "c"), c(1, 1, 2, 1, 2, 1),
            c(12, 12, 12, 3, 12, 12), 1, 2)
  expect_identical(r$unit, c("a", "b", "c"))
  expect_identical(r$const, c(TRUE, TRUE, FALSE))
  expect_identical(rptc(c("b", "b"), c(1, 2), c(3, 3), 1, 2,
                        minrep = 4)$const, FALSE)
  set.seed(25)
  d <- rnorm(40, 0.8)
  w <- rswp(d, gam = c(1, 2))
  p1 <- wilcox.test(d, alternative = "greater", exact = FALSE,
                    correct = FALSE)$p.value
  expect_equal(w$p$pup[1], p1, tolerance = 1e-12)
  expect_equal(w$p$plo[1], p1, tolerance = 1e-12)
  expect_gt(w$p$pup[2], w$p$plo[2])
  expect_gt(w$gcrit, 1)
  h <- rshl(d, gam = c(1, 3))
  hl <- wilcox.test(d, conf.int = TRUE)$estimate
  expect_equal(h$lo[1], unname(hl), tolerance = 1e-12)
  expect_equal(h$hi[1], unname(hl), tolerance = 1e-12)
  expect_lt(h$lo[2], h$hi[2])
  x <- runif(80)
  y <- sin(6 * x) + rnorm(80, sd = 0.2)
  l <- stats::lowess(x, y)
  expect_equal(lwss(x, y)[order(x)], l$y, tolerance = 1e-12)
  expect_equal(lwss(x, y, f = 0.3, iter = 0, delta = 0)[order(x)],
               stats::lowess(x, y, f = 0.3, iter = 0, delta = 0)$y,
               tolerance = 1e-12)
})

test_that("trajectory models agree with the likelihood written out", {
  set.seed(26)
  x <- cbind(1, (1:5 - 3) / 2)
  b <- cbind(c(0.5, -1), c(1.5, 0.5))
  y <- rbind(matrix(rpois(30 * 5, exp(x %*% b[, 1])), 30, byrow = TRUE),
             matrix(rpois(30 * 5, exp(x %*% b[, 2])), 30, byrow = TRUE))
  y[3, 2] <- NA
  pr <- c(0.4, 0.6)
  lk <- sapply(1:2, function(k) {
    m <- exp(x %*% b[, k])
    apply(y, 1, function(v) sum(dpois(v, m, log = TRUE), na.rm = TRUE))
  })
  ll <- sum(log(exp(lk) %*% pr))
  z <- zgll(y, x, par = b, prob = pr)
  expect_equal(z$ll, ll, tolerance = 1e-12)
  po <- sweep(exp(lk), 2, pr, "*")
  expect_equal(z$post, po / rowSums(po), tolerance = 1e-12)
  f <- zgem(y, x, par = b + 0.2, prob = c(0.5, 0.5))
  expect_gte(f$ll, zgll(y, x, par = b + 0.2, prob = c(0.5, 0.5))$ll)
  expect_equal(f$ll, zgll(y, x, par = f$par, prob = f$prob)$ll,
               tolerance = 1e-10)
  expect_lt(max(abs(f$par - b)), 0.5)
  zi <- zgll(y, x, z = matrix(1, 5, 1), par = rbind(b, -20), prob = pr,
             off = matrix(0, 60, 5))
  expect_equal(zi$ll, ll, tolerance = 1e-6)
  g <- zgem(y, x, z = matrix(1, 5, 1), par = rbind(b, -2), prob = pr,
            maxit = 50)
  expect_identical(dim(g$par), c(3L, 2L))
  # the multivariate model
  yy <- array(rnorm(20 * 4 * 2), c(20, 4, 2))
  yy[1:10, , ] <- yy[1:10, , ] + 2
  yy[2, 3, 1] <- NA
  xm <- cbind(1, 1:4)
  be <- array(c(2, 0, 2, 0, 0, 0, 0, 0), c(2, 2, 2))
  sg <- array(c(1, 0.2, 0.2, 1, 1, 0, 0, 1), c(2, 2, 2))
  lk2 <- sapply(1:2, function(k) {
    s <- sg[, , k]
    l <- chol(s)
    sapply(1:20, function(i) {
      tot <- 0
      for (t in 1:4) {
        v <- yy[i, t, ]
        if (anyNA(v)) next
        e <- backsolve(l, v - as.vector(xm[t, ] %*% be[, , k]),
                       transpose = TRUE)
        tot <- tot - sum(log(diag(l))) - 0.5 * sum(e^2) - log(2 * pi)
      }
      tot
    })
  })
  m <- mvtll(yy, xm, be, sg, c(0.5, 0.5))
  expect_equal(m$ll, sum(log(exp(lk2) %*% c(0.5, 0.5))), tolerance = 1e-12)
  e <- mvtem(yy, xm, be, sg, c(0.5, 0.5))
  expect_gte(e$ll, m$ll)
  expect_identical(dim(e$beta), c(2L, 2L, 2L))
})

test_that("glass transfer, tool marks, grid Gi* and the pattern test", {
  par <- c(80, 50, 0.5, 0.05, 0.8, 0.9, 0.1, 0.2, 0.4, 0.6, 0.05, 0.1,
           0.5, 0.7, 1, 1, 0.5)
  sm <- par
  sm[1:2] <- c(10, 5)
  p <- tfrp(sm, ymax = 60, nsim = 50, seed = 1)
  expect_equal(sum(p$prob$p), 1, tolerance = 1e-6)
  fx <- par
  fx[5:14] <- rep(c(0.85, 0.15, 0.5, 0.075, 0.6), each = 2)
  q <- tfrp(fx, deffect = FALSE, ymax = 5, nsim = 1)
  expect_true(all(q$prob$se == 0))
  expect_identical(tfrp(fx, ymax = 5, nsim = 1, tdist = "unif")$prob$y,
                   0:5 + 0)
  set.seed(27)
  mk <- cumsum(rnorm(900))
  a <- mk[1:700] + rnorm(700, sd = 0.1)
  b <- mk[101:800] + rnorm(700, sd = 0.1)
  same <- chmb(a, b, wopt = 120, wval = 30)
  diff <- chmb(a, cumsum(rnorm(700)), wopt = 120, wval = 30)
  expect_lt(same[["p"]], 0.01)
  expect_gt(diff[["p"]], same[["p"]])
  # the classic Gi* against spdep's localG with the cell in its own
  # neighbourhood (queen, binary weights)
  set.seed(28)
  m <- matrix(rpois(48, 3), 6, 8)
  g <- gstr(m, rad = 1.5, nsim = 99, seed = 1)
  expect_identical(g$row[1:7], c(1:6, 1) + 0)
  skip_if_not_installed("spdep")
  nb <- spdep::include.self(spdep::cell2nb(6, 8, type = "queen"))
  # cell2nb numbers cells across rows first: row r, column c is
  # (r - 1) * 8 + c; gstr goes down the columns
  ix <- as.vector(outer(1:6, 1:8, function(r, k) (r - 1) * 8 + k))
  lg <- spdep::localG(as.vector(t(m)), spdep::nb2listw(nb, style = "B"))
  expect_equal(g$gz, as.vector(lg)[ix], tolerance = 1e-10)
  set.seed(29)
  ca <- rpois(40, 8)
  s <- sptst(ca, ca, nboot = 99, seed = 1)
  expect_identical(s$glob[["flagged"]], 0)
  expect_identical(s$glob[["s"]], 1)
  t2 <- sptst(c(ca, NA), c(rev(ca), 3), nboot = 99, seed = 1)
  expect_identical(t2$glob[["left"]], 1)
  expect_true(is.nan(t2$units$p[41]))
})

test_that("dossier functions refuse what they cannot use", {
  errs("fair", list(
    list("one label per value", 1:3, c(0, 1, 1), c("a", "b")),
    list("with no missing values", 1:2, c(0, 1), c("a", NA)),
    list("must be 0 or 1", 1:2, c(0, 2), c("a", "b")),
    list("must be 0 or 1", 1:2, "a", c("a", "b")),
    list("`cut` must be one number", 1:2, c(0, 1), 1:2, cut = NA),
    list("`conf` must be one number strictly between", 1:2, c(0, 1), 1:2,
         conf = 1)
  ))
  x <- matrix(1:8, 4)
  errs("tlvc", list(
    list("at least 2 items", x, rep(1, 4)),
    list("must be a numeric matrix", 1:4, 1:4),
    list("must hold finite numbers", x + c(NA, 0, 0, 0), c(1, 1, 2, 2))
  ))
  bg <- list(mu = 1:2, u = diag(2), c = diag(2))
  errs("tlvlr", list(
    list("the list tlvc", 1:2, 1:2, 1, 1, list(1)),
    list("must be square", 1:2, 1:2, 1, 1, list(1:2, diag(3), diag(2))),
    list("not positive definite", 1:2, 1:2, 1, 1,
         list(1:2, diag(2), -diag(2)))
  ))
  errs("aorwt", list(
    list("numbers of seconds or date-times", "a", 1),
    list("must not be infinite", Inf, 1),
    list("the length of `from`", 1, 1:2),
    list("one instant", 1, 2, origin = NA_real_),
    list("both `tt` and `toff`", 1, 2, tt = 1),
    list("0 <= lo < hi", 1, 2, lo = 5, hi = 1),
    list("strictly increasing", 1, 2, tt = c(2, 1), toff = c(0, 0))
  ))
  errs("bdig", list(
    list("`x` must be numbers", "1"),
    list("from 1 to 9", 1, nd = 10),
    list("from -1 to 22", 1, dp = 23)
  ))
  errs("bnft", list(
    list("from 1 to 6", 1, nd = 7),
    list("2-digit numbers", c(10, 5), nd = 2)
  ))
  errs("rptc", list(
    list("too many units times years", "a", 1, 1, 0, 1.5e9),
    list("must be whole numbers", "a", 1.5, 1, 0, 2)
  ))
  errs("rshl", list(list("`gam` must be at least 1", 1:3, 0.5)))
  errs("rswp", list(list("strictly between 0 and 0.5", 1:3, alpha = 0.6)))
  errs("lwss", list(list("`delta` must be one number from 0", 1:3, 1:3,
                         delta = -1)))
  errs("chmb", list(
    list("at most `wopt`", 1:10, 1:10, 3, 4),
    list("shorter than `wopt`", 1:10, 1:10, 50, 4),
    list("every window is flat", rep(1, 300), rep(1, 300), 50, 10)
  ))
  m <- matrix(1, 2, 2)
  errs("gstr", list(
    list("the dimensions of `cnt`", m, matrix(1, 3, 2)),
    list("at most 46340 rows", matrix(0, 46341, 1)),
    list("from 0 to 1000", m, rad = -1),
    list("whole numbers from 0", m + 0.5)
  ))
  errs("sptst", list(
    list("must be counts", "a", 1),
    list("the same length", 1:2, 1:3),
    list("whole numbers from 0", c(1, 2.5), c(1, 2))
  ))
  errs("tfrp", list(list("outside the model's domain", c(-1, rep(0.5, 16)))))
  y <- matrix(1, 3, 2)
  x <- matrix(1, 2, 1)
  errs("zgll", list(
    list("one row per period", y, matrix(1, 3, 1), par = 1, prob = 1),
    list("one row per period", y, x, z = matrix(1, 3, 1), par = 1:2,
         prob = 1),
    list("the dimensions of `y`", y, x, par = 1, prob = 1, off = x),
    list("one column per group", y, x, par = matrix(1, 1, 2), prob = 1),
    list("`prob` must sum to 1", y, x, par = matrix(1), prob = 0.5),
    list("whole counts", y + 0.5, x, par = matrix(1), prob = 1)
  ))
  errs("zgem", list(list("`maxit` must be one whole number from 1", y, x,
                         par = matrix(1), prob = 1, maxit = 0)))
  ya <- array(1, c(3, 2, 1))
  errs("mvtll", list(
    list("dimensions must agree", ya, x, array(1, c(1, 1, 1)),
         array(1, c(2, 2, 1)), 1),
    list("array of 3 dimensions", y, x, 1, 1, 1),
    list("not positive definite", ya, x, array(1, c(1, 1, 1)),
         array(-1, c(1, 1, 1)), 1)
  ))
  errs("mvtem", list(
    list("dimensions must agree", ya, x, array(1, c(1, 1, 1)),
         array(1, c(2, 2, 1)), 1),
    list("singular", ya, x, array(1, c(1, 1, 1)), array(1, c(1, 1, 1)), 1)
  ))
})

test_that("every dossier function answers one valid call in its shape", {
  x <- c(0.2, 0.7, 0.4, 0.9)
  calls <- list(
    list(quote(fair(x, c(0, 1, 0, 1), c(1, 1, 2, 2))),
         c("cnt", "met", "rat", "lo", "hi"), "list 5"),
    list(quote(tlvc(matrix(x, 4), c(1, 1, 2, 2))), c("mu", "u", "c"),
         "list 3"),
    list(quote(tlvlr(1, 1.1, 2, 2, list(1, matrix(0.1), matrix(1)))),
         c("ln", "log10"), "numeric 2"),
    list(quote(aorwt(0, 7200, nb = 4)), c("tot", "w", "cnt"), "list 3"),
    list(quote(bdig(x)), c("dig", "cnt"), "list 2"),
    list(quote(bnft(1:9)), c("obs", "stat"), "list 2"),
    list(quote(skp(x, 0.01, nsim = 9, seed = 1)), NULL, "integer 4"),
    list(quote(rptc(1:2, c(1, 1), c(12, 12), 1, 1)), c("unit", "const"),
         "data.frame 2"),
    list(quote(rswp(x)), c("p", "gcrit", "t"), "list 3"),
    list(quote(rshl(x)), c("gam", "lo", "hi"), "data.frame 1"),
    list(quote(zgll(matrix(1:4, 2), matrix(1, 2, 1), par = matrix(0),
                    prob = 1)), c("ll", "post"), "list 2"),
    list(quote(zgem(matrix(1:4, 2), matrix(1, 2, 1), par = matrix(0),
                    prob = 1)), c("par", "prob", "post", "ll", "iter"),
         "list 5"),
    list(quote(mvtll(array(x, c(2, 2, 1)), matrix(1, 2, 1),
                     array(0, c(1, 1, 1)), array(1, c(1, 1, 1)), 1)),
         c("ll", "post"), "list 2"),
    list(quote(mvtem(array(x, c(2, 2, 1)), matrix(1, 2, 1),
                     array(0, c(1, 1, 1)), array(1, c(1, 1, 1)), 1)),
         c("beta", "sig", "prob", "post", "ll", "iter"), "list 6"),
    list(quote(tfrp(c(80, 50, 0.5, 0.05, rep(0.5, 10), 1, 1, 0.5), ymax = 2,
                    nsim = 2, seed = 1)), c("prob", "cut"), "list 2"),
    list(quote(lwss(1:4, x)), NULL, "numeric 4"),
    list(quote(chmb(sin(1:200 / 7), sin(1:200 / 7 + 0.1), 40, 10)),
         c("nsame", "ndiff", "u", "p"), "numeric 4"),
    list(quote(gstr(matrix(1:4, 2))),
         c("row", "col", "gz", "pz", "qz", "pp", "ze", "pe", "qe"),
         "data.frame 4"),
    list(quote(sptst(1:4, 4:1, nboot = 9, seed = 1)), c("units", "glob"),
         "list 2")
  )
  fns <- vapply(calls, function(a) as.character(a[[1L]][[1L]]), "")
  expect_setequal(fns, dsx)
  for (k in seq_along(calls)) {
    v <- eval(calls[[k]][[1L]])
    expect_named(v, calls[[k]][[2L]])
    expect_identical(paste(class(v)[1L], NROW(v)), calls[[k]][[3L]],
                     info = fns[k])
  }
})

test_that("named per-case and per-unit arguments are matched by name", {
  s <- c(a = 0.2, b = 0.7, c = 0.4, d = 0.9)
  y <- c(a = 0, b = 1, c = 0, d = 1)
  g <- c(a = "x", b = "x", c = "y", d = "y")
  o <- c(4, 2, 1, 3)
  expect_identical(fair(s, y[o], g[o]), fair(s, y, g))
  expect_identical(fair(s, (y == 1)[o], g[o]), fair(s, y, g))
  expect_identical(lwss(s, y[o]), lwss(s, y))
  expect_identical(rptc(g, c(a = 1, b = 1, c = 1, d = 1)[o],
                        c(a = 12, b = 12, c = 3, d = 12)[o], 1, 1),
                   rptc(g, c(1, 1, 1, 1), c(12, 12, 3, 12), 1, 1))
  a <- c(u = 3, v = 5, w = 2)
  b <- c(w = 4, u = 1, v = 6)
  expect_identical(sptst(a, b, nboot = 9, seed = 1),
                   sptst(a, b[names(a)], nboot = 9, seed = 1))
  expect_identical(sptst(as.list(a), data.frame(t(b)), nboot = 9, seed = 1),
                   sptst(a, b[names(a)], nboot = 9, seed = 1))
  m <- matrix(b, 1, dimnames = list(NULL, names(b)))
  expect_identical(sptst(a, m, nboot = 9, seed = 1),
                   sptst(a, b[names(a)], nboot = 9, seed = 1))
  od <- diag(3)
  expect_identical(bjtc(list(p = 0.5, q = 0.3, r = 0.2), od, 1:3),
                   bjtc(c(0.5, 0.3, 0.2), od, 1:3))
  errs("sptst", list(
    list("a list of single numbers", list(1:2, 3), 1:2)
  ))
  errs("fair", list(
    list("a vector of numbers or labels", s, c(b = 1i, a = 0i, c = 0i,
                                                d = 1i), g)
  ))
})
