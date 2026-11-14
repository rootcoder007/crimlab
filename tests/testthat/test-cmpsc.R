# cmpsc(): cmpsR's results exactly (ref/cmps.txt), the full = TRUE rule
# written out in R, a shifted copy, and refusals.

test_that("cmpsc() reproduces cmpsR exactly", {
  # cmpsR's extract_feature_cmps() on simulated signatures and edge cases
  # (different lengths, missing runs at the ends and inside, a flat
  # stretch, ties, one level of peaks, two cases it refuses), stored in
  # ref/cmps.txt by dev/ref/cmps.R (version named there)
  ref <- refget("cmps")
  for (k in setdiff(names(ref), "pkgs")) {
    x <- ref[[k]]
    if (x$err) {
      expect_error(cmpsc(x$a, x$b, npeaks = x$npeaks), "^cmpsc: ", info = k)
      next
    }
    r <- cmpsc(x$a, x$b, npeaks = x$npeaks)
    expect_identical(r$score, as.integer(x$score), info = k)
    expect_identical(r$nseg, as.integer(x$nseg), info = k)
    expect_identical(r$pos, as.numeric(x$pos), info = k)
    # each segment's candidate nearest the congruent lag
    near <- vapply(x$ccp, function(z) {
      if (!length(z)) NA_real_ else z[which.min(abs(z - x$pos))]
    }, 0)
    expect_identical(r$segs$ccp, near, info = k)
  }
})

# full = TRUE in plain R: segments of len from the first to the last value
# of a, widened about their centre and cut at the ends of a; lags keeping
# a widened segment inside b; local maxima at least their neighbours; the
# peaks of the last level every level confirms; votes over the candidates
# with repeats, the middle (ceiling) of the most voted
cmps_r <- function(a, b, len, npeaks, tx) {
  f0 <- min(which(!is.na(a)))
  f1 <- max(which(!is.na(a)))
  ns <- (f1 - f0 + 1) %/% len
  pk <- function(lo, hi, k) {
    m <- hi - lo + 1
    lags <- (0:(length(b) - m)) - (lo - 1)
    cc <- vapply(lags, function(d) {
      x <- a[lo:hi]
      y <- b[(lo:hi) + d]
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) < 2 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
        return(NA_real_)
      }
      stats::cor(x[ok], y[ok])
    }, 0)
    i <- which(!is.na(cc))
    v <- cc[i]
    n <- length(v)
    top <- if (n == 1L) {
      i[v >= 1]
    } else {
      i[vapply(seq_len(n), function(j) {
        (j == 1L || v[j] >= v[j - 1L]) && (j == n || v[j] >= v[j + 1L])
      }, NA)]
    }
    lags[utils::head(top[order(-cc[top])], k)]
  }
  cand <- lapply(seq_len(ns), function(s) {
    i0 <- f0 + (s - 1) * len
    lv <- lapply(seq_along(npeaks), function(l) {
      tc <- len * 2^(l - 1) - len
      pk(max(i0 - tc %/% 2, 1), min(i0 + len - 1 + ceiling(tc / 2), length(a)),
         npeaks[l])
    })
    last <- lv[[length(lv)]]
    if (length(npeaks) == 1L) return(last)
    last[vapply(last, function(p) {
      all(vapply(lv[-length(lv)], function(q) any(abs(q - p) <= tx), NA))
    }, NA)]
  })
  all <- sort(unlist(cand))
  votes <- vapply(all, function(p) {
    sum(vapply(cand, function(z) any(abs(z - p) <= tx), NA))
  }, 0)
  tied <- all[votes == max(votes)]
  pos <- tied[ceiling(length(tied) / 2)]
  hit <- vapply(cand, function(z) any(abs(z - pos) <= tx), NA)
  list(score = sum(hit), nseg = ns, pos = pos, cmps = hit)
}

test_that("full = TRUE follows the rule written out in R", {
  set.seed(3)
  for (k in 1:3) {
    a <- cumsum(rnorm(330, sd = 0.5))
    b <- if (k < 3) {
      a[(10 * k + 1):330] + rnorm(330 - 10 * k, sd = 0.3)
    } else {
      cumsum(rnorm(300))
    }
    a[1:4] <- NA
    b[7] <- NA
    for (np in list(c(3, 1), 4)) {
      r <- cmpsc(a, b, len = 30, npeaks = np, tx = 8, full = TRUE)
      w <- cmps_r(a, b, 30, np, 8)
      expect_identical(r$nseg, as.integer(w$nseg))
      expect_identical(r$score, as.integer(w$score))
      expect_identical(r$pos, as.numeric(w$pos))
      expect_identical(r$segs$cmps, w$cmps)
      expect_identical(r$segs$start, as.integer(5 + 30 * (seq_len(w$nseg) - 1)))
    }
  }
})

test_that("a shifted copy is found at its lag", {
  set.seed(1)
  a <- cumsum(rnorm(600, sd = 0.5)) + sin(seq_len(600) / 9)
  b <- a[31:600] + rnorm(570, sd = 0.1)
  r <- cmpsc(a, b)
  expect_identical(r$pos, -30)
  expect_identical(r$nseg, 12L)
  expect_identical(r$score, 11L)
  expect_equal(r$scaled, 11 / 12)
  # with full overlap the first segment cannot start before b
  f <- cmpsc(a, b, full = TRUE)
  expect_false(f$segs$cmps[1])
  expect_output(print(r), "of 12 segments congruent")
})

test_that("cmpsc() refuses bad input", {
  a <- sin(1:200)
  errs("cmpsc", list(
    list("`a` must be numbers", a = c(1, Inf, 3), b = a),
    list("`b` must be numbers", a = a, b = "x"),
    list("one whole number from 2", a = a, b = a, len = 1),
    list("at least `len` points", a = c(NA, a[1:20], NA), b = a, len = 30),
    list("one to eight whole numbers from 1", a = a, b = a, npeaks = 0),
    list("one whole number from 0", a = a, b = a, tx = -1),
    list("TRUE or FALSE", a = a, b = a, full = NA),
    list("shorter than a widened segment", a = a, b = a[1:150])
  ))
})
