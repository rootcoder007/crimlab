# cmpsc(): a shifted copy, the algorithm written out in R on a small case,
# cmpsR's scores (ref/cmps.txt), and refusals.

# the algorithm in plain R: segments of len from the first to the last
# value of a; lags keeping a segment inside b; local maxima; the
# consistent peak; votes
cmps_r <- function(a, b, len, npeaks, tx) {
  f0 <- min(which(!is.na(a)))
  f1 <- max(which(!is.na(a)))
  ns <- (f1 - f0 + 1) %/% len
  pk <- function(lo, m, k) {
    lags <- (0:(length(b) - m)) - lo + 1
    cc <- vapply(
      lags,
      function(d) {
        x <- a[lo:(lo + m - 1)]
        y <- b[(lo:(lo + m - 1)) + d]
        ok <- !is.na(x) & !is.na(y)
        if (sum(ok) < 2 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
          return(NA_real_)
        }
        stats::cor(x[ok], y[ok])
      },
      0
    )
    nj <- length(cc)
    isp <- vapply(
      seq_len(nj),
      function(j) {
        !is.na(cc[j]) &&
          (j == 1 || is.na(cc[j - 1]) || cc[j] > cc[j - 1]) &&
          (j == nj || is.na(cc[j + 1]) || cc[j] >= cc[j + 1])
      },
      NA
    )
    j <- which(isp)
    j <- j[order(-cc[j])]
    lags[utils::head(j, k)]
  }
  cand <- vector("list", ns)
  for (s in seq_len(ns)) {
    i0 <- f0 + (s - 1) * len
    lv <- list()
    for (l in seq_along(npeaks)) {
      m <- len * 2^(l - 1)
      if (m > f1 - f0 + 1) {
        next
      }
      lo <- i0 + len %/% 2 - m %/% 2
      lo <- max(f0, min(lo, f1 + 1 - m))
      lv[[length(lv) + 1L]] <- pk(lo, m, npeaks[l])
    }
    if (length(npeaks) == 1L) {
      cand[[s]] <- lv[[1L]]
    } else {
      top <- lv[[length(lv)]]
      ok <- vapply(
        top,
        function(p) {
          all(vapply(lv[-length(lv)], function(q) any(abs(q - p) <= tx), NA))
        },
        NA
      )
      cand[[s]] <- utils::head(top[ok], 1L)
    }
  }
  all <- sort(unique(unlist(cand)))
  votes <- vapply(
    all,
    function(p) sum(vapply(cand, function(cs) any(abs(cs - p) <= tx), NA)),
    0
  )
  tied <- all[votes == max(votes)]
  pos <- tied[(length(tied) + 1L) %/% 2L]
  hit <- vapply(cand, function(cs) any(abs(cs - pos) <= tx), NA)
  list(score = sum(hit), nseg = ns, pos = pos, cmps = hit)
}

test_that("a shifted copy is found at its lag", {
  set.seed(1)
  a <- cumsum(rnorm(600, sd = 0.5)) + sin(seq_len(600) / 9)
  b <- a[31:600] + rnorm(570, sd = 0.1)
  r <- cmpsc(a, b)
  expect_identical(r$pos, -30L)
  expect_identical(r$nseg, 12L)
  # the first segment would start before b at lag -30
  expect_false(r$segs$cmps[1])
  expect_gte(r$score, 10L)
  expect_equal(r$scaled, r$score / 12)
  expect_output(print(r), "of 12 segments congruent")
})

test_that("cmpsc() follows the algorithm written out in R", {
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
      r <- cmpsc(a, b, len = 30, npeaks = np, tx = 8)
      w <- cmps_r(a, b, 30, np, 8)
      expect_identical(r$nseg, as.integer(w$nseg))
      expect_identical(r$score, as.integer(w$score))
      expect_identical(r$pos, as.integer(w$pos))
      expect_identical(r$segs$cmps, w$cmps)
      expect_identical(r$segs$start, as.integer(5 + 30 * (seq_len(w$nseg) - 1)))
    }
  }
})

test_that("cmpsc() is within one segment of cmpsR", {
  # cmpsR's extract_feature_cmps() at its defaults on simulated signatures
  # (pairs 1 to 4 from one source), stored in ref/cmps.txt by
  # dev/ref/cmps.R (version named there). The two differ in how far a
  # segment may hang over the end of the reference.
  ref <- refget("cmps")
  for (k in 1:8) {
    x <- ref[[paste0("p", k)]]
    r <- cmpsc(x$a, x$b)
    expect_identical(r$nseg, as.integer(x$nseg), info = k)
    expect_lte(abs(r$score - x$score), 1, label = paste("score", k))
    if (k <= 4) expect_lte(abs(r$pos - x$pos), 25, label = paste("lag", k))
  }
})

test_that("cmpsc() refuses bad input", {
  a <- sin(1:200)
  errs(
    "cmpsc",
    list(
      list("`a` must be numbers", a = c(1, Inf, 3), b = a),
      list("`b` must be numbers", a = a, b = "x"),
      list("one whole number from 2", a = a, b = a, len = 1),
      list("at least `len` points", a = c(NA, a[1:20], NA), b = a, len = 30),
      list("one to eight whole numbers from 1", a = a, b = a, npeaks = 0),
      list("one whole number from 0", a = a, b = a, tx = -1)
    )
  )
})
