# gmixlr(): the gamma model by enumeration in R (dgamma, pgamma) at given
# parameters, an exact ratio, the fit on a clean mixture, and refusals.

# one marker: kept alleles p (the rest, Q, last), heights y (0: no peak),
# the known contributors' genotypes (index pairs), nu unknowns, the
# proportions m (knowns first); unknowns drawn with coancestry theta given
# the counts cnt of ncnt conditioning alleles
lik1 <- function(
  p,
  y,
  known,
  nu,
  m,
  mu,
  om,
  at,
  pc,
  lam,
  theta = 0,
  cnt = 0 * p,
  ncnt = 0
) {
  n <- length(p)
  gs <- which(upper.tri(diag(n), diag = TRUE), arr.ind = TRUE)
  cb <- if (nu) {
    as.matrix(expand.grid(rep(list(seq_len(nrow(gs))), nu)))
  } else {
    matrix(0L, 1, 0)
  }
  sc <- mu * om^2
  tot <- 0
  for (r in seq_len(nrow(cb))) {
    g <- c(known, lapply(cb[r, ], function(k) gs[k, ]))
    dose <- numeric(n)
    for (k in seq_along(g)) {
      for (a in g[[k]]) {
        dose[a] <- dose[a] + m[k]
      }
    }
    sh <- dose / om^2
    l <- 1
    nd <- 0
    for (a in seq_len(n)) {
      if (y[a] > 0 && sh[a] > 0) {
        l <- l * dgamma(y[a], sh[a], scale = sc)
      } else if (y[a] > 0) {
        l <- l * pc * p[a] * lam * exp(-lam * (y[a] - at))
        nd <- nd + 1
      } else if (sh[a] > 0) {
        l <- l * pgamma(at, sh[a], scale = sc)
      }
    }
    if (nd == 0) {
      l <- l * (1 - pc)
    }
    pr <- 1
    cc <- cnt
    nn <- ncnt
    for (k in seq_len(nu)) {
      al <- gs[cb[r, k], ]
      for (a in al) {
        pr <- pr * (cc[a] * theta + (1 - theta) * p[a]) / (1 + (nn - 1) * theta)
        cc[a] <- cc[a] + 1
        nn <- nn + 1
      }
      if (al[1] != al[2]) pr <- 2 * pr
    }
    tot <- tot + l * pr
  }
  tot
}

fr <- data.frame(locus = "L", allele = 1:4, freq = c(0.1, 0.2, 0.3, 0.25))

test_that("a single source without drop-in gives 1 / (2 p q), fitted or not", {
  pk <- data.frame(locus = "L", allele = c(1, 3), height = c(800, 900))
  ref <- data.frame(id = "S", locus = "L", a1 = 1, a2 = 3)
  r <- gmixlr(pk, fr, ref, hp = "S", nc = 1, pc = 0)
  expect_equal(r$lr, 1 / (2 * 0.1 * 0.3), tolerance = 1e-8)
  expect_true(all(r$fit$conv))
  # the fitted mu and omega: those of two gamma peaks of one dose
  expect_equal(r$fit$mu[1], 850, tolerance = 1e-3)
  r2 <- gmixlr(
    pk,
    fr,
    ref,
    hp = "S",
    nc = 1,
    pc = 0,
    par = c(mu = 1000, omega = 0.2, m1 = 1)
  )
  expect_equal(r2$lr, 1 / (2 * 0.1 * 0.3), tolerance = 1e-12)
  expect_true(all(is.na(r2$fit$conv)))
})

test_that("drop-out and drop-in match the enumeration", {
  pk <- data.frame(locus = "L", allele = 1, height = 900)
  ref <- data.frame(id = "S", locus = "L", a1 = 1, a2 = 2)
  par <- c(mu = 1000, omega = 0.3, m1 = 1)
  r <- gmixlr(pk, fr, ref, hp = "S", nc = 1, par = par)
  # kept: 1 (seen), 2 (the suspect's); the rest Q = 0.7
  p <- c(0.1, 0.2, 0.7)
  hp <- lik1(p, c(900, 0, 0), list(c(1, 2)), 0, 1, 1000, 0.3, 50, 0.05, 0.01)
  hd <- lik1(p, c(900, 0, 0), list(), 1, 1, 1000, 0.3, 50, 0.05, 0.01)
  expect_equal(r$fit$loglik, log(c(hp, hd)), tolerance = 1e-13)
  expect_equal(r$lr, hp / hd, tolerance = 1e-12)
  expect_equal(sum(r$locus$log10lr), r$log10lr, tolerance = 1e-14)
})

test_that("two contributors with coancestry match the enumeration", {
  pk <- data.frame(locus = "L", allele = 1:3, height = c(700, 650, 500))
  ref <- data.frame(id = c("V", "S"), locus = "L", a1 = c(1, 3), a2 = c(2, 3))
  par <- c(mu = 900, omega = 0.25, m1 = 0.6, m2 = 0.4)
  r <- gmixlr(
    pk,
    fr,
    ref,
    hp = c("V", "S"),
    nc = 2,
    hd = "V",
    theta = 0.02,
    par = par
  )
  p <- c(0.1, 0.2, 0.3, 0.4)
  hp <- lik1(
    p,
    c(700, 650, 500, 0),
    list(c(1, 2), c(3, 3)),
    0,
    c(0.6, 0.4),
    900,
    0.25,
    50,
    0.05,
    0.01
  )
  hd <- lik1(
    p,
    c(700, 650, 500, 0),
    list(c(1, 2)),
    1,
    c(0.6, 0.4),
    900,
    0.25,
    50,
    0.05,
    0.01,
    theta = 0.02,
    cnt = c(1, 1, 2, 0),
    ncnt = 4
  )
  expect_equal(r$lr, hp / hd, tolerance = 1e-12)
  expect_identical(r$mix$who, c("V", "S", "V", "U1"))
})

test_that("stutter and degradation match the formula", {
  f2 <- data.frame(
    locus = "L",
    allele = 9:11,
    freq = c(0.2, 0.4, 0.3),
    size = c(146, 150, 154)
  )
  pk <- data.frame(
    locus = "L",
    allele = c(10, 11, 9),
    height = c(1500, 300, 120)
  )
  ref <- data.frame(
    id = c("A", "B"),
    locus = "L",
    a1 = c(10, 10),
    a2 = c(10, 11)
  )
  mu <- 1000
  om <- 0.25
  be <- 0.7
  xi <- 0.08
  r <- gmixlr(
    pk,
    f2,
    ref,
    hp = c("A", "B"),
    nc = 2,
    hd = c("A", "B"),
    stut = TRUE,
    deg = TRUE,
    par = c(mu = mu, omega = om, beta = be, xi = xi, m1 = 0.7, m2 = 0.3)
  )
  c10 <- (0.7 * 2 + 0.3) * be^0.25
  c11 <- 0.3 * be^0.29
  sh <- c((1 - xi) * c10 + xi * c11, (1 - xi) * c11, xi * c10) / om^2
  want <- sum(dgamma(c(1500, 300, 120), sh, scale = mu * om^2, log = TRUE)) +
    log(0.95)
  expect_equal(r$fit$loglik, c(want, want), tolerance = 1e-13)
  expect_identical(r$lr, 1)
})

test_that("the fit of a clean two-person mixture finds the minor contributor", {
  peaks <- data.frame(
    locus = rep(c("D3", "vWA"), c(4, 3)),
    allele = c(14, 15, 16, 17, 16, 17, 18),
    height = c(1210, 380, 1150, 420, 1400, 1960, 610)
  )
  freq <- data.frame(
    locus = rep(c("D3", "vWA"), c(5, 5)),
    allele = c(13, 14, 15, 16, 17, 15, 16, 17, 18, 19),
    freq = c(0.05, 0.1, 0.25, 0.25, 0.2, 0.1, 0.2, 0.3, 0.2, 0.1)
  )
  ref <- data.frame(
    id = rep(c("V", "S"), each = 2),
    locus = rep(c("D3", "vWA"), 2),
    a1 = c(14, 16, 15, 17),
    a2 = c(16, 17, 17, 18)
  )
  r <- gmixlr(peaks, freq, ref, hp = c("V", "S"), nc = 2, hd = "V")
  expect_true(all(r$fit$conv))
  # the minor peaks are about 0.27 of the total
  expect_equal(r$mix$m[2], 0.27, tolerance = 0.03)
  # with peaks this clear the unknown must be the suspect's genotype:
  # the ratio is about 1 / (2 * 0.25 * 0.2 * 2 * 0.3 * 0.2)
  expect_equal(r$log10lr, -log10(0.1 * 0.12), tolerance = 0.01)
  expect_output(print(r), "gmixlr: 2 markers")
})

test_that("gmixlr() refuses bad input", {
  pk <- data.frame(locus = "L", allele = 1, height = 900)
  ref <- data.frame(id = "S", locus = "L", a1 = 1, a2 = 2)
  errs(
    "gmixlr",
    list(
      list(
        "columns locus, allele and height",
        peaks = 1,
        freq = fr,
        ref = ref,
        hp = "S",
        nc = 1
      ),
      list(
        "peaks\\$height` must be numbers",
        peaks = transform(pk, height = "a"),
        freq = fr,
        nc = 1
      ),
      list(
        "columns locus, allele and freq",
        peaks = pk,
        freq = fr[1:2],
        nc = 1
      ),
      list(
        "freq\\$size` must be numbers",
        peaks = pk,
        freq = transform(fr, size = "x"),
        nc = 1
      ),
      list(
        "columns id, locus, a1 and a2",
        peaks = pk,
        freq = fr,
        ref = ref[1:2],
        nc = 1
      ),
      list("TRUE or FALSE", peaks = pk, freq = fr, nc = 1, stut = NA),
      list("named numbers", peaks = pk, freq = fr, nc = 1, par = 1),
      list(
        "ids of profiles",
        peaks = pk,
        freq = fr,
        ref = ref,
        hp = list("S"),
        nc = 1
      ),
      list("whole number from 1 to 5", peaks = pk, freq = fr, nc = 6),
      list(
        "`at` must be one number above 0",
        peaks = pk,
        freq = fr,
        nc = 1,
        at = 0
      ),
      list("`pc` must be one number", peaks = pk, freq = fr, nc = 1, pc = 1),
      list(
        "`theta` must be one number",
        peaks = pk,
        freq = fr,
        nc = 1,
        theta = -1
      ),
      list("needs a column size", peaks = pk, freq = fr, nc = 1, deg = TRUE),
      list("below `at`", peaks = transform(pk, height = 20), freq = fr, nc = 1),
      list(
        "which `freq` does not list",
        peaks = transform(pk, allele = 9),
        freq = fr,
        nc = 1
      ),
      list(
        "no profile in `ref`",
        peaks = pk,
        freq = fr,
        ref = ref,
        hp = "T",
        nc = 1
      ),
      list(
        "more contributors than `nc`",
        peaks = pk,
        freq = fr,
        ref = rbind(ref, transform(ref, id = "T")),
        hp = c("S", "T"),
        nc = 1
      ),
      list(
        "must name mu, omega, m1",
        peaks = pk,
        freq = fr,
        nc = 1,
        par = c(mu = 1)
      ),
      list(
        "sum to 1",
        peaks = pk,
        freq = fr,
        nc = 2,
        par = c(mu = 1, omega = 1, m1 = 0.5, m2 = 0.4)
      ),
      list(
        "numeric alleles",
        peaks = transform(pk, allele = "x"),
        freq = transform(fr, allele = c("1", "2", "3", "x")),
        nc = 1,
        stut = TRUE
      ),
      list("must have a peak", peaks = pk[0, ], freq = fr, nc = 1),
      list(
        "no genotype for S",
        peaks = pk,
        freq = rbind(fr, transform(fr, locus = "M")),
        ref = ref,
        hp = "S",
        nc = 1
      )
    )
  )
})

test_that("gmixlr() equals EuroForMix's likelihood and its maximum", {
  # EuroForMix (its version and commit in ref/efm.txt): the log likelihood
  # of each account at given parameters and its maximum, stored by
  # dev/ref/efm.R. EuroForMix keeps its genotype weights in single
  # precision, about 1e-7 of each, so the two agree to 1e-6, not 1e-12.
  ref <- refget("efm")
  for (k in setdiff(names(ref), c("pkgs", "commit"))) {
    x <- ref[[k]]
    a <- list(x$peaks, x$freq, x$ref, hp = x$hp, nc = x$nc,
              hd = if (length(x$hd)) x$hd, pc = x$pc, lambda = x$lambda,
              theta = x$theta, stut = x$stut, deg = x$deg)
    r <- do.call(gmixlr, c(a, list(par = x$par)))
    expect_lt(max(abs(r$fit$loglik - x$ll)), 1e-6, label = k)
    if (!anyNA(x$mle)) {
      f <- do.call(gmixlr, a)
      expect_lt(max(abs(f$fit$loglik - x$mle)), 1e-6, label = k)
    }
  }
})
