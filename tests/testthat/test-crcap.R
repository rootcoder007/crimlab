# crcap() and dkfig(): against the formulas written out, a log-linear
# Poisson fit by glm() for three lists, a grid search over the error rates
# for the survey bounds, and the refusals.

test_that("two lists give Petersen, Chapman and the dependence box", {
  r <- crcap(c(n1 = 500, n2 = 400, m = 80), kappa = 1.5)
  e <- r$est
  expect_identical(r$type, "two")
  expect_identical(e[["petersen"]], 2500)
  expect_equal(e[["lower"]], 2500 / 1.5, tolerance = 1e-15)
  expect_equal(e[["upper"]], 3750, tolerance = 1e-15)
  expect_equal(e[["chapman"]], 501 * 401 / 81 - 1, tolerance = 1e-15)
  expect_equal(
    e[["se"]],
    sqrt(501 * 401 * 420 * 320 / (81^2 * 82)),
    tolerance = 1e-15
  )
  expect_identical(e[["floor"]], 820)
  # names in another order, or none
  expect_identical(
    crcap(c(m = 80, n2 = 400, n1 = 500))$est,
    crcap(c(500, 400, 80))$est
  )
  expect_null(r$pairs)
  expect_output(print(r), "two lists")
})

test_that("three lists match a log-linear fit with no three-way term", {
  cells <- c(
    "100" = 120,
    "010" = 90,
    "001" = 60,
    "110" = 30,
    "101" = 20,
    "011" = 15,
    "111" = 8
  )
  r <- crcap(rev(cells), int = c(-0.5, 0, 1))
  d <- data.frame(
    a = c(1, 0, 0, 1, 1, 0, 1),
    b = c(0, 1, 0, 1, 0, 1, 1),
    c = c(0, 0, 1, 0, 1, 1, 1),
    y = unname(cells)
  )
  f <- glm(
    y ~ (a + b + c)^2,
    family = poisson,
    data = d,
    control = glm.control(epsilon = 1e-14, maxit = 100)
  )
  m0 <- exp(unname(coef(f)[1]))
  expect_equal(r$est[["missing"]], m0, tolerance = 1e-10)
  expect_identical(r$est[["observed"]], sum(cells))
  expect_equal(
    r$sens$missing,
    r$est[["missing"]] * exp(-c(-0.5, 0, 1)),
    tolerance = 1e-15
  )
  expect_equal(r$sens$total, sum(cells) + r$sens$missing, tolerance = 1e-15)
  # the pairs: list 1 holds 100, 110, 101, 111, and so on
  n1 <- 120 + 30 + 20 + 8
  n2 <- 90 + 30 + 15 + 8
  m <- 30 + 8
  p <- r$pairs
  expect_identical(p$lists, c("1-2", "1-3", "2-3"))
  expect_equal(p$petersen[1], n1 * n2 / m, tolerance = 1e-15)
  expect_equal(
    p$chapman[1],
    (n1 + 1) * (n2 + 1) / (m + 1) - 1,
    tolerance = 1e-15
  )
  expect_identical(p$floor[1], n1 + n2 - m)
  expect_true(all(p$petersen >= p$floor))
  expect_output(print(r), "three lists")
  expect_null(crcap(cells)$sens)
})

test_that("crcap() refuses bad input in words", {
  cells <- c(
    "100" = 1,
    "010" = 1,
    "001" = 1,
    "110" = 1,
    "101" = 1,
    "011" = 1,
    "111" = 1
  )
  errs(
    "crcap",
    list(
      list("c\\(n1, n2, m\\) for two lists", c(1, 2)),
      list("`x` must be at least 0", c(1, -2, 1)),
      list("`kappa` must be one number from 1", c(5, 5, 1), kappa = 0.5),
      list("overlap m above 0", c(5, 5, 0)),
      list("overlap m above 0", c(5, 5, 6)),
      list("named n1, n2 and m", c(a = 5, n2 = 5, m = 1)),
      list(
        "`x` names n1 twice",
        stats::setNames(c(5, 5, 1), c("n1", "n1", "m"))
      ),
      list("must name its seven cells", unname(cells)),
      list(
        "must name its seven cells",
        setNames(cells, c(names(cells)[-7], "000"))
      ),
      list("110, 101, 011\\) must be above 0", replace(cells, 5, 0)),
      list("`int` goes with three lists", c(5, 5, 1), int = 0),
      list("`int` must be finite", cells, int = NA)
    )
  )
})

test_that("the survey bounds hold every admissible error and only those", {
  r <- dkfig(
    v = c(0.12, 0.08, 0.02),
    r = c(0.03, 0.03, 0.05),
    amax = 0.02,
    bmax = 0.25,
    thr = 0.15
  )
  b <- r$bnd
  lo <- pmax(b$r, (b$v - 0.02) / 0.98)
  hi <- pmax(lo, pmin(1, b$v / 0.75))
  expect_equal(b$lo, lo, tolerance = 1e-15)
  expect_equal(b$hi, hi, tolerance = 1e-15)
  expect_equal(b$dlo, lo - b$r, tolerance = 1e-15)
  expect_equal(b$ratio, hi / b$r, tolerance = 1e-15)
  expect_equal(b$brk, 1 - b$v / 0.15, tolerance = 1e-15)
  # the true rates v with v (1 - beta) + alpha (1 - v) = vobs over a grid of
  # the error rates: their range is the interval (where r does not bind)
  g <- expand.grid(
    a = seq(0, 0.02, length.out = 201),
    b = seq(0, 0.25, length.out = 201)
  )
  t <- (0.12 - g$a) / (1 - g$a - g$b)
  expect_equal(range(t), c(b$lo[1], b$hi[1]), tolerance = 1e-12)
  expect_null(r$cmp)
  expect_null(r$hier)
  # no false positives or negatives: the interval is the survey rate
  p <- dkfig(v = 0.1, r = 0)$bnd
  expect_identical(c(p$lo, p$hi), c(0.1, 0.1))
  expect_true(is.na(p$ratio))
  expect_false("brk" %in% names(p))
  expect_output(print(r), "bounds on the dark figure")
})

test_that("a partial survey bounds the total; incidents bound offences", {
  r <- dkfig(va = 900, ra = 300, rb = 450, lam = 0.5)$cmp
  k <- 0.5 * 300 / 900
  expect_equal(r[["rate"]], 1 / 3, tolerance = 1e-15)
  expect_equal(r[["tothi"]], 900 + 450 / k, tolerance = 1e-15)
  expect_identical(r[["totlo"]], 1350)
  expect_identical(r[["dklo"]], 600)
  expect_equal(r[["dkhi"]], 600 + 450 * (1 / k - 1), tolerance = 1e-15)
  expect_equal(r[["width"]], r[["dkhi"]] - r[["dklo"]], tolerance = 1e-14)
  expect_equal(r[["rateup"]], 750 / 1350, tolerance = 1e-15)
  h <- dkfig(k = c(1, 1, 2, 1, 3))$hier
  expect_equal(unname(h), c(5, 8, 1.6, 0.6, 5, 15), tolerance = 1e-15)
  expect_named(h, c("incid", "offen", "ratio", "extra", "lower", "upper"))
})

test_that("dkfig() refuses bad input in words", {
  errs(
    "dkfig",
    list(
      list("give `v` and `r`"),
      list("`r` must give the recorded rate", v = 0.1),
      list("`v` must be rates from 0 to 1", v = 1.5, r = 0),
      list("`r` must be rates from 0 to 1", v = 0.5, r = 2),
      list(
        "length 1 or the length of `v`, 3",
        v = c(0.1, 0.2, 0.3),
        r = c(0, 0)
      ),
      list("`amax` must be one number from 0 to 1", v = 0.1, r = 0, amax = -1),
      list(
        "`amax \\+ bmax` must be below 1",
        v = 0.1,
        r = 0,
        amax = 0.5,
        bmax = 0.5
      ),
      list("`thr` must lie above every `v`", v = 0.1, r = 0, thr = 0.05),
      list("`thr` must lie above every `v`", v = 0, r = 0, thr = 0.05),
      list("`ra` and `rb` must come with `va`", va = 10),
      list("`ra` must not exceed `va`", va = 10, ra = 20, rb = 1),
      list("`lam` must be above 0", va = 10, ra = 2, rb = 1, lam = 0),
      list("`rb` must be one number from 0", va = 10, ra = 2, rb = -1),
      list("whole numbers of offences", k = c(1, 1.5)),
      list("`k` must be at least 1", k = 0)
    )
  )
})
