# stgdid(): the group-time effects against 'did' (values stored),
# Wooldridge's cells against lm() with cell dummies, DID_M by hand, the
# bootstrap's seed and its errors against the Monte Carlo spread, and the
# refusals.

pd <- pepanel()

test_that("the group-time effects and their averages are did's", {
  # did 2.5.1, dev/ref-pe.R: att_gt(control_group = "notyettreated"), the
  # cells in did's order (cohort, then period)
  s <- stgdid(pd$rate, pd$town, pd$year, pd$first, reps = 0)
  want <- c(0.260336113194107, 0.150498601071501, -1.75809852807957,
            -3.10019012209789, -5.23579490387038, -6.08157947755966,
            -6.95640032485601, -7.60408538599005, -8.83412225340697,
            0.499322414698696, -0.60175672721392, 0.707812023184348,
            -0.692547352536337, -2.06324636359564, -3.46253515469553,
            -4.12104476213233, -6.1377535168619, -5.74978251375019,
            -0.339726442892985, 0.0443917915405688, 0.21417708570105,
            0.960864616246274, -0.14266670157196, -2.03270962146849,
            0.0214861980633371, -1.34107853231906, -2.60720506966173)
  expect_equal(s$att$est, want, tolerance = 1e-12)
  expect_identical(s$att$period[1:2], c(2016, 2017))
  # the simple and dynamic averages, by cohort size
  post <- s$att$evtm >= 0
  expect_equal(s$overall$est, mean(s$att$est[post]), tolerance = 1e-14)
  expect_equal(s$dynamic$est[s$dynamic$evtm == 0],
               mean(s$att$est[s$att$evtm == 0]), tolerance = 1e-14)
  expect_true(all(is.na(s$att$se)))
  never <- stgdid(pd$rate, pd$town, pd$year, pd$first, ctrl = "never",
                  reps = 0)
  expect_equal(
    never$att$est[1:6],
    c(0.376241724555213, -0.0423476529835387, -1.41235261224755,
      -2.65382523237464, -4.84080174899979, -6.3641561965119),
    tolerance = 1e-12
  )
  expect_output(print(never), "against units never")
})

test_that("covariates give did's doubly robust, IPW and regression cells", {
  # did 2.5.1: att_gt(xformla = ~x, est_method = ...), the first six cells;
  # the logistic fits converge to 1e-10 here and to glm()'s 1e-8 in did
  want <- list(
    dr = c(0.328729715647284, 0.746126415832205, -1.88110444560276,
           -2.94999808307653, -4.82916952649365, -5.3577828593876),
    ipw = c(0.32424216413549, 0.494870140421211, -1.819094493116,
            -2.96234804713758, -4.96138360422102, -5.63295764504193),
    reg = c(0.272953567622261, 0.766551182897008, -1.92723624856659,
            -3.06650448478301, -4.88462336309318, -5.35069245612971)
  )
  for (m in names(want)) {
    s <- stgdid(pd$rate, pd$town, pd$year, pd$first, adj = pd["x"], est = m,
                reps = 0)
    expect_equal(s$att$est[1:6], want[[m]], tolerance = 1e-7, info = m)
  }
  # a covariate that separates the cohort: the propensity fails, NA cells
  sep <- as.numeric(!is.na(pd$first) & pd$first == 2018)
  s <- stgdid(pd$rate, pd$town, pd$year, pd$first, adj = sep, est = "ipw",
              reps = 0)
  expect_true(all(is.na(s$att$est[s$att$cohort == 2018])))
  s <- stgdid(pd$rate, pd$town, pd$year, pd$first, adj = sep, est = "reg",
              reps = 0)
  expect_true(all(is.na(s$att$est[s$att$cohort == 2018])))
})

test_that("Wooldridge's cells are the coefficients of the cell dummies", {
  s <- stgdid(pd$rate, pd$town, pd$year, pd$first, meth = "wool", reps = 0)
  cell <- ifelse(pd$on == 1, paste(pd$first, pd$year), "0")
  m <- lm(rate ~ factor(cell) + factor(town) + factor(year), pd)
  b <- coef(m)[paste0("factor(cell)", s$att$cohort, " ", s$att$period)]
  expect_equal(s$att$est, unname(b), tolerance = 1e-10)
  # without never-treated units the periods from the last adoption go
  k <- !is.na(pd$first)
  s2 <- stgdid(pd$rate[k], pd$town[k], pd$year[k], pd$first[k],
               meth = "wool", reps = 0)
  expect_true(all(is.na(s2$att$est[s2$att$period >= 2022])))
  d <- pd[k & pd$year < 2022, ]
  cell <- ifelse(d$on == 1, paste(d$first, d$year), "0")
  m2 <- lm(rate ~ factor(cell) + factor(town) + factor(year), d)
  expect_equal(s2$att$est[1], unname(coef(m2)["factor(cell)2018 2018"]),
               tolerance = 1e-10)
})

test_that("DID_M and its placebo are the switchers' differences", {
  law <- pd$on * !(pd$town %in% c(1, 2) & pd$year >= 2022)
  s <- stgdid(pd$rate, pd$town, pd$year, trt = law, meth = "dcdh", reps = 0)
  # by hand, period by period
  yw <- matrix(pd$rate, 10)
  dw <- matrix(law, 10)
  num <- den <- pn <- pdn <- 0
  for (a in 2:10) {
    for (dir in 0:1) {
      st <- dw[a - 1, ] == dir
      sw <- st & dw[a, ] != dir
      cc <- st & dw[a, ] == dir
      if (any(sw) && any(cc)) {
        dd <- mean(yw[a, sw] - yw[a - 1, sw]) - mean(yw[a, cc] - yw[a - 1, cc])
        num <- num + sum(sw) * if (dir == 1) -dd else dd
        den <- den + sum(sw)
      }
      if (a < 3) next
      st <- st & dw[a - 2, ] == dir
      sw <- st & dw[a, ] != dir
      cc <- st & dw[a, ] == dir
      if (any(sw) && any(cc)) {
        dd <- mean(yw[a - 1, sw] - yw[a - 2, sw]) -
          mean(yw[a - 1, cc] - yw[a - 2, cc])
        pn <- pn + sum(sw) * if (dir == 1) -dd else dd
        pdn <- pdn + sum(sw)
      }
    }
  }
  expect_equal(s$overall$est, c(num / den, pn / pdn), tolerance = 1e-13)
  expect_identical(s$overall$n, c(den, pdn))
  expect_null(s$att)
  expect_output(print(s), "dcdh")
  # from first: the same as from the absorbing treatment
  a <- stgdid(pd$rate, pd$town, pd$year, pd$first, meth = "dcdh", reps = 0)
  b <- stgdid(pd$rate, pd$town, pd$year, trt = pd$on, meth = "dcdh",
              reps = 0)
  expect_equal(a$overall, b$overall, tolerance = testthat_tolerance())
})

test_that("the bootstrap is seeded, clustered, and its errors fit", {
  a <- stgdid(pd$rate, pd$town, pd$year, pd$first, reps = 30, seed = 4)
  b <- stgdid(pd$rate, pd$town, pd$year, pd$first, reps = 30, seed = 4)
  expect_identical(a, b)
  set.seed(9)
  c1 <- stgdid(pd$rate, pd$town, pd$year, pd$first, reps = 30)
  set.seed(9)
  c2 <- stgdid(pd$rate, pd$town, pd$year, pd$first, reps = 30)
  expect_identical(c1, c2)
  r <- stgdid(pd$rate, pd$town, pd$year, pd$first, clus = pd$region,
              reps = 20, seed = 1)
  expect_true(all(r$att$se > 0))
  z <- qnorm(0.975)
  expect_equal(r$overall$hi, r$overall$est + z * r$overall$se,
               tolerance = testthat_tolerance())
  # one sample: no spread
  r1 <- stgdid(pd$rate, pd$town, pd$year, pd$first, reps = 1, seed = 1)
  expect_true(is.na(r1$overall$se))
  # two clusters, many samples: some lack a cohort
  r2 <- stgdid(pd$rate, pd$town, pd$year, pd$first, meth = "wool",
               clus = as.integer(is.na(pd$first)), reps = 40, seed = 2)
  expect_true(is.finite(r2$overall$se))
})

test_that("the effect is recovered and the bootstrap error fits the spread", {
  set.seed(21)
  est <- se <- numeric(150)
  for (k in seq_along(est)) {
    d <- expand.grid(t = 1:6, u = 1:30)
    g <- c(3, 5, NA)[1 + (d$u - 1) %% 3]
    on <- !is.na(g) & d$t >= g
    y <- rnorm(30)[d$u] + 0.3 * d$t + 2 * on + rnorm(180)
    s <- stgdid(y, d$u, d$t, g, reps = 100, seed = k)
    est[k] <- s$overall$est
    se[k] <- s$overall$se
  }
  expect_lt(abs(mean(est) - 2), 4 * sd(est) / sqrt(150))
  # 150 draws: the sd's own error is about 6%; the bootstrap's about 7%
  expect_lt(abs(mean(se) / sd(est) - 1), 0.15)
})

test_that("stgdid() refuses what it cannot estimate", {
  y <- pd$rate
  u <- pd$town
  t <- pd$year
  f <- pd$first
  expect_error(stgdid(y[-1], u[-1], t[-1], f[-1]), "balanced")
  expect_error(stgdid(y, u, t), "one of `first` and `trt`")
  expect_error(stgdid(y, u, t, f, trt = pd$on), "one of `first` and `trt`")
  expect_error(stgdid(y, u, t, f + 0.5), "one of the periods")
  expect_error(stgdid(y, u, t, f, adj = pd$x, meth = "wool"), "csa")
  expect_error(stgdid(y, u, t, f, clus = t), "constant within")
  expect_error(stgdid(y, u, t, f, clus = rep(1, 300)), "one cluster")
  expect_error(stgdid(y, u, t, rep(2015, 300)), "after the first period")
  expect_error(stgdid(y, u, t, trt = rep(0, 300), meth = "dcdh"),
               "no unit switches")
  # first beyond the window is never; at the start, always (left out)
  f2 <- as.numeric(f)
  f2[u == 2] <- 2010
  f2[u == 1] <- 2030
  s <- stgdid(y, u, t, f2, reps = 0)
  expect_identical(s$att$ncoh[1], 5)
})
