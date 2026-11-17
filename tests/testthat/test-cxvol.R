# Placement volatility, alert complexity and the analyses built on them,
# on small made-up tables with the columns of the restrictive confinement
# release: every expected value is worked out by hand in the comments.
# Nothing here reaches the network (otpyr() by key reads a copy of the
# service in a temporary folder, helper-od.R).

a01 <- function() {
  data.frame(
    EndFiscalYear = c(2025, 2025, 2025, 2025, 2026, 2026, 2026),
    UniqueIndividual_ID = c("a", "a", "a", "b", "a", "c", "c"),
    Region_AtTimeOfPlacement = c(
      "West",
      "East",
      "West",
      "East",
      "North",
      "East",
      "East"
    ),
    Region_MostRecentPlacement = c(
      "West",
      "West",
      "West",
      "East",
      "East",
      "East",
      "East"
    ),
    Gender = c("Male", "Male", "Male", "Female", "Male", "Male", "Male"),
    Age_Category = "25 to 49",
    MentalHealth_Alert = c("Yes", "No", "No", "No", "No", "Yes", "Yes"),
    SuicideRisk_Alert = c("No", "No", "No", "No", "No", "Yes", "Yes"),
    SuicideWatch_Alert = c("No", "No", "No", "No", "No", "Yes", "Yes"),
    Number_Of_Placements = c(3, 1, 2, 5, 4, 1, 1)
  )
}

test_that("plcvol() counts changes within and between records", {
  # person 1: rows (W,W), (E,W), (W,W). Within: (E,W) = 1. Given order
  # W, E, W: two changes; sorted by region E, W, W: one.
  ra <- c("W", "E", "W", "N", "N")
  rb <- c("W", "W", "W", "N", "N")
  g <- plcvol(ra, rb, c(1, 1, 1, 2, 2), ord = "given")
  expect_identical(g$within, c(1L, 0L))
  expect_identical(g$across, c(2L, 0L))
  expect_identical(g$vm, c(3L, 0L))
  r <- plcvol(ra, rb, c(1, 1, 1, 2, 2))
  expect_identical(r$vm, c(2L, 0L))
  # sorted by region, across is the distinct regions at placement less one
  expect_identical(r$across, c(1L, 0L))
  # the bound 2R - 1 is reached by alternating records
  b <- plcvol(c("A", "B", "A"), c("B", "A", "B"), rep(1, 3), ord = "given")
  expect_identical(b$vm, 5L)
  # an injective relabelling changes nothing; merging regions can
  m <- c(W = "x", E = "y", N = "z")
  expect_identical(plcvol(m[ra], m[rb], c(1, 1, 1, 2, 2))$vm, r$vm)
  expect_identical(
    plcvol(rep("x", 5), rep("x", 5), c(1, 1, 1, 2, 2))$vm,
    c(0L, 0L)
  )
  # rows in another order: within unchanged; region order unchanged
  o <- c(3, 1, 2, 5, 4)
  expect_identical(plcvol(ra[o], rb[o], c(1, 1, 1, 2, 2)[o])$vm, r$vm)
  # a list id: person and year
  py <- plcvol(ra, rb, list(id = c(1, 1, 1, 1, 1), yr = c(1, 1, 1, 2, 2)))
  expect_identical(names(py)[1:2], c("id", "yr"))
  expect_identical(py$vm, c(2L, 0L))
  errs(
    "plcvol",
    list(
      list("one value per record", ra, rb[-1], 1:5),
      list("one label per record", ra, rb, 1:2),
      list("all present", NA, NA, 1),
      list("must be one of", ra, rb, 1:5, ord = "x")
    )
  )
})

test_that("alrtcx() counts distinct configurations, not alerts at once", {
  f <- data.frame(
    mh = c("Yes", "No", "No", "Yes", "Yes", NA),
    sr = c("No", "No", "Yes", "Yes", "Yes", "No"),
    sw = c("No", "No", "No", "Yes", "Yes", "No")
  )
  a <- alrtcx(f, c(1, 1, 2, 3, 3, 4))
  # 1: 100 and 000, two configurations, though never two alerts at once;
  # 2: 010 alone; 3: 111 twice, one configuration though three alerts;
  # 4: a missing flag, dropped
  expect_identical(a$id, c(1, 2, 3))
  expect_identical(a$ac, c(2L, 1L, 1L))
  expect_identical(a$cx, c(1L, 0L, 0L))
  expect_identical(a$anyal, c(1L, 1L, 1L))
  # keep counts only the configurations given
  expect_identical(
    alrtcx(f, c(1, 1, 2, 3, 3, 4), keep = "100")$ac,
    c(1L, 0L, 0L)
  )
  # logical and 0/1 flags give the same; order of rows does not matter
  l <- data.frame(x = c(TRUE, FALSE, TRUE), y = c(0, 0, 1))
  expect_identical(alrtcx(l, c(1, 1, 1))$ac, 3L)
  expect_identical(alrtcx(l[3:1, ], c(1, 1, 1))$ac, 3L)
  errs(
    "alrtcx",
    list(
      list("equal length", list(a = 1:2, b = 1), 1:2),
      list("codes of 2 digits", l, 1:3, keep = "1"),
      list("every flag present", data.frame(x = NA), 1)
    )
  )
})

test_that("otpyr() builds the person-year table from the release", {
  p <- otpyr(a01())
  # a 2025: rows W/W, E/W, W/W; sorted E, W, W: within 1, across 1;
  # configurations 100, 000: ac 2; placements 6; cluster "East + West"
  expect_identical(as.character(p$id), c("a", "b", "a", "c"))
  expect_identical(p$vm, c(2L, 0L, 1L, 0L))
  expect_identical(p$cx, c(1L, 0L, 0L, 0L))
  expect_identical(p$np, c(6, 5, 4, 2))
  expect_identical(p$rc, c("East + West", "East", "East + North", "East"))
  expect_identical(attr(p, "years"), c(2025, 2026))
  # Ontario's own spelling of the columns is found
  d <- a01()
  names(d) <- c(
    "end_fiscal_year",
    "unique_individual_id",
    "region_at_time_of_placement",
    "region_most_recent_placement",
    "gender",
    "age_category",
    "mental_health_alert",
    "suicide_risk_alert",
    "suicide_watch_alert",
    "number_of_placements"
  )
  expect_identical(otpyr(d)$vm, p$vm)
  # the years kept; a missing value dropped and counted
  expect_identical(nrow(otpyr(a01(), years = 2026)), 2L)
  d <- a01()
  d$Gender[2] <- NA
  expect_identical(attr(otpyr(d), "drop"), 1L)
  errs(
    "otpyr",
    list(
      list("no column", a01()[-1]),
      list("must be named", a01(), cols = c("x")),
      list("data frame or a table key", 1:3)
    )
  )
})

test_that("otpyr() reads a table of the service by its key", {
  dd <- odfix()
  gz <- file.path(dd, "db1", "a01.csv.gz")
  con <- gzfile(gz, "w")
  utils::write.csv(a01(), con, row.names = FALSE)
  close(con)
  odman(
    dd,
    sprintf(
      r"({"datasets": [{"db": "db1", "table": "a01", "key": "otis/a01",
    "rows": 7, "columns": [], "bytes_gz": %.0f, "sha256": "%s"}]})",
      file.size(gz),
      .odhash(gz, "test")
    )
  )
  odrun(dd, expect_identical(otpyr("otis/a01")$vm, otpyr(a01())$vm))
})

test_that("cntmix() finds the maximum of the marginal likelihood", {
  set.seed(3)
  g <- rep(1:4, each = 6)
  d <- rep(0:1, 12)
  y <- rpois(24, exp(-0.3 + 0.5 * d + c(-0.6, 0, 0.4, 0.8)[g]))
  f <- cntmix(y, data.frame(d = d), g, fam = "poisson")
  # the exact marginal likelihood by quadrature, at (b0, b1, sigma)
  ex <- function(b0, b1, s) {
    sum(vapply(
      1:4,
      function(j) {
        i <- g == j
        log(
          stats::integrate(
            function(u) {
              vapply(
                u,
                function(v) {
                  exp(sum(stats::dpois(
                    y[i],
                    exp(b0 + b1 * d[i] + v),
                    log = TRUE
                  )))
                },
                0
              ) *
                stats::dnorm(u, 0, s)
            },
            -8 * s,
            8 * s
          )$value
        )
      },
      0
    ))
  }
  b <- f$coef$est
  l0 <- ex(b[1], b[2], f$sigma)
  # the Laplace value is close to the exact one at the estimate
  expect_lt(abs(f$loglik - l0), 0.05)
  # and moving the treatment coefficient lowers the exact likelihood
  expect_lt(ex(b[1], b[2] + 0.1, f$sigma), l0)
  expect_lt(ex(b[1], b[2] - 0.1, f$sigma), l0)
  # relabelling the clusters changes nothing
  f2 <- cntmix(y, data.frame(d = d), letters[5 - g], fam = "poisson")
  expect_equal(f2$coef$est, f$coef$est, tolerance = 1e-6)
  errs(
    "cntmix",
    list(
      list("two clusters or more", y, data.frame(d = d), rep(1, 24)),
      list("0 throughout", 0 * y, data.frame(d = d), g),
      list("one complete row", y, data.frame(d = d[-1]), g),
      list("must be one of", y, data.frame(d = d), g, fam = "x")
    )
  )
})

pyfix <- function(years = 2025:2026, n = 600) {
  set.seed(7)
  p <- data.frame(
    year = sample(years, n, TRUE),
    gender = sample(c("Male", "Female"), n, TRUE, c(0.7, 0.3)),
    age = sample(c("18 to 24", "25 to 49", "50+"), n, TRUE),
    rc = sample(c("A", "B", "A + B", "C"), n, TRUE),
    cx = stats::rbinom(n, 1, 0.3)
  )
  p$vm <- stats::rpois(n, exp(-1 + 0.4 * p$cx + (p$rc == "A + B")))
  p
}

test_that("cxvol() gives the descriptives, matching, model and bounds", {
  p <- pyfix()
  r <- cxvol(p, fam = "poisson")$res
  expect_identical(r$n, 600L)
  expect_equal(r$vm1, mean(p$vm[p$cx == 1]))
  expect_equal(r$vm0, mean(p$vm[p$cx == 0]))
  # one control per complex person-year
  expect_identical(r$nm, 2L * sum(p$cx == 1))
  # the crude ratio is the ratio of the group means
  expect_equal(r$crude, r$vm1 / r$vm0, tolerance = 1e-8)
  # under both monotonicity assumptions the bound is [0, vm1m - vm0m]
  expect_equal(r$bhi, r$vm1m - r$vm0m, tolerance = 1e-12)
  # the Manski interval has width the outcome's range on the matched frame
  expect_equal(r$mhi - r$mlo, r$ymax, tolerance = 1e-12)
  expect_lte(r$ymax, max(p$vm))
  y <- cxvol(p, by = "year", fam = "poisson")$res
  expect_identical(y$sample, c("2025", "2026"))
  errs(
    "cxvol",
    list(
      list("cannot identify a year term", p[p$year == 2026, ]),
      list("no person-year with cx = 1", transform(p, cx = 0)),
      list("must be one of", p, by = "x")
    )
  )
  expect_silent(cxvol(p[p$year == 2026, ], yterm = FALSE, fam = "poisson"))
})

test_that("cxdml() gives the ATE and ATTE from one cross-fit", {
  p <- pyfix()
  r <- cxdml(p, lrn = "glm")
  expect_identical(r$estimand, c("ATE", "ATTE"))
  expect_true(all(is.finite(r$est) & r$se > 0))
  expect_identical(nrow(cxdml(p, by = "year", lrn = "glm")), 4L)
  errs(
    "cxdml",
    list(
      list("cannot identify a year term", p[p$year == 2025, ]),
      list("too few", transform(p[1:12, ], cx = c(1, rep(0, 11))), by = "year")
    )
  )
})

test_that("cxyoy() gives the rates and their changes", {
  p <- data.frame(
    year = c(1, 1, 2, 2, 2),
    vm = c(1, 1, 0, 3, 0),
    np = c(10, 10, 20, 30, 50)
  )
  r <- cxyoy(p)
  # year 1: 2 changes, 20 placements, 2 person-years; year 2: 3, 100, 3
  expect_equal(r$rnp, c(100, 30))
  expect_equal(r$rpy, c(1000, 1000))
  expect_equal(r$share, c(0.4, 0.6))
  expect_equal(r$dvm, c(NA, 50))
  expect_equal(r$drnp, c(NA, -70))
  expect_equal(r$rnplo[1], 1000 * stats::qgamma(0.025, 2) / 20)
  expect_identical(cxyoy(transform(p, vm = 0))$rnplo, c(0, 0))
})

test_that("cxsim() reports bias, spread and coverage", {
  s <- cxsim(reps = 2, n = 300, ncl = 6, seed = 1)
  expect_identical(s$res$est, c("crude", "mixed", "matched"))
  expect_true(all(s$res$ok <= 2))
  expect_true(is.na(s$res$cover[1]))
  s2 <- cxsim(reps = 2, n = 300, ncl = 6, conf = 1, seed = 1)
  expect_identical(s2$conf, 1)
  errs("cxsim", list(list("at least 1", reps = 0)))
})
