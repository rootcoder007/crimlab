# Alert complexity and placement volatility (cxvol(), cxdml(), cxyoy()) on
# small made-up person-year tables and releases: every expected value is
# worked out by hand in the comments.

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

# A release of 200 person-years with one to three rows each.
rel <- function(n = 200) {
  set.seed(9)
  yr <- sample(2025:2026, n, TRUE)
  i <- rep(seq_len(n), sample(1:3, n, TRUE))
  m <- length(i)
  reg <- c("East", "West", "North")
  yn <- function(p) sample(c("Yes", "No"), m, TRUE, c(p, 1 - p))
  data.frame(
    EndFiscalYear = yr[i],
    UniqueIndividual_ID = paste0("p", i),
    Region_AtTimeOfPlacement = sample(reg, m, TRUE, c(6, 3, 1)),
    Region_MostRecentPlacement = sample(reg, m, TRUE, c(6, 3, 1)),
    Gender = sample(c("Male", "Female"), n, TRUE, c(0.8, 0.2))[i],
    Age_Category = sample(c("18 to 24", "25 to 49", "50+"), n, TRUE)[i],
    MentalHealth_Alert = yn(0.2),
    SuicideRisk_Alert = yn(0.1),
    SuicideWatch_Alert = "No",
    Number_Of_Placements = sample(1:4, m, TRUE)
  )
}

test_that("cxvol() and cxyoy() take the release, refuse gaps and print", {
  r <- rel()
  p <- otpyr(r)
  v <- cxvol(p, fam = "poisson")
  expect_identical(cxvol(r, fam = "poisson"), v)
  expect_identical(cxyoy(r), cxyoy(p))
  q <- pyfix()
  errs(
    "cxvol",
    list(
      list("none missing", replace(q, "vm", list(c(NA, q$vm[-1])))),
      list("no person-year with cx = 0", transform(q, cx = 1))
    )
  )
  prt(
    v,
    paste0(
      "^cxvol: alert complexity \\(cx\\) and placement volatility \\(vm\\), ",
      "years pooled\n"
    )
  )
  expect_output(
    print(v),
    "irr: rate ratio of the matched count model \\(poisson\\)"
  )
  prt(cxvol(p, yterm = FALSE, fam = "poisson"), "years pooled, no year term\n")
  prt(
    cxvol(p, by = "year", fam = "poisson"),
    "each fiscal year \\(no year term\\)\n.*\n *2025 .*\n *2026 "
  )
  prt(cxyoy(p), "^cxyoy: regional changes per fiscal year\n")
  prt(
    cxdml(q, lrn = "glm"),
    "^cxdml: interactive model, effect of cx on vm \\(additive scale\\)\n"
  )
})

test_that("a covariate with one value in a sample is left out", {
  set.seed(7)
  n <- 300
  p <- data.frame(
    year = 2025, gender = "Male", age = sample(c("a", "b"), n, TRUE),
    rc = sample(c("A", "B", "C"), n, TRUE), cx = rbinom(n, 1, 0.3),
    stringsAsFactors = FALSE
  )
  p$vm <- rpois(n, exp(-1 + 0.4 * p$cx))
  one <- cxvol(p, yterm = FALSE, fam = "poisson")
  # the same as with a gender that varies but is no use to the matching:
  # the estimate is a finite rate ratio
  est <- one$fits[[1L]]$coef
  expect_true(all(is.finite(est$est)))
  expect_false(any(grepl("gender", est$term)))
  dm <- cxdml(p, yterm = FALSE, ntree = 20, nfold = 2, seed = 1)
  expect_true(all(is.finite(dm$est)))
  p$age <- "a"
  errs("cxvol", list(
    list("one value of every covariate", p, yterm = FALSE, fam = "poisson")
  ))
})
