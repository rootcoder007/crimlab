# vrfpol() and polrpttxt(): the reports of the reference implementation
# (rmoriebricklayer 1.1.1, computed once on these simulated inputs and
# printed with 17 significant digits), its text, and the error reports.

vp_check <- function(r, inputs, rr, paf, pafci, dd, att, ci = NULL) {
  expect_s3_class(r, "vrfpol")
  expect_identical(r$status, "ok")
  expect_identical(attr(r, "exitstat"), 0L)
  expect_equal(unname(unlist(r$inputs)), inputs, tolerance = 1e-12)
  expect_true(all(vapply(r$assumptions, `[[`, NA, "ok")))
  p <- r$pipeline
  expect_equal(p$crf$rr, rr, tolerance = 1e-12)
  expect_equal(p$paf, paf, tolerance = 1e-12)
  expect_equal(unname(p$pafci), pafci, tolerance = 1e-12)
  expect_equal(p$displaced$deadisp, dd, tolerance = 1e-12)
  expect_equal(p$burden$attrcases, att, tolerance = 1e-12)
  if (is.null(ci)) {
    expect_null(p$equity)
  } else {
    expect_equal(p$equity$concidx, ci, tolerance = 1e-12)
  }
}

test_that("the demo data give the reference reports", {
  vp_check(
    vrfpol("no2", demo = TRUE),
    c(24.261137432351465, 0.98099999999999998, 500, 1000000, 10, 0),
    1.0286433508290656, 0.027331145821359856,
    c(0.013826085361389974, 0.053419353782373266), 136.58342874945922,
    136.65572910679927, -0.0646924587654783
  )
  vp_check(
    vrfpol("pm25", demo = TRUE),
    c(10.111016119469536, 0.90900000000000003, 500, 1000000,
      5.7999999999999998, 0),
    1.0337345574772916, 0.029752365019980683,
    c(0.022600530176668458, 0.03326134312880151), 148.32005240151668,
    148.7618250999034, -0.066174373400658293
  )
  set.seed(5)
  s <- .Random.seed
  vrfpol("no2", demo = TRUE)
  expect_identical(.Random.seed, s)
})

test_that("scalar inputs give the reference report and its text", {
  s <- vrfpol("pm25", expomean = 9, expoprev = 1, brate = 750, pop = 2.8e6,
              reg = "Toronto", years = "2019-2021")
  vp_check(s, c(9, 1, 750, 2800000, 5.7999999999999998, 0),
           1.0249332957580797, 0.024326749712661116,
           c(0.018473288441899549, 0.027200092442166115),
           510.86174396588382, 510.86174396588342)
  expect_identical(
    vapply(s$assumptions, `[[`, "", "note")[1:6],
    c("mean 9 vs ref 5.8", "exposure_mean=9", "exposure_prevalence=1",
      "brate=750 per 100k per year", "population=2,800,000",
      "reference=5.8")
  )
  want <- paste0(
      "==================================================================\n",
      "  verify-pollution -- PM25 -> allcause\n",
      "  region: Toronto\n",
      "  years:  2019-2021\n",
      "  data:   scalar arguments\n",
      "==================================================================\n",
      "\n",
      "Inputs\n",
      "  exposure mean:      9.000\n",
      "  exposure prevalence:1.000\n",
      "  baseline rate/100k: 750.00\n",
      "  population:         2,800,000\n",
      "  reference conc:     5.8\n",
      "\n",
      "Assumption log\n",
      "  [PASS] exposure finite and non-negative -- mean 9 vs ref 5.8\n",
      "  [PASS] exposure within the concentration-response range (at most",
      " 10,000 ug/m3) -- exposure_mean=9\n",
      "  [PASS] prevalence in [0,1] -- exposure_prevalence=1\n",
      "  [PASS] brate finite: zero, or between 1 and 10,000 per 100k --",
      " brate=750 per 100k per year\n",
      "  [PASS] population finite, positive and at most 10 billion --",
      " population=2,800,000\n",
      "  [PASS] reference finite and non-negative -- reference=5.8\n",
      "  [PASS] pollutant supported by the CRFs -- Current CRFs: NO2",
      " (log-linear), PM2.5 (log-linear all-cause; Burnett IER for IHD and",
      " stroke). Other pollutants reject.\n",
      "\n",
      "Concentration-response\n",
      "  RR:       1.0249\n",
      "  source:   Chen & Hoek (2020) Environ Int 143:105974; WHO (2021)",
      " Global AQ Guidelines\n",
      "\n",
      "Attributable fraction (PAF): 0.0243  (95% CI 0.0185 to 0.0272, from",
      " the risk's limits)\n",
      "\n",
      "Mortality displaced\n",
      "  expected avoided deaths: 510.9\n",
      "\n",
      "Burden of pollution\n",
      "  attributable deaths:   510.9\n",
      "\n",
      "STATUS: ok\n"
    )
  expect_identical(polrpttxt(s), want)
  expect_output(print(s), "STATUS: ok")
})

test_that("files are read as the reference reads them", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(exposure = c(12, 18, 25, 31, 9), income = 1:5),
                   f, row.names = FALSE)
  vp_check(
    vrfpol("no2", expcsv = f),
    c(21.5, 0.80000000000000004, 500, 1000000, 10, 0),
    1.0230343062880205, 0.018094018499144087,
    c(0.0091228872072107604, 0.03559537561195459), 90.062693485219555,
    90.470092495720436, 0.029473684210526308
  )
  g <- tempfile(fileext = ".csv.gz")
  con <- gzfile(g, "w")
  writeLines(c("value,unit", "10,ppb", "20,ppb", ",ppb", "15,ug/m3"), con)
  close(con)
  expect_message(
    r <- vrfpol("no2", expcsv = g),
    "2 NO2 row\\(s\\) in ppb converted"
  )
  vp_check(
    r, c(23.799999999999997, 1, 500, 1000000, 10, 1),
    1.0277044499710215, 0.026957604369430042,
    c(0.013637610167880405, 0.052685921170979449), 134.78802184715033,
    134.78802184715022
  )
  h <- tempfile(fileext = ".csv")
  writeLines(c("exposure", "12", "<2", "30"), h)
  expect_identical(
    vrfpol("no2", expcsv = h)$error,
    paste(
      "1 exposure value(s) are text, not numbers (\"<2\"): a censored value",
      "such as \"<2\" is a measurement, and dropping it would bias the mean",
      "and the prevalence; replace them (for example with half the detection",
      "limit) or remove those rows"
    )
  )
  expect_match(vrfpol("no2", expcsv = tempdir())$error, "is a directory")
})

test_that("failed assumptions and unusable arguments are reports", {
  b <- vrfpol("no2", expomean = 20, expoprev = 1, brate = 0.008)
  expect_identical(b$status, "assumption_failure")
  expect_identical(attr(b, "exitstat"), 1L)
  expect_identical(vapply(b$assumptions, `[[`, NA, "ok"),
                   c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_match(polrpttxt(b),
               "\\[FAIL\\] brate finite.*STATUS: assumption_failure")
  e <- vrfpol("no2")
  expect_identical(
    unclass(e)[c("status", "error")],
    list(status = "error",
         error = paste("no exposure data: pass expcsv, or expomean and",
                       "expoprev, or demo = TRUE"))
  )
  expect_identical(attr(e, "exitstat"), 2L)
  expect_identical(polrpttxt(e), paste0("ERROR: ", e$error, "\n"))
  expect_identical(
    vrfpol("no2", expomean = 20)$error,
    paste("expomean and expoprev go together: pass both (a zero is a",
          "measurement), or neither")
  )
  z <- vrfpol("no2", expomean = 8, expoprev = 0.5)
  expect_identical(z$pipeline$paf, 0)
  expect_identical(z$pipeline$burden$attrcases, 0)
  expect_error(vrfpol(3), "^vrfpol: `poll` must be one pollutant name")
  expect_error(polrpttxt(list(1)),
               "^polrpttxt: `report` must be a report from vrfpol\\(\\)\\.$")
  expect_error(polrpttxt(list(status = "ok")), "not a complete vrfpol")
})
