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

test_that("each unusable argument is its own error report", {
  vp_err <- function(...) vrfpol("no2", ...)$error
  s <- list(expomean = 20, expoprev = 1)
  vs <- function(...) do.call(vp_err, c(s, list(...)))
  expect_identical(vp_err(demo = "yes"), "demo must be TRUE or FALSE")
  expect_identical(
    vp_err(demo = TRUE, expcsv = "a.csv"),
    paste("expcsv and demo = TRUE are two sources of exposure: pass one",
          "(the demo was used silently before)")
  )
  expect_identical(
    vp_err(demo = TRUE, expomean = 1),
    paste("expomean / expoprev cannot be combined with expcsv or demo:",
          "pass one source of exposure")
  )
  expect_identical(
    vs(pop = NULL),
    "population must be a single number (persons at risk), not NULL"
  )
  expect_identical(vs(brate = "5"),
                   "brate must be a single number, not character")
  expect_identical(vs(brate = NULL), "brate must be a single number, not NULL")
  expect_identical(vs(rare = NA), "rare must be TRUE or FALSE")
  expect_identical(vs(hirate = 1), "hirate must be TRUE or FALSE")
  expect_identical(
    vs(basrateper = 0),
    paste("basrateper must be a single positive number",
          "(100000 for a rate per 100,000, 1 per person-year)")
  )
  expect_identical(vs(brate = c(1, 2)),
                   "brate must be a single number, not 2 values")
  expect_identical(vs(y = 3), "outcome must be a single outcome name")
  # an outcome the curve does not have: the stage's refusal, as a report
  r <- vrfpol("no2", expomean = 20, expoprev = 1, y = "asthma")
  expect_identical(attr(r, "exitstat"), 2L)
  expect_identical(
    r$error,
    paste("crfno2: Unknown outcome 'asthma'. Available: allcause,",
          "respiratory, chiast or pass betaper10 explicitly.")
  )
  expect_identical(vp_err(expcsv = 3), "expcsv must be a single file path")
  f <- file.path(tempdir(), "vrfp-none.csv")
  expect_identical(vp_err(expcsv = f),
                   sprintf("exposure CSV not found: %s", f))
})

test_that("each assumption that fails names why in its note", {
  note <- function(r, i) {
    expect_identical(r$status, "assumption_failure")
    expect_false(r$assumptions[[i]]$ok)
    r$assumptions[[i]]$note
  }
  o <- vrfpol("o3", expomean = 20, expoprev = 1)
  expect_identical(o$pollutant, "o3")
  expect_identical(vapply(o$assumptions, `[[`, NA, "ok"),
                   c(rep(TRUE, 6), FALSE))
  expect_identical(
    note(vrfpol("no2", expomean = 20, expoprev = 1, pop = 2e10), 5),
    "population=20,000,000,000 is more people than live on Earth"
  )
  expect_identical(
    note(vrfpol("no2", expomean = -1, expoprev = 1), 1),
    "expomean=-1 is not a finite non-negative concentration"
  )
  tiny <- vrfpol("no2", expomean = 20, expoprev = 1, brate = 1e-7,
                 rare = TRUE)
  expect_identical(
    tiny$assumptions[[4]]$assumption,
    "brate finite: zero, or between 1e-6 and 10,000 per 100k (a rare outcome)"
  )
  expect_identical(
    note(tiny, 4),
    paste("brate=1e-07 per 100k is below one case per 100 billion",
          "person-years: not a usable rate (check basrateper)")
  )
  expect_identical(
    note(vrfpol("no2", expomean = 20, expoprev = 1, brate = 2e4), 4),
    paste("brate=20,000 per 100k a year is above a whole population's",
          "all-cause death rate (an age-specific rate for the oldest ages",
          "can be this high): if the outcome really is this common, say so",
          "with hirate = TRUE")
  )
  # a measured zero rate passes and gives no cases
  z <- vrfpol("no2", expomean = 20, expoprev = 1, brate = 0)
  expect_identical(z$status, "ok")
  expect_identical(z$assumptions[[4]]$note,
                   "brate=0: a measured zero rate, so no cases")
  expect_identical(z$pipeline$burden$attrcases, 0)
})

test_that("the text names cases for asthma and leaves out a missing CI", {
  # the IHD curve (Burnett IER) has no published interval: no CI is printed
  i <- vrfpol("pm25", y = "ihd", expomean = 12, expoprev = 1)
  expect_null(i$pipeline$pafci)
  # the IER curve: alpha 1.91, gamma 0.14, delta 0.49, 6.2 above the floor
  rr <- 1 + 1.91 * (1 - exp(-0.14 * 6.2^0.49))
  expect_equal(i$pipeline$crf$rr, rr, tolerance = 1e-12)
  expect_match(polrpttxt(i),
               "\nAttributable fraction (PAF): 0.3563\n\nMortality",
               fixed = TRUE)
  # childhood asthma counts incident cases: 1.05^(2.5 * 20 / 10) = 1.05^5
  ch <- vrfpol("no2", y = "chiast", expomean = 30, expoprev = 1)
  expect_equal(ch$pipeline$crf$rr, 1.05^5, tolerance = 1e-12)
  expect_identical(ch$pipeline$displaced$unit, "incident cases")
  txt <- polrpttxt(ch)
  expect_match(
    txt,
    "\nCases displaced\n  expected avoided incident cases: 1082.4\n",
    fixed = TRUE
  )
  expect_match(txt, "  attributable incident cases:   1082.4\n", fixed = TRUE)
  # income data add the equity section
  expect_match(
    polrpttxt(vrfpol("no2", demo = TRUE)),
    "\n\nEquity analysis\n  concentration index: -0.0647\n\nSTATUS: ok\n",
    fixed = TRUE
  )
})

test_that("rows and columns of a CSV are checked before they are used", {
  fs <- character()
  on.exit(unlink(fs), add = TRUE)
  w <- function(lines) {
    f <- tempfile(fileext = ".csv")
    fs <<- c(fs, f)
    writeLines(lines, f)
    f
  }
  ve <- function(f) vrfpol("no2", expcsv = f)$error
  expect_identical(
    ve(w(c("pm", "1"))),
    paste0("CSV missing 'exposure' column (ug/m3); a NAPS pull's 'value' ",
           "column also works.")
  )
  expect_identical(
    ve(w(c("value", "1"))),
    paste0("a 'value' column needs a 'unit' column (NAPS pulls carry ",
           "one; NO2 is in ppb)")
  )
  expect_identical(ve(w(c("value,unit", "1,ppm"))),
                   "exposure unit ppm is not ug/m3 for no2")
  # a missing-value token is dropped, not refused, and the text says so
  expect_message(
    r <- vrfpol("no2", expcsv = w(c("exposure", "12", "n/a", "30"))),
    "note: 1 row(s) have no exposure value and were dropped", fixed = TRUE
  )
  expect_identical(r$inputs$expomean, 21)
  expect_identical(r$inputs$drorows, 1L)
  expect_match(polrpttxt(r),
               "\n  rows dropped:       1 (no exposure value)\n",
               fixed = TRUE)
  expect_message(
    r <- vrfpol("no2", expcsv = w(c("exposure,income", "12,", "30,"))),
    "note: the income column is empty: no equity analysis", fixed = TRUE
  )
  expect_identical(r$status, "ok")
  expect_null(r$pipeline$equity)
  expect_identical(
    ve(w(c("exposure,income", "12,a", "30,b"))),
    paste0("the income column must be numeric (a quintile or an amount); ",
           "text cannot be ranked")
  )
  # one complete row, or no exposure above zero: no index, but a report
  r <- vrfpol("no2", expcsv = w(c("exposure,income", "12,3", "30,NA")))
  expect_identical(r$status, "ok")
  expect_null(r$pipeline$equity)
  r <- vrfpol("no2", expcsv = w(c("exposure,income", "0,1", "0,2")))
  expect_identical(r$status, "ok")
  expect_null(r$pipeline$equity)
  # one income level: the index is undefined and the stage refuses
  expect_identical(
    ve(w(c("exposure,income", "12,3", "30,3"))),
    paste("expconcidx: `income` must take at least two values: with one",
          "income level there is no income gradient and the index is",
          "undefined.")
  )
  # nothing above the reference: the mean of all rows, prevalence 0
  r <- vrfpol("no2", expcsv = w(c("exposure,income", "4,1", "6,2")))
  expect_identical(unlist(r$inputs[c("expomean", "expoprev")]),
                   c(expomean = 5, expoprev = 0))
  expect_identical(r$pipeline$paf, 0)
  f <- tempfile(fileext = ".csv")
  fs <- c(fs, f)
  file.create(f)
  emp <- function(f) {
    sprintf(paste0("%s is empty (expected a column exposure in ug/m3, or a ",
                   "NAPS pull's value + unit)"), f)
  }
  expect_identical(ve(f), emp(f))
  f <- w(c("  ", ""))
  expect_identical(ve(f), emp(f))
  f <- w(c("exposure,Exposure", "1,2"))
  expect_identical(
    ve(f),
    paste0(f, ": has the column name 'Exposure' more than once: say which ",
           "column is meant")
  )
  f <- w(c("exposure,income", "1,2", "3"))
  expect_identical(
    ve(f),
    paste0(f, ": line 3 has 1 fields where the header has 2: no rows are ",
           "used from a file whose rows do not line up")
  )
  # a quote left open past the first lines: the parser stops part-way
  f <- w(c("exposure,income", paste0(11:18, ",", 1:8), "19,\"9", "20,10"))
  expect_identical(
    ve(f),
    paste0(f, ": could not be read cleanly (EOF within quoted string): no ",
           "rows are used from a file the parser stopped part-way through")
  )
  # no newline at the end is not a problem
  f <- tempfile(fileext = ".csv")
  fs <- c(fs, f)
  writeChar("exposure\n12\n30", f, eos = NULL)
  expect_identical(vrfpol("no2", expcsv = f)$inputs$expomean, 21)
})

test_that("BUG: a quote left open in the first rows is not an empty file", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("exposure,income", "1,2", "3,\"4", "5,6"), f)
  expect_match(vrfpol("no2", expcsv = f)$error,
               "could not be read cleanly", fixed = TRUE)
})

test_that("BUG: the negative-value refusal is a whole sentence", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("exposure", "12", "-3", "30"), f)
  expect_identical(
    vrfpol("no2", expcsv = f)$error,
    "1 exposure value(s) are negative: a concentration cannot be negative"
  )
})

test_that("compressed files are read whole or refused", {
  fs <- character()
  on.exit(unlink(fs), add = TRUE)
  tf <- function(ext = ".csv") {
    f <- tempfile(fileext = ext)
    fs <<- c(fs, f)
    f
  }
  wb <- function(raw) {
    f <- tf()
    writeBin(raw, f)
    f
  }
  rd <- function(f) readBin(f, "raw", file.size(f))
  ve <- function(f) vrfpol("no2", expcsv = f)$error
  lines <- c("exposure", "12", "30")
  for (open in list(bzfile, xzfile)) {
    f <- tf(".z")
    con <- open(f, "w")
    writeLines(lines, con)
    close(con)
    expect_identical(vrfpol("no2", expcsv = f)$inputs$expomean, 21)
  }
  # the xz file above, one byte damaged: the decoder's words are kept
  x <- rd(f)
  x[30] <- as.raw(255)
  f <- wb(x)
  expect_identical(
    ve(f),
    paste0(f, ": is truncated or damaged (the decoder said: lzma decoder ",
           "corrupt data)")
  )
  # a bzip2 file without its end-of-stream marker, long and short
  b <- tf(".bz2")
  con <- bzfile(b, "w")
  writeLines(lines, con)
  close(con)
  b <- rd(b)
  trunc <- ": is truncated: the bzip2 stream has no end-of-stream marker"
  f <- wb(b[seq_len(length(b) - 4L)])
  expect_identical(ve(f), paste0(f, trunc))
  f <- wb(charToRaw("BZh91"))
  expect_identical(ve(f), paste0(f, trunc))
  # a gzip file whose eight-byte trailer is cut
  g <- tf(".gz")
  con <- gzfile(g, "w")
  writeLines(c("exposure", 1:2000), con)
  close(con)
  g <- rd(g)
  f <- wb(g[seq_len(length(g) - 4L)])
  expect_identical(
    ve(f),
    paste0(f, ": is truncated or damaged (the decoder said: invalid or ",
           "incomplete compressed data)")
  )
  # a gzip file that holds nothing
  f <- tf(".gz")
  close(gzfile(f, "w"))
  expect_match(ve(f), "is empty (expected", fixed = TRUE)
})

test_that("BUG: a gzip file cut inside its data is refused", {
  f <- tempfile(fileext = ".gz")
  on.exit(unlink(f), add = TRUE)
  con <- gzfile(f, "w")
  writeLines(c("exposure", 1:2000), con)
  close(con)
  g <- readBin(f, "raw", file.size(f))
  writeBin(g[seq_len(length(g) - 40L)], f)
  # read silently, the file gives 1990 rows or so of the 2000
  r <- vrfpol("no2", expcsv = f)
  expect_identical(r$status, "error")
  expect_match(paste(r$error), "is truncated or damaged", fixed = TRUE)
})

test_that("a file that decompresses past the limit is refused", {
  # 68 MiB of blanks in about 70 KB of gzip: the limit is then 64 MiB
  f <- tempfile(fileext = ".gz")
  on.exit(unlink(f), add = TRUE)
  con <- gzfile(f, "wb")
  blank <- charToRaw(strrep(" ", 4194304L))
  for (i in 1:17) writeBin(blank, con)
  close(con)
  expect_lt(200 * file.size(f), 64 * 2^20)
  expect_identical(
    vrfpol("no2", expcsv = f)$error,
    paste0(f, ": decompresses to more than the size limit")
  )
})

test_that("UTF-16, UTF-32, a UTF-8 mark and Latin-1 are decoded", {
  fs <- character()
  on.exit(unlink(fs), add = TRUE)
  wb <- function(raw) {
    f <- tempfile(fileext = ".csv")
    fs <<- c(fs, f)
    writeBin(raw, f)
    f
  }
  txt <- "exposure\n12\n30\n"
  enc <- function(to) iconv(txt, "UTF-8", to, toRaw = TRUE)[[1L]]
  mean21 <- function(f) {
    expect_identical(vrfpol("no2", expcsv = f)$inputs$expomean, 21)
  }
  mean21(wb(c(as.raw(c(0xff, 0xfe)), enc("UTF-16LE"))))
  mean21(wb(c(as.raw(c(0xfe, 0xff)), enc("UTF-16BE"))))
  mean21(wb(c(as.raw(c(0xff, 0xfe, 0, 0)), enc("UTF-32LE"))))
  mean21(wb(c(as.raw(c(0xef, 0xbb, 0xbf)), charToRaw(txt))))
  for (to in c("UTF-16LE", "UTF-16BE")) {
    f <- wb(enc(to))
    expect_message(
      mean21(f),
      sprintf("note: %s has no byte-order mark; read as %s", basename(f), to),
      fixed = TRUE
    )
  }
  # a UTF-32 mark before three bytes that are no UTF-32 character
  f <- wb(c(as.raw(c(0xff, 0xfe, 0, 0)), charToRaw("ABC")))
  expect_identical(vrfpol("no2", expcsv = f)$error,
                   paste0(f, ": cannot be decoded as UTF-32LE"))
  f <- wb(c(charToRaw("exposure\n1"), as.raw(0), charToRaw("2\n")))
  expect_identical(
    vrfpol("no2", expcsv = f)$error,
    paste0(f, ": contains NUL bytes and is not UTF-16 text: not a text CSV")
  )
  # "caf\xe9" in Latin-1 is the header cafe with an acute accent
  f <- wb(c(charToRaw("exposure,caf"), as.raw(0xe9),
            charToRaw("\n12,1\n30,2\n")))
  expect_message(
    r <- vrfpol("no2", expcsv = f),
    sprintf("note: %s is not UTF-8; read as Latin-1 (ISO-8859-1)",
            basename(f)),
    fixed = TRUE
  )
  expect_identical(r$inputs$expomean, 21)
})
