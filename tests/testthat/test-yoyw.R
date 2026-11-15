# The change-table writers (R/yoyw.R, R/yoyh.R). The files in ref/yoyw
# were written once by rmoriebricklayer 1.1.1 (yoycsv(digits = 2),
# yoytsv(digits = 2, meta = FALSE), yoymd(), yoyhtml() with its version
# line replaced by GENERATOR, and the cells of its PDF table) on the
# simulated table below; the same table, laid out as yoypx() lays it out,
# must give the same bytes.

d <- data.frame(prd = rep(2019:2022, 2),
                group = rep(c("Assault", "Robbery"), each = 4),
                value = c(812, 845, 790, 901, 120, 131, 98, 98))
y <- yoypx(d$value, d$prd, d$group)
rf <- function(f) test_path("ref", "yoyw", f)
bytes <- function(f) rawToChar(readBin(f, "raw", file.size(f)))
meta <- paste("value: value; period: prd; by: group; lag: 1; units: count;",
              "interval: 95% exact rate ratio")

# the table as rmoriebricklayer names its columns
bl <- function(y) {
  x <- as.data.frame(unclass(y))[c("group", "prd", "value", "prev", "change",
                                   "pctchg", "flag", "pctlo", "pctup")]
  names(x)[4L] <- "previous"
  x$flag[!is.na(x$flag)] <- "no comparison period"
  x$verdict <- c("down", "unchanged", "up")[sign(x$change) + 2L]
  x
}

test_that("the delimited writers give rmoriebricklayer's bytes", {
  expect_identical(yoycsv(bl(y), digits = 2, meta = meta), bytes(rf("bl.csv")))
  expect_identical(yoytsv(bl(y), digits = 2), bytes(rf("bl.tsv")))
})

test_that("yoymd() gives rmoriebricklayer's bytes", {
  z <- y
  z$flag[!is.na(z$flag)] <- "no comparison period"
  # bl.md is the file rmoriebricklayer wrote with writeLines(), which adds
  # one newline after the text its yoymd() returns; the text is the table
  ref <- sub("\n$", "", bytes(rf("bl.md")))
  expect_identical(yoymd(z, meta = meta), ref)
  f <- tempfile(fileext = ".md")
  on.exit(unlink(f))
  yoymd(z, f, meta = meta)
  expect_identical(bytes(f), ref)
})

test_that("yoyhtml() gives rmoriebricklayer's bytes", {
  f <- tempfile(fileext = ".html")
  on.exit(unlink(f))
  sub <- "lag 1 period \u00b7 units: count \u00b7 95% exact rate-ratio interval"
  yoyhtml(y, f, title = "Reported offences (simulated)", sub = sub)
  h <- sub("crimlab [0-9.]+", "GENERATOR", bytes(f))
  expect_identical(h, bytes(rf("bl.html")))
})

test_that("yoypdf() writes rmoriebricklayer's cells in a PDF of its own", {
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f))
  expect_identical(yoypdf(y, f, notes = "simulated"), f)
  p <- bytes(f)
  expect_match(p, "^%PDF-1.4\n")
  expect_match(p, "%%EOF\n$")
  for (s in unique(readLines(rf("bl.pdftext")))) {
    expect_true(grepl(paste0("(", s, ") Tj"), p, fixed = TRUE), info = s)
  }
  big <- yoypx(rep(c(10, 12), 60), seq_len(120))
  yoypdf(big, f, height = 4, better = "higher", pal = "safe")
  expect_match(bytes(f), "(page 3 of", fixed = TRUE)
  expect_error(yoypdf(y, f, width = 1), "from 2 to 200")
})

test_that("the summary equals rmoriebricklayer's", {
  s <- yoysmry(y)
  expect_identical(s$group, c("Assault", "Robbery"))
  expect_identical(s$from, c(2019, 2019))
  expect_equal(c(s$totalpct, s$cagrpct),
               c(10.960591133004915, -18.333333333333336, 3.5276256728592248,
                 -6.5279839465936291), tolerance = 1e-14)
  expect_identical(c(s$up, s$down, s$flat), c(2L, 1L, 1L, 1L, 0L, 1L))
  expect_output(print(s), "2 series")
  z <- yoypx(c(0, 5, 10), 1:3, units = "continuous")
  expect_true(is.na(yoysmry(z)$totalpct))
  expect_error(yoysmry(d), "must be a change table")
  expect_error(yoysmry(y, lag = 0), "one positive number")
})

test_that("labels change the periods shown, never the periods", {
  expect_identical(fylbl(c(2023, 2001, NA)), c("2022/23", "2000/01", NA))
  expect_identical(fylbl("2023", "-", FALSE), "2022-2023")
  expect_error(fylbl(2023.5), "whole years")
  expect_error(fylbl(-1), "positive years")
  expect_error(fylbl(TRUE), "numeric or character")
  z <- yoylbl(y, fylbl)
  expect_identical(z$prd, y$prd)
  expect_identical(z$perlbl[1:2], c("2018/19", "2019/20"))
  expect_match(yoymd(z), "2018/19", fixed = TRUE)
  expect_identical(yoylbl(y, c("2019" = "first"))$perlbl[1:2],
                   c("first", "2020"))
  expect_identical(yoylbl(y, letters[1:4])$perlbl[5], "a")
  expect_error(yoylbl(y, letters[1:3]), "3 entries for 4")
  expect_error(yoylbl(y, function(p) stop("no")), "failed on the periods")
  expect_error(yoylbl(y, list(1)), "function of the periods")
  expect_error(yoylbl(y, function(p) "a"), "one label per row")
})

test_that("a table written and read back is the table", {
  a <- data.frame(code = c("046", "", "112"), n = c(1 / 3, 2, 1e-20),
                  ok = c(TRUE, FALSE, NA), stringsAsFactors = FALSE)
  for (w in list(yoycsv, yoytsv)) {
    f <- tempfile()
    w(a, f, na = "NA")
    b <- utils::read.table(f, header = TRUE, sep = if (identical(w, yoycsv))
      "," else "\t", colClasses = c("character", "numeric", "logical"),
      na.strings = "NA", quote = "\"")
    unlink(f)
    expect_identical(b, a)
  }
  j <- yoyjs(a, pretty = FALSE)
  expect_match(j, "\"code\":\"046\"", fixed = TRUE)
  expect_match(j, "\"n\":0.3333333333333333", fixed = TRUE)
  expect_match(j, "\"ok\":null", fixed = TRUE)
  expect_match(yoyjs(y, meta = "simulated"), "\"meta\": \"simulated\"")
  expect_match(yoycsv(data.frame(a = "x,\"y\"")), "\"x,\"\"y\"\"\"",
               fixed = TRUE)
  expect_match(yoytsv(data.frame(a = "x\ty")), "x y", fixed = TRUE)
  expect_match(yoyjs(data.frame(a = "q\"\n\001", b = Inf)),
               "\"q\\\"\\n\\u0001\"", fixed = TRUE)
  expect_error(yoycsv(1), "must be a data frame")
  expect_error(yoycsv(a, digits = 20), "from 0 to 15")
  expect_error(yoycsv(a, meta = NA), "NULL or text")
})

test_that("yoywr() picks the writer from the file name", {
  for (e in c("html", "csv", "tsv", "json", "md", "pdf")) {
    f <- tempfile(fileext = paste0(".", e))
    expect_identical(yoywr(y, f), f)
    expect_gt(file.size(f), 100)
    unlink(f)
  }
  f <- tempfile(fileext = ".txt")
  expect_error(yoywr(y, f), "cannot tell the format")
  expect_identical(yoywr(y, f, format = "markdown"), f)
  unlink(f)
  expect_error(yoywr(y, f, format = "doc"), "unknown format")
  expect_error(yoyhtml(y, NA), "one file path")
  expect_error(yoyhtml(y, f, pal = "x"), "`pal` must be one of")
  expect_error(yoyhtml(y, f, notes = 1), "`notes` must be text")
  expect_match(yoymd(yoypx(c(5, 6), 1:2, units = "percent")), "points")
  expect_match(yoymd(y, align = FALSE), "| --- |", fixed = TRUE)
})

test_that("the colour follows the better direction; control bytes show", {
  f <- tempfile(fileext = ".html")
  on.exit(unlink(f))
  z <- y
  z$group[1] <- "A<b>\033"
  yoyhtml(z, f, better = "lower", notes = "n", bars = FALSE, intv = FALSE)
  h <- bytes(f)
  expect_match(h, "<tr class=\"bad\">", fixed = TRUE)
  expect_match(h, "<tr class=\"good\">", fixed = TRUE)
  expect_match(h, "A&lt;b&gt;\\033", fixed = TRUE)
  yoyhtml(z, f, better = "higher")
  expect_match(bytes(f), "<tr class=\"good\">", fixed = TRUE)
  # as rmoriebricklayer writes it: the control byte as the text \033, whose
  # backslash Markdown leaves alone (a backslash escapes only punctuation)
  expect_match(yoymd(z), "| A&lt;b&gt;\\033 |", fixed = TRUE)
  z$group[1] <- "046 & A\\*b"
  m <- yoymd(z)
  expect_match(m, "| 046 &amp; A\\\\\\*b |", fixed = TRUE)
  expect_match(m, "| 2019 |   812 |", fixed = TRUE)
  yoyhtml(z, f)
  expect_match(bytes(f), ">046 &amp; A\\*b<", fixed = TRUE)
})
