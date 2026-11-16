# The change-table renderers of R/yoyh.R: the HTML page, the PDF the
# package writes itself, and the writer that picks the format.

y <- yoypx(c(812, 845, 790), 2021:2023)

test_that("the HTML page escapes its text and draws the direction", {
  f <- tempfile(fileext = ".html")
  on.exit(unlink(f), add = TRUE)
  expect_identical(yoyhtml(y, f, title = "A & B", better = "lower"), f)
  h <- readLines(f, warn = FALSE, encoding = "UTF-8")
  expect_true(any(grepl("<title>A &amp; B</title>", h, fixed = TRUE)))
  expect_true(any(grepl("<tr class=\"bad\">", h, fixed = TRUE)))
  expect_true(any(grepl("<tr class=\"good\">", h, fixed = TRUE)))
  expect_true(any(grepl("class=\"bar ", h, fixed = TRUE)))
  yoyhtml(y, f, bars = FALSE)
  expect_false(any(grepl("class=\"bar ", readLines(f, warn = FALSE),
                         fixed = TRUE)))
  expect_error(yoyhtml(y, c(f, f)), "yoyhtml: `file` must be one file path.",
               fixed = TRUE)
})

test_that("the PDF is a file of its own, paginated when it is long", {
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f), add = TRUE)
  expect_identical(yoypdf(y, f), f)
  expect_identical(readBin(f, "raw", 8L), charToRaw("%PDF-1.4"))
  yoypdf(yoypx(seq_len(80) + 100, seq_len(80)), f, height = 4)
  p <- rawToChar(readBin(f, "raw", file.size(f)))
  expect_match(p, "(page 1 of 7)", fixed = TRUE)
  expect_match(p, "/Count 7 >>", fixed = TRUE)
  expect_error(yoypdf(y, f, width = 1),
               "yoypdf: `width` must be one number from 2 to 200.",
               fixed = TRUE)
  expect_error(yoypdf(y, f, height = NA), "`height` must be one number")
  expect_error(yoypdf(y, NA), "yoypdf: `file` must be one file path.",
               fixed = TRUE)
})

test_that("the writer refuses a file or a format it cannot take", {
  f <- tempfile(fileext = ".txt")
  expect_error(yoywr(y, f, format = "doc"), "yoywr: unknown format: doc.",
               fixed = TRUE)
  expect_error(yoywr(y, 1), "yoywr: `file` must be one file path.",
               fixed = TRUE)
  expect_identical(.yasc("a—b ▲"), "a-b +")
})

test_that("yoywr() refuses in its own name", {
  expect_error(yoywr(NULL, tempfile(fileext = ".csv")), "^yoywr: ")
})

test_that("yoypdf() checks lag as yoysmry() does", {
  y <- yoypx(c(10, 12), 1:2)
  expect_error(yoypdf(y, tempfile(fileext = ".pdf"), lag = NULL),
               "^yoypdf: `lag` must be one positive number\\.$")
})
