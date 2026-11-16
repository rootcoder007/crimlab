test_that("ovrtxt() keeps blocks on lines of their own", {
  h <- paste0("<html><head><title>R</title><style>p {color: red}</style>",
              "</head><body><h2>The Team</h2>\n<p>Number of SIU ",
              "Investigators\nassigned: 3</p><!-- a <p> note -->",
              "<p>WO\n#1 <b>Interviewed</b></p><script>var a = '<p>';",
              "</script><p>WO&nbsp;#2</p></body></html>")
  # written out: the title, the heading and each paragraph on a line,
  # one empty line between blocks (a block ends a line where it starts
  # and where it ends), the page's own breaks spaces
  expect_identical(
    ovrtxt(h),
    paste("R", "", "The Team", "", "Number of SIU Investigators assigned: 3",
          "", "WO #1 Interviewed", "", "WO #2", sep = "\n")
  )
  expect_identical(ovrtxt("a<br/>b<BR>c"), "a\nb\nc")
})

test_that("ovrtxt() decodes entities and leaves plain text its lines", {
  expect_identical(
    ovrtxt("<p>&eacute;t&#233; &amp; &#x2019;s &bogus; &lt;b&gt;</p>"),
    "\u00e9t\u00e9 & 's &bogus; <b>"
  )
  expect_identical(ovrtxt("a  b\r\nc\n\n\n\nd ", html = FALSE),
                   "a b\nc\n\nd")
  # no markup: a '<' that opens no tag
  expect_identical(ovrtxt("1 < 2\nand 3 > 2"), "1 < 2\nand 3 > 2")
  expect_identical(ovrtxt("<p>a</p>", html = FALSE), "<p>a</p>")
  expect_identical(ovrtxt(c(NA, "x")), c(NA, "x"))
  expect_identical(ovrtxt(NA), NA_character_)
  expect_identical(ovrtxt(character(0)), character(0))
})

test_that("ovrtxt() refuses what is not text", {
  expect_error(ovrtxt(1),
               "ovrtxt: `x` must be text (a character vector), not numeric.",
               fixed = TRUE)
  expect_error(ovrtxt(list("a")),
               "ovrtxt: `x` must be text (a character vector), not list.",
               fixed = TRUE)
  expect_error(ovrtxt("a", html = "yes"),
               "ovrtxt: `html` must be NA, TRUE or FALSE.", fixed = TRUE)
  expect_error(ovrtxt("a", html = c(TRUE, FALSE)),
               "ovrtxt: `html` must be NA, TRUE or FALSE.", fixed = TRUE)
})

test_that("ovrdte() reads dates written out in English and French", {
  x <- c("January 5, 2023", "Jan. 5th 2023", "5 January 2023",
         "5th of January, 2023", "5 janvier 2023", "2023-01-05",
         "le 1er février 2021", "3 AOÛT 2017", "SEPT. 9, 2019")
  expect_identical(ovrdte(x),
                   c(rep("2023-01-05", 6), "2021-02-01", "2017-08-03",
                     "2019-09-09"))
  # a day that does not exist is passed over; the next date is read
  expect_identical(ovrdte("February 30, 2023, then March 2, 2023"),
                   "2023-03-02")
  expect_identical(ovrdte("February 29, 2024"), "2024-02-29")
  expect_identical(ovrdte(c("May 2023", "no date", NA, "2023-13-01")),
                   rep(NA_character_, 4))
  expect_identical(ovrdte(character(0)), character(0))
})

test_that("ovrdte() refuses what is not text", {
  expect_error(ovrdte(20230105),
               "ovrdte: `x` must be text (a character vector), not numeric.",
               fixed = TRUE)
})

test_that("ovrsch() gives the fields of the preset and their rules", {
  s <- ovrsch()
  expect_identical(names(s), c("name", "rule", "cue", "alt", "sect", "desc"))
  expect_identical(nrow(s), 17L)
  expect_true(all(s$rule %in% c("after", "num", "date", "tag", "sent",
                                "prev", "lines", "name", "flag", "before",
                                "word")))
  expect_identical(s$cue[s$name == "number_of_witness_officials"], "WO")
  expect_identical(anyDuplicated(s$name), 0L)
  expect_identical(ovrsch("siu"), s)
})

test_that("ovrsch() refuses a preset it does not have", {
  expect_error(ovrsch("iio"),
               "ovrsch: `preset` must be one of \"siu\", not \"iio\".",
               fixed = TRUE)
})
