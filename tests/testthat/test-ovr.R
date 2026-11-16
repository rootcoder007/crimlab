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
