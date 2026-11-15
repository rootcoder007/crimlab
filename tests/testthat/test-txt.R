# Text to plain ASCII (R/txt.R). The values were computed once by the
# table of rmoriebricklayer 1.1.1 (its transliteration without stringi,
# then its strip to printable ASCII) on the same text.

x <- c("Prof. \u00c1ngela Zorro Medina",
       "Stra\u00dfe \u00c6sir \u0152uvre \u00de\u00fe",
       "\u0141\u00f3d\u017a \u0160koda \u0158\u00edm",
       "\u041c\u043e\u0441\u043a\u0432\u0430 caf\u00e9", NA, "a\tb\nc",
       rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9, 0x21))),
       "\u00d7\u00f7\u00ff")

test_that("toascii() folds as rmoriebricklayer's table does", {
  expect_identical(
    toascii(x),
    c("Prof. Angela Zorro Medina", "Strasse AEsir OEuvre Thth",
      "Lodz Skoda Rim", " cafe", "", "abc", "caf!", "y")
  )
  expect_true(all(unlist(lapply(toascii(x), function(s) {
    as.integer(charToRaw(s)) %in% 32:126
  }))))
  expect_identical(toascii(c(1.5, TRUE)), c("1.5", "TRUE"))
  expect_identical(toascii(character()), character())
  expect_error(toascii(list("a")), "toascii: `x` must be text")
})

test_that("asciifb() keeps valid UTF-8 and repairs the rest", {
  y <- asciifb(x)
  expect_identical(y[-7], x[-7])
  expect_identical(y[7], "caf!")
  expect_identical(asciifb(x, force = TRUE), toascii(x))
  expect_identical(asciifb(NA), NA_character_)
  expect_identical(asciifb(NA, force = TRUE), "")
  expect_identical(asciifb("plain"), "plain")
})

test_that("wrtxtfb() writes the UTF-8 bytes with one newline per line", {
  p <- tempfile()
  on.exit(unlink(p))
  expect_identical(wrtxtfb(c("\u00c1ngela", NA, "z"), p), p)
  expect_identical(
    readBin(p, "raw", 100),
    as.raw(c(0xc3, 0x81, 0x6e, 0x67, 0x65, 0x6c, 0x61, 0x0a, 0x4e, 0x41,
             0x0a, 0x7a, 0x0a))
  )
  wrtxtfb(x[7], p)
  expect_identical(readLines(p), "caf!")
  wrtxtfb(character(), p)
  expect_identical(readBin(p, "raw", 10), as.raw(0x0a))
  expect_error(wrtxtfb("x", c("a", "b")), "one non-empty file path")
  expect_error(wrtxtfb("x", ""), "one non-empty file path")
  expect_error(wrtxtfb(mean, p), "`text` must be text")
})
