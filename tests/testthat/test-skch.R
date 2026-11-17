# One-pass summaries: the values of the reference implementation
# (rmoriebricklayer 1.1.1, computed once on these inputs and printed with
# 17 significant digits), exactness against base R, refusals and printing.

ref <- local({
  mo1 <- c(
    6, 33333335.5, 13333332516666722, 4.4444440361110829e+23,
    4.4444439000000542e+31
  )
  mo2 <- c(
    9, 22222224.902777776, 15555554533333416, 8.6419744567901082e+23,
    7.4897109497943084e+31
  )
  mo3 <- c(
    9, 22222224.902777776, 1944444316666677, 44095853.735546127,
    1.3363062095621068, -0.21428571428572996
  )
  mo4 <- c(
    2, 1.5, 0.5, 0.70710678118654757, NaN, NaN
  )
  rv1 <- c(
    35, 69, 306, 632, 640, 786, 798
  )
  rv2 <- c(
    12, 16, 22, 32, 37
  )
  rv3 <- numeric(0)
  hl1 <- c(
    0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0
  )
  hl2 <- c(
    9, 7, 5, 7, 7, 8, 6, 11, 7, 5, 7, 7, 5, 8, 7, 6, 7, 7, 7, 6, 6, 12, 8, 8,
    7, 5, 7, 7, 7, 7, 5, 6, 6, 6, 6, 12, 5, 11, 10, 8, 7, 6, 12, 6, 6, 6, 8, 6,
    7, 9, 6, 7, 5, 8, 7, 4, 10, 10, 9, 6, 7, 7, 8, 12, 4024.0345548037885
  )
  hl3 <- 2.0652326764022773
  environment()
})

test_that("the one-pass summaries equal the reference implementation", {
  x <- c(3.5, -1, 2.25, 8, 1e8 + 0.5, 1e8 - 0.25, NA, 4, 0.125, 7)
  a <- smryupd(onlsmry(x[1:4]), x[5:7])
  b <- onlsmry(x[8:10])
  expect_equal(a$state, ref$mo1, tolerance = 1e-12)
  expect_equal(smrymrg(a, b)$state, ref$mo2, tolerance = 1e-12)
  expect_equal(unname(smrystat(smrymrg(a, b))), ref$mo3, tolerance = 1e-10)
  expect_equal(unname(smrystat(onlsmry(c(1, 2)))), ref$mo4,
               tolerance = testthat_tolerance())
  expect_identical(rsvidx(1000, 7, seed = 1), ref$rv1)
  expect_identical(rsvidx(50, 5, seed = 123456789), ref$rv2)
  expect_identical(rsvidx(10, 0), ref$rv3)
  expect_identical(
    as.numeric(dstnskch(c("a", "b", "c", "a", "1", 1, 2.5, TRUE, NA), p = 4)),
    ref$hl1
  )
  s <- dstnskch(as.character(1:3000), p = 6)
  # the registers are whole numbers; the estimate ends in a logarithm,
  # whose last digit differs between systems
  expect_equal(c(as.numeric(s), dstncnt(s)), ref$hl2, tolerance = 1e-14)
  expect_identical(dstncnt(dstnskch(c("café", "x"), p = 5)), ref$hl3)
  # the same text in latin1 is the same value
  lat <- iconv("café", "UTF-8", "latin1")
  expect_identical(dstnskch(lat, p = 5), dstnskch("café", p = 5))
})

test_that("block accumulation is exact against base R", {
  set.seed(41)
  x <- stats::rnorm(1000, mean = 7, sd = 3)
  acc <- onlsmry()
  for (i in seq(1, 1000, by = 100)) acc <- smryupd(acc, x[i:(i + 99)])
  s <- smrystat(acc)
  expect_named(s, c("n", "mean", "variance", "sd", "skewness", "kurtosis"))
  expect_identical(s[["n"]], 1000)
  expect_equal(s[["mean"]], mean(x), tolerance = testthat_tolerance())
  expect_equal(s[["variance"]], stats::var(x), tolerance = 1e-12)
  for (bs in c(1L, 7L, 333L)) {
    a2 <- onlsmry()
    for (i in seq(1, 1000, by = bs)) {
      a2 <- smryupd(a2, x[i:min(1000, i + bs - 1L)])
    }
    expect_equal(smrystat(a2)[["variance"]], stats::var(x),
                 tolerance = testthat_tolerance())
  }
  expect_equal(smrystat(onlsmry(x)), s, tolerance = testthat_tolerance())
  # far from zero the digits survive
  y <- 1e8 + stats::rnorm(1000)
  b <- onlsmry()
  for (i in seq(1, 1000, by = 100)) b <- smryupd(b, y[i:(i + 99)])
  expect_equal(smrystat(b)[["variance"]], stats::var(y), tolerance = 1e-9)
  # merging in either order, with empty sides
  a <- onlsmry(x[1:250])
  b <- onlsmry(x[251:1000])
  expect_equal(smrystat(smrymrg(a, b)), smrystat(smrymrg(b, a)),
               tolerance = testthat_tolerance())
  expect_identical(smrymrg(onlsmry(), a)$state, a$state)
  expect_identical(smrymrg(a, onlsmry())$state, a$state)
  expect_identical(smryupd(a, c(NA, NA))$state, a$state)
  e <- smrystat(onlsmry())
  expect_true(all(is.nan(e[-1L])))
  expect_true(is.nan(smrystat(onlsmry(c(2, 2, 2, 2)))[["skewness"]]))
  expect_equal(smrystat(onlsmry(c("1", "3")))[["mean"]], 2,
               tolerance = testthat_tolerance())
})

test_that("the summaries refuse what they cannot hold", {
  expect_error(onlsmry(c("a", "b")), "onlsmry: `x` must be numeric, not text")
  expect_error(onlsmry(TRUE), "must be a number, not TRUE")
  expect_error(onlsmry(list(1)), "got list")
  expect_error(onlsmry(c(1, Inf)), "infinite")
  expect_error(smryupd(list(), 1), "smryupd: `acc` must come from onlsmry")
  expect_error(smrymrg(onlsmry(), 1), "both arguments")
  expect_error(smrystat(1), "smrystat: `acc`")
})

test_that("the reservoir is uniform and reproducible", {
  expect_identical(rsvidx(1000, 5, seed = 1), rsvidx(1000, 5, seed = 1))
  expect_identical(rsvidx(3, 10), c(1, 2, 3))
  hits <- unlist(lapply(1:400, function(s) rsvidx(50, 5, s)))
  expect_lt(abs(mean(hits) - 25.5), 2)
  expect_identical(rsvsmp(letters, 4, seed = 2), letters[rsvidx(26, 4, 2)])
  expect_error(rsvidx(10, -1), "rsvidx: `k` must be a single non-negative")
  expect_error(rsvidx(10.5, 1), "`n` must be a single")
  expect_error(rsvidx(3e9, 1), "at most 2147483647")
  expect_error(rsvidx(10, 1, seed = -1), "`seed` must be a whole number")
  expect_error(.Call(C_clab_r_rsv, 10, 2e8, 1), "rsvidx:")
})

test_that("the sketch counts near the truth and refuses bad registers", {
  set.seed(1)
  x <- sample(1:5000, 50000, replace = TRUE)
  est <- dstncnt(dstnskch(x))
  expect_lt(abs(est / length(unique(x)) - 1), 0.05)
  s <- dstnskch(x[1:25000])
  s <- dstnskch(x[25001:50000], regstr = s)
  expect_identical(s, skchmrg(dstnskch(x[1:25000]), dstnskch(x[25001:50000])))
  expect_identical(dstncnt(NULL), 0)
  expect_identical(dstncnt(dstnskch(character(0))), 0)
  expect_error(dstnskch(1, p = 3), "between 4 and 20")
  expect_error(dstnskch(1, p = 4, regstr = 1:3), "2\\^p entries")
  expect_error(dstncnt(c(1, 99)), "non-negative integer registers")
  expect_error(dstncnt("a"), "got character")
  expect_error(skchmrg(dstnskch(1:10, p = 4), dstnskch(1:10, p = 5)), "same")
})

test_that("an accumulator prints", {
  out <- capture.output(print(onlsmry(c(1, 2, 3, 10))))
  expect_match(out, "Accumulated summary: n 4")
})

test_that("rsvsmp() refuses in its own name", {
  expect_error(rsvsmp(1:5, k = NULL), "^rsvsmp: `k` must be numeric")
})

test_that("dstnskch() and rsvsmp() refuse what is not a vector", {
  expect_error(dstnskch(sum), "^dstnskch: `x` must be a vector")
  expect_error(rsvsmp(new.env(), 2), "^rsvsmp: `x` must be a vector")
  bad <- rawToChar(as.raw(c(0xff, 0xfe, 0x41)))
  expect_error(rsvidx(10, bad), "^rsvidx: `k` holds text that is not valid")
})

test_that("an accumulator without its state is refused in words", {
  a <- structure(list(1), class = c("onlsmry", "list"))
  expect_error(smrystat(a), "^smrystat: `acc` must come from onlsmry")
  expect_error(smryupd(a, 1), "^smryupd: `acc` must come from onlsmry")
  expect_error(smrymrg(a, onlsmry()), "^smrymrg: both arguments")
})

test_that("SHA-256 gives the FIPS 180-4 test vectors", {
  h <- function(x) paste(.Call(C_clab_r_sha256, charToRaw(x)), collapse = "")
  expect_identical(
    h(""),
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  )
  expect_identical(
    h("abc"),
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
  )
  # 56 bytes: the length no longer fits the first padding block
  expect_identical(
    h("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
    "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
  )
  expect_identical(
    h(paste0(
      "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn",
      "hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu"
    )),
    "cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1"
  )
  expect_identical(
    h(strrep("a", 1e6)),
    "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
  )
})
