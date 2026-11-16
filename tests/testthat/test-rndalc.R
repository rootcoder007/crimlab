test_that("rndalc() fills every complete block in the ratio", {
  a <- rndalc(60, arms = c("x", "y", "z"), ratio = c(2, 1, 1),
              block = c(4, 8), seed = 5)
  d <- a$alloc
  expect_identical(names(d), c("id", "block", "arm"))
  expect_identical(d$id[c(1L, 60L)], c("P01", "P60"))
  b <- split(d$arm, d$block)
  expect_identical(names(b), as.character(seq_along(b)))
  # every block but the last holds 4 or 8 places, 2:1:1
  for (i in seq_len(length(b) - 1L)) {
    m <- length(b[[i]])
    expect_true(m %in% c(4L, 8L))
    k <- table(factor(b[[i]], c("x", "y", "z")))
    expect_identical(as.vector(k), as.integer(c(2, 1, 1) * m / 4))
  }
  expect_identical(sum(lengths(b)), 60L)
  want <- matrix(as.vector(table(factor(d$arm, c("x", "y", "z")))), 1L,
                 dimnames = list("all", c("x", "y", "z")))
  expect_identical(a$tab, want)
  expect_identical(a$seed, 5)
})

test_that("rndalc() with fixed blocks of four gives two of each arm", {
  a <- rndalc(40, block = 4, seed = 2)
  d <- a$alloc
  expect_identical(d$block, rep(1:10, each = 4))
  k <- tapply(d$arm == "A", d$block, sum)
  expect_identical(as.vector(k), rep(2L, 10))
  expect_identical(a$tab, matrix(c(20L, 20L), 1L,
                                 dimnames = list("all", c("A", "B"))))
})

test_that("rndalc() runs the blocks of each stratum on their own", {
  s <- rndalc(arms = c("drug", "placebo"), ratio = c(2, 1),
              strata = c(Toronto = 30, Ottawa = 24), seed = 7)
  d <- s$alloc
  expect_identical(names(d), c("id", "stratum", "block", "arm"))
  expect_identical(d$stratum, rep(c("Toronto", "Ottawa"), c(30, 24)))
  expect_identical(as.vector(tapply(d$block, d$stratum, min)), c(1L, 1L))
  expect_identical(rownames(s$tab), c("Toronto", "Ottawa"))
  expect_identical(rowSums(s$tab), c(Toronto = 30, Ottawa = 24))
  # the default sizes are 6 and 9; complete blocks hold the arms 2:1
  for (st in c("Toronto", "Ottawa")) {
    e <- d[d$stratum == st, ]
    b <- split(e$arm, e$block)
    for (i in seq_len(length(b) - 1L)) {
      m <- length(b[[i]])
      expect_true(m %in% c(6L, 9L))
      expect_identical(sum(b[[i]] == "drug"), as.integer(2 * m / 3))
    }
  }
  expect_identical(rndalc(54, arms = c("drug", "placebo"), ratio = c(2, 1),
                          strata = c(Toronto = 30, Ottawa = 24), seed = 7),
                   s)
})

test_that("rndalc() randomises simply with block 0 and keeps its seed", {
  r <- rndalc(25, block = 0, seed = 2)
  expect_true(all(is.na(r$alloc$block)))
  expect_true(all(r$alloc$arm %in% c("A", "B")))
  expect_identical(sum(r$tab), 25L)
  expect_identical(rndalc(25, block = 0, seed = 2), r)
  set.seed(11)
  a <- rndalc(30, prefix = "S")
  expect_identical(a$alloc$id[1L], "S01")
  expect_identical(rndalc(30, prefix = "S", seed = a$seed), a)
  expect_false(identical(rndalc(30, seed = 1)$alloc$arm,
                         rndalc(30, seed = 2)$alloc$arm))
  out <- utils::capture.output(v <- print(r))
  expect_identical(out[1L], "rndalc: 25 participants, seed 2")
  expect_identical(v, r)
})

test_that("rndalc() refuses arms, ratios, blocks and strata it cannot use", {
  expect_error(rndalc(10, arms = "A", seed = 1),
               "rndalc: `arms` must name two or more arms, each once.",
               fixed = TRUE)
  expect_error(rndalc(10, arms = c("A", "A"), seed = 1),
               "rndalc: `arms` must name two or more arms, each once.",
               fixed = TRUE)
  expect_error(rndalc(10, ratio = c(1, 1.5), seed = 1),
               "rndalc: `ratio` must be a whole number.", fixed = TRUE)
  expect_error(rndalc(10, ratio = 1, seed = 1),
               "rndalc: `ratio` must give one number per arm.", fixed = TRUE)
  expect_error(rndalc(10, block = 5, seed = 1),
               paste("rndalc: `block` must be 0, or sizes that are multiples",
                     "of sum(ratio) (2)."), fixed = TRUE)
  expect_error(rndalc(10, block = c(0, 4), seed = 1),
               paste("rndalc: `block` must be 0, or sizes that are multiples",
                     "of sum(ratio) (2)."), fixed = TRUE)
  expect_error(rndalc(seed = 1), "rndalc: `n` must be numeric, got NULL.",
               fixed = TRUE)
  expect_error(rndalc(0, seed = 1),
               "rndalc: `n` must be between 1 and 1e+07.", fixed = TRUE)
  expect_error(rndalc(10, strata = c(a = 4, b = 5), seed = 1),
               "rndalc: `n` must equal the places of `strata` (9).",
               fixed = TRUE)
  expect_error(rndalc(strata = c(4, 5), seed = 1),
               "rndalc: `strata` must be named, each stratum once.",
               fixed = TRUE)
  expect_error(rndalc(strata = c(a = 0, b = 0), seed = 1),
               "rndalc: `strata` must hold from 1 to 10000000 places in all.",
               fixed = TRUE)
  expect_error(rndalc(10, seed = -1),
               "rndalc: `seed` must be between 0 and 2147483647.",
               fixed = TRUE)
  expect_error(rndalc(10, seed = 1, prefix = c("a", "b")),
               "rndalc: `prefix` must be one string.", fixed = TRUE)
})
