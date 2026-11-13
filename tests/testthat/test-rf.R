# The forest against ranger on known functions, with ranger's errors from
# dev/ref-rf.R; the kernel's own checks are in src/t_rf.cpp.
ref <- readRDS(test_path("ref-rf.rds"))

test_that("rfor's error is within 15 per cent of ranger's", {
  e <- mean(vapply(1:5, function(s) {
    mean((rfor(ref$x, ref$y, ref$xt, seed = s)$pred - ref$ft)^2)
  }, 0))
  expect_lt(e, 1.15 * ref$reg)
  e <- mean(vapply(1:5, function(s) {
    mean((rfor(ref$x, ref$yb, ref$xt, minnode = 10, seed = s)$pred -
            ref$pt)^2)
  }, 0))
  expect_lt(e, 1.15 * ref$prob)
})

test_that("rfor is fixed by its seed and blind to increasing transforms", {
  a <- rfor(ref$x, ref$y, ref$xt, ntree = 50, seed = 3)
  expect_identical(rfor(ref$x, ref$y, ref$xt, ntree = 50, seed = 3), a)
  expect_false(identical(rfor(ref$x, ref$y, ref$xt, ntree = 50,
                              seed = 4)$pred, a$pred))
  b <- rfor(exp(ref$x), ref$y, exp(ref$xt), ntree = 50, seed = 3)
  expect_identical(b, a)
  # rows in another order: another forest, the same error up to noise
  o <- sample(nrow(ref$x))
  c1 <- rfor(ref$x[o, ], ref$y[o], ref$xt, ntree = 200, seed = 3)
  expect_lt(abs(mean((c1$pred - ref$ft)^2) -
                  mean((rfor(ref$x, ref$y, ref$xt, ntree = 200,
                             seed = 3)$pred - ref$ft)^2)), 0.1)
  expect_length(a$oob, nrow(ref$x))
  s <- rfor(ref$x, ref$y, ntree = 20, mtry = 6, frac = 0.5, repl = FALSE)
  expect_length(s$pred, nrow(ref$x))
})

test_that("rfor stops on bad input, naming the argument", {
  x <- ref$x[1:20, ]
  y <- ref$y[1:20]
  errs("rfor", list(
    list("`x` must be a numeric matrix", y, y),
    list("at least two rows", x[1, , drop = FALSE], y[1]),
    list("`y` must have length 20", x, y[-1]),
    list("`newx` must have the 6 columns", x, y, x[, -1]),
    list("`ntree`", x, y, ntree = 0),
    list("`mtry` must be one whole", x, y, mtry = 0),
    list("`mtry` must be at most", x, y, mtry = 7),
    list("`minnode`", x, y, minnode = 0),
    list("`frac` must be greater than 0", x, y, frac = 0),
    list("`frac` must be at most 1", x, y, frac = 2),
    list("`frac` leaves no row", x, y, frac = 0.01),
    list("`repl`", x, y, repl = NA),
    list("`seed`", x, y, seed = -2)
  ))
})
