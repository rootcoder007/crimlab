# bndxw() and synexp(): against sums and a logistic fit computed here in
# base R, and the refusals.

test_that("bndxw() carries values through the correspondence", {
  set.seed(13)
  from <- c("a", "b", "c", "c", "d", "e", "e")
  to <- c("A", "A", "A", "B", "B", "C", "B")
  w <- c(1, 1, 0.25, 0.75, 1, 0.6, 0.6)
  unit <- sample(c(letters[1:5], "z"), 60, TRUE)
  x <- rpois(60, 4)
  out <- bndxw(from, to, x, unit, w)
  tot <- tapply(x, factor(unit, letters[1:5]), sum)
  tot[is.na(tot)] <- 0
  ref <- tapply(tot[from] * w, to, sum)
  expect_identical(out$tab$to, c("A", "B", "C"))
  expect_equal(out$tab$value, unname(c(ref)), tolerance = 1e-15)
  expect_identical(out$map$from, letters[1:5])
  expect_identical(out$map$to, c("A", "A", "B", "B", "C"))
  expect_identical(out$map$share, c(1, 1, 0.75, 1, 0.6))
  expect_equal(out$map$wsum, c(1, 1, 1, 1, 1.2), tolerance = 1e-15)
  expect_identical(out$lost, as.numeric(sum(x[unit == "z"])))
  expect_identical(out$unmatched, if (any(unit == "z")) "z" else character())
  one <- bndxw(c(1, 2), c(9, 9), c(2, 3), c(2, 2))
  expect_identical(one$tab$value, 5)
  expect_identical(one$unmatched, character())
})

test_that("synexp() applies the survey fit to the areas", {
  set.seed(14)
  s <- data.frame(age = runif(500, 18, 80), male = rbinom(500, 1, 0.5))
  y <- rbinom(500, 1, plogis(-3 + 0.03 * s$age + 0.4 * s$male))
  g <- glm(y ~ age + male, family = binomial(), data = s,
           control = glm.control(epsilon = 1e-14, maxit = 100))
  ax <- data.frame(age = c(35, 52, 61), male = c(0.48, 0.51, 0.4))
  pop <- c(12000, 8000, 0)
  out <- synexp(y, s, ax, pop)
  b <- unname(coef(g))
  expect_identical(out$coef$term, c("(Intercept)", "age", "male"))
  expect_equal(out$coef$est, b, tolerance = 1e-8)
  r <- plogis(drop(cbind(1, as.matrix(ax)) %*% b))
  expect_equal(out$areas$rate, r, tolerance = 1e-8)
  expect_equal(out$areas$expo, r * pop, tolerance = 1e-8)
  expect_equal(out$total, sum(r * pop), tolerance = 1e-8)
  expect_equal(out$dev, deviance(g), tolerance = 1e-8)
  u <- synexp(y, unname(as.matrix(s)), as.matrix(ax), pop)
  expect_identical(u$coef$term, c("(Intercept)", "x1", "x2"))
})

test_that("bndxw() and synexp() refuse bad input", {
  errs("bndxw", list(
    list("`from` must hold at least one link", character(), character(), 1,
         "a"),
    list("`to` must have one label per value", c("a", "b"), "A", 1, "a"),
    list("`w` must be at least 0", "a", "A", 1, "a", w = -1),
    list("`unit` must have one label per value, 2", "a", "A", 1:2, "a"),
    list("`unit` must be labels", "a", "A", 1, NA)
  ))
  x <- cbind(1:6, c(2, 1, 4, 3, 6, 6))
  errs("synexp", list(
    list("`y` must be 0 or 1", c(0, 1, 2, 0, 1, 0), x, x, rep(1, 6)),
    list("`y` must hold both 0 and 1", rep(0, 6), x, x, rep(1, 6)),
    list("`x` needs more rows than columns plus one", c(0, 1, 1), x[1:3, ],
         x, rep(1, 6)),
    list("`ax` must be a numeric matrix of 2 columns", c(0, 1, 0, 1, 0, 1),
         x, 1:2, 1),
    list("`pop` must be at least 0", c(0, 1, 0, 1, 0, 1), x,
         x[1, , drop = FALSE], -1),
    list("no fit", c(0, 0, 0, 1, 1, 1), x, x, rep(1, 6))
  ))
})
