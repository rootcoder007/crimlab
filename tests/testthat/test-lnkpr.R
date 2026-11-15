# Crime linkage. The kernel is checked by hand in src/t_lnk.cpp; here the
# binding against dist(), glm() and the Mann-Whitney statistic of
# wilcox.test().

lksim <- function() {
  set.seed(9)
  s <- rep(1:8, each = 5)
  cx <- runif(8, 0, 5000)[s] + rnorm(40, 0, 600)
  cy <- runif(8, 0, 5000)[s] + rnorm(40, 0, 600)
  day <- runif(8, 0, 300)[s] + rnorm(40, 0, 30)
  mo <- matrix(rbinom(40 * 5, 1, 0.3), 40)
  s[c(3, 17, 30)] <- NA
  list(x = cbind(x = cx, y = cy), t = day, s = s, mo = mo)
}

test_that("pair measures, the logistic model and the areas follow base R", {
  d <- lksim()
  r <- lnkpr(d$x, t = d$t, series = d$s, mo = d$mo)
  expect_s3_class(r, "lnkpr")
  p <- r$pairs
  expect_identical(nrow(p), 780L)
  dm <- as.matrix(dist(d$x))
  expect_equal(p$dist, dm[cbind(p$i, p$j)], tolerance = 1e-13)
  expect_equal(p$dt, abs(d$t[p$i] - d$t[p$j]), tolerance = 1e-14)
  jac <- function(a, b) if (any(a | b)) sum(a & b) / sum(a | b) else 0
  expect_equal(
    p$jac,
    mapply(function(i, j) jac(d$mo[i, ] == 1, d$mo[j, ] == 1), p$i, p$j),
    tolerance = 1e-15
  )
  expect_identical(
    p$linked,
    ifelse(is.na(d$s[p$i]) | is.na(d$s[p$j]), NA, d$s[p$i] == d$s[p$j])
  )
  tr <- p[!is.na(p$linked), ]
  g <- glm(
    linked ~ dist + dt + jac,
    family = binomial,
    data = tr,
    control = glm.control(epsilon = 1e-14, maxit = 100)
  )
  expect_equal(r$coef$est, unname(coef(g)), tolerance = 1e-7)
  expect_equal(r$coef$se, unname(sqrt(diag(vcov(g)))), tolerance = 1e-6)
  expect_equal(
    p$prob,
    unname(predict(g, p, type = "response")),
    tolerance = 1e-7
  )
  mw <- function(s, l) {
    w <- suppressWarnings(wilcox.test(s[l], s[!l], exact = FALSE))$statistic
    unname(w) / (sum(l) * sum(!l))
  }
  expect_equal(r$auc[["dist"]], mw(-tr$dist, tr$linked), tolerance = 1e-12)
  expect_equal(r$auc[["jac"]], mw(tr$jac, tr$linked), tolerance = 1e-12)
  expect_output(print(r), "780 pairs")
})

test_that("without behaviours the model has distance and time only", {
  d <- lksim()
  r <- lnkpr(d$x, t = d$t, series = d$s)
  expect_identical(r$coef$term, c("(Intercept)", "dist", "dt"))
  expect_true(is.na(r$auc[["jac"]]))
  expect_identical(r$nmo, 0L)
})

test_that("lnkpr() refuses bad input in words", {
  x <- cbind(c(0, 1, 2), c(0, 1, 2))
  errs(
    "lnkpr",
    list(
      list("`series` must have one label", x, t = 1:3, series = 1:2),
      list("linked and unlinked pairs", x, t = 1:3, series = c(1, 1, 1)),
      list(
        "`mo` must hold 0 or 1",
        x,
        t = 1:3,
        series = c(1, 1, 2),
        mo = matrix(2, 3, 1)
      ),
      list("`t`", x, t = 1:2, series = c(1, 1, 2))
    )
  )
})
