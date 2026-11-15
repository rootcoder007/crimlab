# The plot methods of the exposure-and-crime family (R/apxplt.R); those
# of the regression discontinuity are in test-rdplt.R, and keeps() and
# chk() in helper-plot.R. Nothing here looks at pixels:
# each method is drawn on a null device and checked for what it returns,
# that it leaves par() as it found it, that it takes main, col and pal and
# refuses bad ones; the positions the methods compute are checked against
# numbers worked by hand. The data are simulated and small.

# a small panel: 12 areas over 30 days (ten clusters at least, so that
# no fit warns), two instruments, a moderator
axpan <- function() {
  set.seed(20261115)
  d <- expand.grid(day = 1:30, area = 1:12)
  n <- nrow(d)
  d$z1 <- rnorm(n)
  d$z2 <- rnorm(n)
  d$w <- rnorm(n)
  d$pm <- 20 + 3 * d$z1 - 2 * d$z2 + rnorm(n)
  d$o3 <- 30 + 0.5 * d$pm + rnorm(n, 0, 3)
  d$crimes <- rpois(n, exp(0.5 + 0.03 * d$pm + 0.1 * d$w))
  d$theft <- rpois(n, 2)
  d$yl <- 1 + 0.5 * d$pm + 0.2 * d$w + rnorm(n)
  d$smoke <- rbinom(n, 1, 0.4)
  d$grp <- ifelse(d$area <= 6, "poor", "not poor")
  d
}

test_that("the positions the plots compute, by hand", {
  # z = qnorm(0.975) and qnorm(0.95)
  ci <- .axci(c(1, -2), c(0.5, 1), 0.95)
  expect_equal(ci$lo, c(1 - 0.5 * 1.959963984540054, -2 - 1.959963984540054),
               tolerance = 1e-14)
  expect_equal(.axci(0, 2, 0.9)$hi, 2 * 1.644853626951472, tolerance = 1e-14)
  # the rings as cut() names them: [0,500] and (500,1600]
  g <- .axring(c("down:[0,500]", "down:(500,1.6e+03]"))
  expect_identical(g$lo, c(0, 500))
  expect_identical(g$hi, c(500, 1600))
  expect_identical(g$mid, c(250, 1050))
  expect_identical(.axrun(c("run1", "run2", "run3p")), c(1, 2, 3))
})

test_that("ppmlfe, ivfe, ppmlcf and expgmm draw their estimates", {
  d <- axpan()
  f <- ppmlfe(d, "crimes", c("pm", "w"), fe = c("area", "day"))
  p <- chk(f)
  expect_identical(p$label, f$coef$term)
  expect_equal(p$est, exp(f$coef$est), tolerance = 1e-14)
  expect_equal(p$hi, exp(f$coef$hi), tolerance = 1e-14)
  expect_identical(drw(plot(f, terms = "pm"))$label, "pm")
  expect_error(drw(plot(f, terms = "zz")),
               "plot: `terms` must be NULL or names of terms of the fit: pm, w")
  iv <- ivfe(d, "yl", endo = "pm", inst = c("z1", "z2"), x = "w",
             fe = c("area", "day"))
  q <- chk(iv)
  expect_named(q, c("coef", "fstab", "first"))
  expect_identical(q$coef$est, iv$coef$est)
  z <- stats::qnorm(0.975)
  expect_equal(q$fstab$lo, iv$fstab$est - z * iv$fstab$se, tolerance = 1e-14)
  fs <- iv$first
  expect_identical(q$first$est, c(fs$wald, fs$kp, fs$eff, fs$cv))
  ls <- ivfe(d, "yl", x = c("pm", "w"), fe = c("area", "day"))
  expect_named(drw(plot(ls)), "coef")
  cf <- ppmlcf(d, "crimes", endo = "pm", inst = c("z1", "z2"), x = "w",
               fe = c("area", "day"), how = "plug", se = "naive")
  r <- chk(cf)
  expect_named(r, c("coef", "first"))
  expect_equal(r$coef$est, exp(cf$coef$est), tolerance = 1e-14)
  g <- expgmm(d, "crimes", endo = "pm", inst = c("z1", "z2"), fe = "day",
              se = "hc1")
  expect_equal(chk(g)$lo, exp(g$coef$lo), tolerance = 1e-14)
})

test_that("plcday, runfx, hetfx, dosebin and typefx draw their profiles", {
  d <- axpan()
  pc <- plcday(d, "pm", "area", "day", shifts = c(2, -2), y = "crimes",
               x = "pm", fe = c("area", "day"))
  p <- chk(pc)
  expect_identical(p$label, c("2 days before", "same day", "2 days after"))
  expect_identical(p$est, pc$tab$est[c(3L, 1L, 2L)])
  expect_identical(p$group, c("other days", "same day", "other days"))
  rf <- runfx(d, "crimes", "smoke", "area", "day", fe = c("area", "day"))
  r <- chk(rf)
  expect_identical(r$x, c(1, 2, 3))
  expect_identical(r$est, rf$coef$est)
  h <- hetfx(d, "pm", "grp", y = "crimes", x = "pm", fe = c("area", "day"))
  expect_identical(chk(h)$est, h$tab$est)
  db <- dosebin(d, "crimes", "pm", nq = 3, fe = c("area", "day"))
  b <- chk(db)
  expect_identical(b$x, db$tab$mid)
  expect_identical(b$mark, c(TRUE, FALSE, FALSE))
  tf <- typefx(d, c("crimes", "theft"), "pm", x = "pm", fe = c("area", "day"))
  t <- chk(tf)
  expect_identical(t$label, paste0(tf$tab$y, " (adjusted p ",
                                   signif(tf$tab$padj, 2), ")"))
  expect_identical(t$hi, tf$tab$hi)
})

test_that("jntx draws its conditional effects and its surface", {
  d <- axpan()
  j <- jntx(d, "crimes", "pm", "o3", form = "linear", ngrid = 4,
            fe = c("area", "day"))
  p <- chk(j)
  expect_named(p, c("cond", "surf"))
  expect_equal(p$cond$est, exp(j$cond$est), tolerance = 1e-14)
  expect_identical(p$cond$label[2L],
                   paste0("x2 at 50% (", signif(j$cond$x2[2L], 3), ")"))
  expect_identical(p$surf, j$surf)
})

test_that("a case-crossover fit is drawn as the fit of its estimator", {
  set.seed(9)
  d <- expand.grid(day = as.Date("2023-01-01") + 0:89, area = 1:3)
  d$pm <- rnorm(nrow(d), 15, 4)
  d$crimes <- rpois(nrow(d), exp(1 + 0.02 * d$pm))
  h <- ccross(d, "crimes", "area", "day", ma = "pm", x = "pm", se = "hc1")
  expect_s3_class(h, "ppmlfe")
  expect_equal(drw(plot(h, terms = "pm"))$est,
               exp(h$coef$est[h$coef$term == "pm"]), tolerance = 1e-14)
})

test_that("updwn draws the rings, the placebo segments and the rows", {
  set.seed(5)
  a <- data.frame(area = 1:30, ax = runif(30, 0, 5000),
                  ay = c(rep(c(-800, 800), 6), rep(c(-300, 300), 2),
                         runif(14, -4000, 4000)))
  d <- merge(a, data.frame(day = 1:60))
  d$wind <- runif(60, 0, 360)[d$day]
  d$crimes <- rpois(nrow(d), 3)
  road <- data.frame(x0 = 0, y0 = 0, x1 = 5000, y1 = 0)
  u <- updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
             offs = c(-1600, 0, 1600))
  p <- chk(u)
  expect_named(p, c("coef", "placebo", "class"))
  expect_identical(p$coef$est, u$coef$est)
  expect_identical(p$placebo$x, c(-1600, 0, 1600))
  expect_identical(p$placebo$mark, c(FALSE, TRUE, FALSE))
  expect_identical(p$class$value, as.numeric(u$class))
  ur <- updwn(d, "crimes", "ax", "ay", "day", "wind", road, band = 1600,
              dist = c(0, 500, 1600))
  q <- drw(plot(ur))
  expect_named(q, c("coef", "class"))
  expect_identical(q$coef$x, c(250, 1050))
  expect_identical(q$coef$hi, ur$coef$hi)
})
