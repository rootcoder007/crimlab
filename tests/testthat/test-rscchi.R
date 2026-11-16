# The reference values are those of survey::svychisq() (survey 4.5) on
# these data: its "Chisq" p-value, which divides Pearson's statistic by
# the mean design effect, and its "F" statistic and degrees of freedom,
# with chisq.test() on svytable(Ntotal = 60) for Pearson's statistic.
rs_dat <- function() {
  i <- 1:60
  data.frame(
    x = c("a", "b", "c")[(i * 7 + i %/% 5) %% 3 + 1],
    y = c("no", "yes")[((i * 5 + i %/% 7) %% 2) + 1],
    z = c("lo", "mid", "hi")[(i * 11 + i %/% 4) %% 3 + 1],
    w = 1 + (i * 13) %% 9,
    st = rep(c("s1", "s2", "s3"), each = 20),
    psu = rep(1:12, each = 5),
    N = rep(c(10, 8, 12), each = 20),
    stringsAsFactors = FALSE
  )
}

rs_check <- function(r, x2, p2, rs1p, f, ndf, ddf, fp, k) {
  t <- r$test
  expect_identical(t$test, c("pearson", "rs1", "rs2", "f"))
  expect_equal(t$stat[1L], x2, tolerance = 1e-12)
  expect_equal(t$p[1L], p2, tolerance = 1e-10)
  expect_identical(t$df1[1:2], c(k, k))
  # first order: X^2 over the mean design effect, which is X^2 / (F K)
  expect_equal(t$stat[2L], f * k, tolerance = 1e-9)
  expect_equal(t$p[2L], rs1p, tolerance = 1e-9)
  # second order: X^2 d0 / tr(Delta) = F d0 on d0
  expect_equal(t$stat[3L], f * ndf, tolerance = 1e-9)
  expect_equal(t$df1[3L], ndf, tolerance = 1e-9)
  expect_equal(t$p[3L], stats::pchisq(f * ndf, ndf, lower.tail = FALSE),
               tolerance = 1e-9)
  expect_equal(c(t$stat[4L], t$df1[4L], t$df2[4L], t$p[4L]),
               c(f, ndf, ddf, fp), tolerance = 1e-9)
  expect_true(all(is.na(t$df2[1:3])))
  expect_equal(sum(r$deff), x2 / f, tolerance = 1e-9)
  expect_equal(r$dbar, x2 / f / k, tolerance = 1e-9)
}

test_that("rscchi() agrees with svychisq() with strata, clusters, weights", {
  d <- rs_dat()
  r <- rscchi(d$x, d$y, d$w, d$st, d$psu)
  rs_check(r, 1.1169054738919768, 0.57209355862601141, 0.76602231334240189,
           0.26654397997598717, 1.6812477000741737, 15.131229300667563,
           0.73212806783300088, 2)
  expect_identical(r$des, c(n = 60, strata = 3, psu = 12, df = 9))
  tw <- tapply(d$w, list(d$x, d$y), sum)
  expect_equal(r$tab, tw / sum(tw), tolerance = 1e-14)
})

test_that("rscchi() agrees with svychisq() with weights alone", {
  d <- rs_dat()
  r <- rscchi(d$x, d$z, d$w)
  rs_check(r, 2.1123358348968111, 0.71510642495067211, 0.78528903370924796,
           0.43246171404254102, 3.9721902857918954, 234.35922686172182,
           0.78380486574471098, 4)
  expect_length(r$deff, 4L)
  expect_identical(r$des[["df"]], 59)
})

test_that("rscchi() agrees with svychisq() with a finite population", {
  d <- rs_dat()
  r <- rscchi(d$x, d$z, d$w, d$st, d$psu, d$N)
  rs_check(r, 2.1123358348968111, 0.71510642495067211, 0.67778480426476706,
           0.57909998156172837, 2.8492036631321547, 25.642832968189392,
           0.62563216751166972, 4)
  u <- rscchi(d$z, d$y, strata = d$st, psu = d$psu)
  rs_check(u, 0.13574660633484165, 0.93437884891875178, 0.94276904081808865,
           0.058933945945440304, 1.351260661531311, 12.161345953781799,
           0.87909507710300094, 2)
  out <- utils::capture.output(v <- print(u))
  expect_identical(out[1L], paste("rscchi: 3 x 2 table, 60 respondents,",
                                  "3 strata, 12 first-stage units"))
  expect_identical(v, u)
})

test_that("rscchi() refuses tables and designs it cannot test", {
  d <- rs_dat()
  expect_error(rscchi(d$x, d$y[1:59]),
               "rscchi: `y` must have one label per value, 60.", fixed = TRUE)
  expect_error(rscchi(d$x, rep("a", 60)),
               "rscchi: `x` and `y` must each hold two labels or more.",
               fixed = TRUE)
  expect_error(rscchi(d$x, d$y, strata = d$st, psu = d$st),
               paste("rscchi: a stratum has one first-stage unit: its",
                     "variance cannot be estimated."), fixed = TRUE)
  expect_error(rscchi(d$x, d$y, strata = d$st, psu = d$psu, fpc = d$psu),
               "rscchi: `fpc` must be constant within each stratum.",
               fixed = TRUE)
  expect_error(rscchi(d$x, d$y, w = -d$w),
               "rscchi: `w` must be greater than 0.", fixed = TRUE)
})
