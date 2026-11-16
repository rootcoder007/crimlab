# Air pollution: the values of the reference implementation
# (rmoriebricklayer 1.1.1, computed once on these inputs and printed with
# 17 significant digits), closed forms, refusals and printing.

ref <- local({
  sig <- c(
    10.972602697183962, 86.291099460800964, 578.85829274263926, 10, 80, 600,
    15.842360687626789, 118.84501643331319, 647.23186796471236,
    12.296340919151518, 113.58873183551262, 1440, 7.9800746888610634,
    62.757163244218887, 420.98784926737403, 6, 48, 360, 15.842360687626789,
    118.84501643331319, 647.23186796471236, 12.296340919151518,
    113.58873183551262, 1440, 5.4863013485919812, 43.145549730400482,
    289.42914637131963, 3.9801487608399566, 30.79201435678004,
    189.73665961010275, 10.891622972743416, 81.705948797902821,
    444.97190922573975, 10, 80, 600, 3.9900373444305317, 31.378581622109444,
    210.49392463368702, 2.8934569330224731, 18.973665961010276,
    76.752257888019756, 7.9211803438133943, 59.422508216656595,
    323.61593398235618, 6.948083337796513, 52.915026221291818,
    304.70002504620498, 2.9925280083228989, 23.533936216582084,
    157.87044347526526, 1.4778325123152711, 10.714285714285714,
    47.368421052631575, 5.4458114863717082, 40.85297439895141,
    222.48595461286988, 3.8579425773632972, 25.298221281347033,
    102.33634385069301, 1.9950186722152659, 15.689290811054722,
    105.24696231684351, 0.78817733990147798, 5.7142857142857144,
    25.263157894736842, 5.4458114863717082, 40.85297439895141,
    222.48595461286988, 3.8579425773632972, 25.298221281347033,
    102.33634385069301
  )
  brg1 <- c(
    9.4559400000000018, 28.882805133335584, 199.5315252754838,
    3.9261018394952139, 28.882805133335584, 28.882805133335584
  )
  brg2 <- c(
    189.74919599999998, 300.33772257152515, 970.10543819147802,
    14.225210753454133, 300.33772257152515
  )
  brg3 <- c(
    9.4559400000000018, 49.427678724356873, 446.70197379615104,
    3.9261018394952139, 49.427678724356873
  )
  gp1 <- c(
    7.8077308399621045e-06, 0.00010790278628667551, 2.7955357636728076e-05, 0
  )
  gp2 <- c(
    0.00011551886670548154, 1.1546882591268064e-05, 1.8749515633743157e-06, 0
  )
  gs1 <- c(
    5.5484159572262007e-09, 5.3892454092814018e-27, 0, 9.5658543729456589e-16
  )
  ad1 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0.013225000000000001, 0.045689687499999999,
    0.12824406250000001, 0.27661328125000001, 0.29719828124999992, 0, 0,
    0.010748375000000001, 0.044073234374999999, 0.12847568749999999,
    0.50168143750000005, 0.35123428125, 0, 0, 0.0053898125000000005,
    0.018742218749999994, 0.050208500000000003, 0.097479421875,
    0.11238187499999998, 0, 0, 0.0013365625000000001, 0.0052085000000000013,
    0.011822437500000001, 0.022892781250000001, 0.024854765624999994, 0, 0, 0,
    0, 0, 0, 0, 0, 2.1475002031249999, 0.55000000000000004,
    0.30000000000000004, 0.85000000000000009
  )
  ad2 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, -8.1315162936412833e-20, 0.00024291280146000017,
    -6.5052130349130266e-19, 0.0011924174835750011, -1.7347234759768071e-18, 0,
    0, 0.00017664377848000012, -7.0473141211557788e-19, 0.0017387162100630015,
    -3.4694469519536142e-18, 0.0042774650026250038, 0, 0,
    -3.2526065174565133e-19, 0.001174602260007001, -2.8189256484623115e-18,
    0.0057657771753750059, -6.9388939039072284e-18, 0, 0,
    0.0003262585512150003, -1.3010426069826053e-18, 0.0032113148782950034,
    -6.0715321659188248e-18, 0.0079000705890000071, 0, 0, 0, 0, 0, 0, 0, 0,
    0.20000000000000001, 5, 1
  )
  lg1 <- c(
    3.9991399718202651, 4.6592924264582889, 3.4960997360493797,
    4.4826158587558025, 2.5671111510631297, 0.47918103646461263,
    1.964331453532248, 3.7189058994646689, 5.1382503007548097,
    3.6989903821456975, 3.8408518288293734, 2.9999318144724074, 0, 0, 0, 0, 0,
    0, 0.050000000000000003, 0.050000000000000003, 0, 0, 0.10000000000000001,
    0, 0, 0, 0, 0.050000000000000003
  )
  lg2 <- c(
    0.43417736858951228, 2.4759339879889599, 1.6527152957941385,
    2.3339020078983106, 2.571395457582466, 3.1956949320746464
  )
  pf1 <- c(
    3.5073539281017183e-53, 8.3842953086539609e-05, 0.00040608322423182669,
    6.5092140855195981e-19, 0.00083655533181385396, 0.00051333729507465958,
    3.5073539281017183e-53, 8.3842953086539609e-05, 0.00040608322423182669,
    950, 100
  )
  cp1 <- c(
    1.0547841431002283, 0.052333507972567306, 1.0410161921565433,
    1.061686049461863, 12.333333333333334
  )
  cp2 <- c(
    1.4072969162536555, 1.9402562405998198
  )
  cp3 <- c(
    1.6322666287534864, 0.489969618659293
  )
  cn1 <- c(
    1.2523673607612227, 0.21345696824126523, 1.0925889237548163,
    1.3746429483847153
  )
  cn2 <- c(
    0.95504196219071469, -0.045999999999999999
  )
  af1 <- c(
    0.053030303030303066, 0.040307101727447253, 0.059266227657572952,
    -0.11111111111111108
  )
  md1 <- c(
    792.74057836309862, -0.015227266976758469
  )
  pb1 <- c(
    0.33256401757488729, 831.41004393721823, 2500, 1.553634483026495
  )
  pb2 <- c(
    0.23535675698634645, 1412.1405419180787, 1.3077994334439491
  )
  pa_order <- c(
    "B", "A", "C"
  )
  pa1 <- c(
    1.0362876182839384, 1.0159682532195475, 1, 0.0350169370391877,
    0.015717275779971489, 0, 9.8047423709725567, 7.5442923743863153, 0, 280,
    480, 50
  )
  ec1 <- c(
    -0.21225071225071224, 19.5, -2.0694444444444442, 0.25641025641025639,
    0.47008547008547008, 0.68376068376068377, 0.81196581196581197,
    0.9145299145299145, 1
  )
  environment()
})

test_that("the dispersion models equal the reference implementation", {
  tol <- 1e-10
  x <- c(50, 400, 3000)
  r <- NULL
  for (s in LETTERS[1:6]) {
    for (st in c("rural", "urban")) r <- c(r, unlist(pgsig(x, s, st)))
  }
  expect_equal(unname(r), ref$sig, tolerance = tol)
  expect_equal(
    unname(unlist(briggs(c(10, 500, 1e5), 4, 1, 12, 420, 285))), ref$brg1,
    tolerance = tol
  )
  expect_equal(
    unname(unlist(briggs(c(10, 1e5), 3, 3, 20, 500, 285, stab = "B"))),
    ref$brg2,
    tolerance = tol
  )
  expect_equal(
    unname(unlist(
      briggs(c(10, 1e5), 4, 1, 12, 420, 285, stab = "F", dthdz = 0.01)
    )),
    ref$brg3,
    tolerance = tol
  )
  rc <- rbind(c(300, 20, 0), c(1500, -40, 2), c(5000, 0, 10), c(-1, 0, 0))
  expect_equal(gplume(10, 3, 50, rc), ref$gp1, tolerance = tol)
  expect_equal(
    gplume(10, 3, 50, rc, stab = "A", set = "urban", mixh = 300, nimg = 2),
    ref$gp2,
    tolerance = tol
  )
  expect_equal(
    gspuff(5, 2, 10, rc, 300, stab = "B", mixh = 200), ref$gs1,
    tolerance = tol
  )
  c0 <- matrix(0, 7, 6)
  c0[3, 3] <- 1
  c0[4, 4] <- 2
  src <- matrix(0, 7, 6)
  src[5, 3] <- 0.5
  a <- advdif2d(c0, 0.7, -0.4, 0.2, 0.1, 1, 1, 0.5, 6, source = src)
  expect_equal(
    c(a$field, a$mass, a$cfl, a$difn, a$stabn), ref$ad1,
    tolerance = tol
  )
  expect_warning(
    a <- advdif2d(c0, 2, 1, 0.5, 0.5, 1, 1, 1, 3),
    "stability bound"
  )
  expect_equal(c(a$field, a$dt, a$substeps, a$stabn), ref$ad2, tolerance = tol)
  p <- lagrprt(
    5, 1, 2, 0.8, 0.1, 0.5, 0.3, 0.7, 4,
    seed = 5, grid = c(-2, 6, -2, 6, 4, 4)
  )
  expect_equal(
    c(p$x, p$y, p$mean_x, p$meany, p$concentration), ref$lg1,
    tolerance = tol
  )
  p <- lagrprt(3, 0, 0, 1, 1, 1, 1, 1, 2, seed = -7)
  expect_equal(c(p$x, p$y), ref$lg2, tolerance = tol)
  f <- plmfld(100, 5, 50, x = c(100, 600, 2000), y = c(-100, 0, 100))
  expect_equal(c(f$field, f$dx, f$dy), ref$pf1, tolerance = tol)
})

test_that("the health functions equal the reference implementation", {
  tol <- 1e-10
  expect_equal(
    unname(unlist(crfpm25(c(5, 12, 20))[
      c("rr", "logrr", "rrlo", "rrup", "expconc")
    ])),
    ref$cp1,
    tolerance = tol
  )
  expect_equal(
    crfpm25(c(8, 30), "ihd", refconc = 5)$extra$rrperunit, ref$cp2,
    tolerance = tol
  )
  expect_equal(
    unname(unlist(crfpm25(17, "stroke")[c("rr", "logrr")])), ref$cp3,
    tolerance = tol
  )
  expect_equal(
    unname(unlist(
      crfno2(c(15, 40), "chiast")[c("rr", "logrr", "rrlo", "rrup")]
    )),
    ref$cn1,
    tolerance = tol
  )
  expect_equal(
    unname(unlist(crfno2(33, betaper10 = -0.02)[c("rr", "logrr")])), ref$cn2,
    tolerance = tol
  )
  a <- attrfrac(1.08, 0.7, 1.06, 1.09)
  expect_equal(
    unname(c(as.numeric(a), attr(a, "ci"), attrfrac(0.8, 0.5))), ref$af1,
    tolerance = tol
  )
  expect_equal(
    c(
      mortdisp(5, 2.8e6, 0.0075, log(1.08) / 10),
      mortdisp(-3, 1e5, 5e-6, 0.01, rare = TRUE)
    ),
    ref$md1,
    tolerance = tol
  )
  b <- polbrdn(12, 0.9, 0.005, 5e5, poll = "PM2.5", y = "ihd")
  expect_equal(
    c(b$paf, b$attrcases, b$basecases, b$extra$rr), ref$pb1,
    tolerance = tol
  )
  b <- polbrdn(30, 1, 0.012, 5e5, y = "chiast", refconc = 8)
  expect_equal(c(b$paf, b$attrcases, b$extra$rr), ref$pb2, tolerance = tol)
  t <- data.frame(
    fsa = c("A", "B", "C"), exposure = c(18, 28, 9),
    population = c(60000, 40000, 5000), brate = c(0.008, 0.007, 0.01),
    stringsAsFactors = FALSE
  )
  pa <- polbrdarea(t)
  expect_identical(pa$fsa, c("B", "A", "C"))
  expect_equal(
    unname(unlist(pa[c("rr", "paf", "attrcases", "basecases")])), ref$pa1,
    tolerance = tol
  )
  e <- expconcidx(
    data.frame(h = c(30, 25, 25, 15, 10, 12), i = c(1, 2, 2, 3, 5, 4)),
    "h", "i"
  )
  expect_equal(
    c(e$concidx, e$expomean, e$extra$cov_h_R, e$extra$curve$exposure),
    ref$ec1,
    tolerance = tol
  )
})

test_that("pgsig follows the Briggs fits for rural and urban classes", {
  x <- c(100, 2500)
  r <- pgsig(x, "C")
  expect_equal(r$sigmay, 0.11 * x / sqrt(1 + 1e-4 * x), tolerance = 1e-12)
  expect_equal(r$sigmaz, 0.08 * x / sqrt(1 + 2e-4 * x), tolerance = 1e-12)
  expect_equal(
    pgsig(x, "A", "urban")$sigmaz, 0.24 * x * sqrt(1 + 1e-3 * x),
    tolerance = 1e-12
  )
  expect_equal(pgsig(x, "a")$sigmay, 0.22 * x / sqrt(1 + 1e-4 * x),
               tolerance = 1e-12)
  expect_error(pgsig(100, "G"), "Pasquill")
  expect_error(pgsig(100, c("A", "B")), "Pasquill")
  expect_error(pgsig(100, "D", "suburban"), "rural")
})

test_that("briggs: buoyant branches, stable branch, transitional rise", {
  F_ <- 9.80616 * 12 * (420 - 285) / (4 * 420)
  r <- briggs(c(10, 1e5), 4, diam = 1, exitvel = 12, stktmp = 420,
              ambtmp = 285)
  expect_equal(r$flux, F_, tolerance = 1e-8)
  expect_equal(
    r$rise, c(1.6 * F_^(1 / 3) * 10^(2 / 3) / 4, 21.425 * F_^0.75 / 4),
    tolerance = 1e-12
  )
  expect_equal(r$xfinal, 49 * F_^(5 / 8), tolerance = 1e-12)
  s <- 9.80616 * 0.02 / 285
  st <- briggs(1e5, 4, 1, 12, 420, 285, stab = "E")
  expect_equal(st$finalrise, 2.6 * (F_ / (4 * s))^(1 / 3), tolerance = 1e-12)
  expect_equal(st$xfinal, 2.0715 * 4 / sqrt(s), tolerance = 1e-12)
  expect_identical(briggs(10, 4, 1, 12, 285.2, 285.4)$flux, 0)
  expect_error(briggs(10, 0, 1, 12, 420, 285), "wind speed")
  expect_error(briggs(10, 4, 1, 12, 155, 20), "kelvin")
  expect_error(briggs(10, 4, 1, 12, 280, 290), "below `ambtmp`")
})

test_that("gplume and gspuff reproduce their closed forms", {
  s <- pgsig(800, "D")
  expect_equal(
    gplume(10, 3, 0, rbind(c(800, 0, 0)))[1],
    10 / (pi * 3 * s$sigmay * s$sigmaz),
    tolerance = 1e-14
  )
  expect_identical(gplume(10, 3, 0, list(c(-5, 0, 0), c(800, 0, 0)))[1], 0)
  expect_equal(
    gplume(10, 3, 0, c(800, 0, 0, 900, 0, 0)),
    gplume(10, 3, 0, rbind(c(800, 0, 0), c(900, 0, 0))),
    tolerance = 1e-15
  )
  expect_equal(
    gplume(10, 3, 0, data.frame(x = 800, y = 0, z = 0)),
    gplume(10, 3, 0, rbind(c(800, 0, 0))),
    tolerance = 1e-15
  )
  near <- gplume(10, 3, 50, rbind(c(200, 0, 0)))
  expect_equal(near, gplume(10, 3, 50, rbind(c(200, 0, 0)), mixh = 300),
               tolerance = 1e-12)
  expect_error(gplume(10, 3, 0, matrix(1:4, ncol = 2)), "rcpt")
  expect_error(gplume(10, 3, 0, c(1, 2, 3, 4)), "rcpt")
  expect_error(gplume(10, 3, 0, list(c(1, 2))), "rcpt")
  expect_error(gplume(10, 3, 0, rbind(c(NA, 0, 0))), "rcpt")
  expect_warning(gplume(10, 3, 0, rbind(c(800, 0, -1))), "below ground")
  expect_warning(
    gplume(10, 3, 0, rbind(c(800, 0, 400)), mixh = 300),
    "above the mixing height"
  )
  expect_warning(gspuff(1, 3, 400, rbind(c(8, 0, 0)), 1, mixh = 300),
                 "above the mixing height 300")
  expect_error(
    gplume(10, 3, 0, rbind(c(1, 0, 0)), nimg = 2000), "at most 1,000"
  )
})

test_that("advdif2d: advection shifts, diffusion spreads, source adds", {
  c0 <- matrix(0, 6, 6)
  c0[3, 4] <- 1
  expect_equal(advdif2d(c0, 1, 0, 0, 0, 1, 1, 1, 1)$field[4, 4], 1,
               tolerance = testthat_tolerance())
  d <- advdif2d(c0, 0, 0, 0.1, 0.2, 1, 1, 1, 1)
  expect_equal(c(d$field[3, 4], d$field[2, 4], d$field[3, 5]), c(0.4, 0.1, 0.2),
               tolerance = 1e-15)
  n <- advdif2d(c0, -1, 0, 0, 0, 1, 1, 1, 1)
  expect_equal(n$field[2, 4], 1, tolerance = testthat_tolerance())
  expect_equal(advdif2d(c0, 0, -1, 0, 0, 1, 1, 1, 1)$field[3, 3], 1,
               tolerance = testthat_tolerance())
  src <- matrix(0, 6, 6)
  src[4, 4] <- 2
  w <- advdif2d(matrix(0, 6, 6), 0, 0, 0, 0, 1, 1, 0.5, 1, source = src)
  expect_equal(w$mass, 1, tolerance = testthat_tolerance())
  expect_identical(advdif2d(c0, 1, 0, 0, 0, 1, 1, 1, 0)$field, c0)
  expect_error(advdif2d(c0, 0, 0, 0, 0, 1, 1, 1, 1, source = matrix(0, 2, 2)),
               "dimensions")
  expect_error(advdif2d(1:3, 0, 0, 0, 0, 1, 1, 1, 1), "numeric matrix")
  expect_error(advdif2d(c0 - 1, 0, 0, 0, 0, 1, 1, 1, 1), "non-negative")
  expect_error(advdif2d(c0 * NA, 0, 0, 0, 0, 1, 1, 1, 1), "finite")
  expect_error(
    advdif2d(c0, 0, 0, 0, 0, 1, 1, 1, 1, source = src - 3), "emission rate"
  )
  expect_error(
    advdif2d(c0, 0, 0, 0, 0, 1, 1, 1, 1, source = src * NA), "finite numeric"
  )
  expect_error(advdif2d(c0, 0, 0, 0, 0, 1, 1, 1, 1e9), "above maxwork")
  expect_error(
    suppressWarnings(advdif2d(c0, 9, 0, 0, 0, 1, 1, 1, 1e6)),
    "from a stability number"
  )
  expect_error(advdif2d(c0, 0, 0, 0, 0, 1, 1, 1, 1, maxwork = 0), "maxwork")
  expect_error(advdif2d(c0, 0, 0, 0, 0, 1, 1, 1, 3e9, maxwork = Inf), "2\\^31")
  expect_error(advdif2d(c0, 1e300, 0, 1e300, 0, 1e-300, 1, 1e300, 1),
               "not finite")
})

test_that("lagrprt: the random walk, the normals seam and the grid", {
  zeros <- function(n, seed, stream) rep(0, n)
  det <- lagrprt(3, 0, 0, 1, 2, 0.5, 0.5, 1, 4, normals = zeros)
  expect_equal(det$x, rep(4, 3), tolerance = testthat_tolerance())
  expect_equal(det$y, rep(8, 3), tolerance = testthat_tolerance())
  g <- lagrprt(10, 0.5, 0.5, 0, 0, 0, 0, 1, 1, grid = c(0, 1, 0, 1, 2, 2),
               normals = zeros)
  expect_equal(g$concentration[2, 2], 10 / (10 * 0.5 * 0.5),
               tolerance = testthat_tolerance())
  expect_equal(sum(g$concentration) * 0.25, 1, tolerance = testthat_tolerance())
  out <- lagrprt(4, 5, 5, 0, 0, 0, 0, 1, 1, grid = c(0, 1, 0, 1, 2, 2))
  expect_equal(sum(out$concentration), 0, tolerance = testthat_tolerance())
  # the default stream's normals are standard normal
  p <- lagrprt(20000, 0, 0, 0, 0, 0.5, 0.5, 1, 1, seed = 1)
  expect_lt(abs(mean(p$x)), 0.03)
  expect_lt(abs(stats::sd(p$x) - 1), 0.03)
  expect_error(lagrprt(0, 0, 0, 1, 1, 1, 1, 1, 1), "positive")
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, grid = 1:3), "grid")
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, grid = c(1, 0, 0, 1, 1, 1)),
               "xmin < xmax")
  expect_error(
    lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, grid = c(0, 1, 0, 1, 1e5, 1e5)),
    "cells, above maxwork"
  )
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, seed = 0.5), "seed")
  expect_error(lagrprt(2e8, 0, 0, 1, 1, 1, 1, 1, 1), "100 million")
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1e9), "above maxwork")
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, maxwork = -1), "maxwork")
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 3e9, maxwork = Inf), "2\\^31")
  expect_error(lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, normals = 3), "a function")
  expect_error(
    lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, normals = function(n, s, k) stop("x")),
    "failed when called"
  )
  expect_error(
    lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, normals = function(n, s, k) 1),
    "1, 1 of them finite"
  )
  expect_error(
    lagrprt(2, 0, 0, 1, 1, 1, 1, 1, 1, normals = function(n, s, k) "a"),
    "a character"
  )
})

test_that("dispersion inputs are refused with the argument named", {
  expect_error(pgsig(TRUE), "pgsig: `x` must be numeric, not logical")
  expect_error(pgsig("a"), "got character")
  expect_error(pgsig(NULL), "got NULL")
  expect_error(gplume(c(1, 2), 1, 1, c(1, 0, 0)), "single number")
  expect_error(pgsig(NA_real_), "finite")
  expect_error(pgsig(-1), "at least 0")
  expect_error(pgsig(1e301), "below 1e300")
  expect_error(advdif2d(matrix(0, 3, 3), 0, 0, 0, 0, 1, 1, 1, 1.5), "whole")
  expect_error(plmfld(1, 1, 1, c(3, 1), 1:2), "strictly increasing")
})

test_that("CRFs: 1 at the reference, the log-linear and IER forms", {
  expect_identical(crfpm25(5.8)$rr, 1)
  r <- crfpm25(12)
  expect_identical(r$rr, exp(log(1.08) * (12 - 5.8) / 10))
  expect_null(r$extra$rrperunit)
  # pow() and exp() differ in the last digit between systems
  expect_equal(
    crfpm25(12, y = "ihd")$rr, 1 + 1.91 * (1 - exp(-0.14 * (12 - 5.8)^0.49)),
    tolerance = 1e-14
  )
  expect_true(is.na(crfpm25(12, y = "ihd")$rrlo))
  expect_identical(crfno2(10)$rr, 1)
  expect_identical(crfno2(30)$logrr, 2 * crfno2(20)$logrr)
  # respiratory deaths: RR 1.03 (1.00 to 1.05) per 10 ug/m3 above 10
  rs <- crfno2(25, y = "respiratory")
  expect_identical(rs$rr, exp(log(1.03) * (25 - 10) / 10))
  expect_identical(c(rs$rrlo, rs$rrup), exp(log(c(1, 1.05)) * 15 / 10))
  expect_equal(crfno2(20, y = "anything", betaper10 = 0.05)$logrr, 0.05)
  expect_match(crfno2(20, y = "chiast")$citation, "Khreis")
  expect_error(crfpm25(10, y = "nope"), "Unknown outcome")
  expect_error(crfno2(20, y = "nope"), "or pass betaper10")
  expect_error(crfno2(20, betaper10 = c(1, 2)), "betaper10")
  expect_error(crfpm25("twelve"), "`expo` must be numeric, got character")
  expect_error(crfpm25(TRUE), "not logical")
  expect_error(crfpm25(NULL), "got NULL")
  expect_error(crfpm25(numeric(0)), "finite concentrations")
  expect_error(crfpm25(-1), "non-negative")
})

test_that("attrfrac, mortdisp and polbrdn by hand and refused", {
  expect_identical(as.numeric(attrfrac(2, 1)), 1 - 1 / 2)
  expect_identical(as.numeric(attrfrac(1.5, 0.4)), 0.4 * 0.5 / (1 + 0.4 * 0.5))
  expect_error(attrfrac(0, 1), "positive")
  expect_error(attrfrac(1.5, 1.2), "must be in")
  expect_error(attrfrac(c(1, 2), 0.5), "single")
  expect_error(attrfrac(1.5, 0.5, rrlo = 1), "go together")
  expect_error(attrfrac(1.5, 0.5, rrlo = 1.6, rrup = 2), "between rrlo")
  expect_error(attrfrac(factor(2), 0.5), "not a factor")
  expect_error(attrfrac(list(2), 0.5), "got list")
  expect_identical(
    mortdisp(10, 1e6, 0.008, 0.0039), 0.008 * 1e6 * (1 - exp(-0.039))
  )
  expect_identical(mortdisp(1, 1e6, 0, 0.1), 0)
  expect_error(mortdisp(1, -1, 0.008, 0.001), "non-negative")
  expect_error(mortdisp(1, 1, 800, 0.001), "looks like a rate per 100,000")
  expect_error(mortdisp(1, 1, 1e-12, 0.001), "100 billion")
  expect_error(mortdisp(1, 1, 1e-6, 0.001), "rare = TRUE")
  expect_error(mortdisp(1, 1, 0.5, 0.001), "50,000 per 100,000")
  expect_error(mortdisp(1, 1, 0.5, 0.001, hirate = NA), "TRUE or FALSE")
  b <- polbrdn(25, 1, 0.008, 1e6, poll = "NO2")
  rr <- exp(log(1.02) * 1.5)
  expect_identical(b$extra$rr, rr)
  expect_identical(b$basecases, 8000)
  expect_identical(polbrdn(4, 1, 0.008, 1e6, poll = "PM2.5")$attrcases, 0)
  expect_identical(polbrdn(30, 1, 0.012, 5e5, y = "chiast")$extra$unit,
                   "incident cases")
  expect_error(polbrdn(4, 1, 0.008, 1e6, poll = "O3"), "Unknown pollutant")
  expect_error(polbrdn(4, 2, 0.008, 1e6), "must be in")
})

test_that("polbrdarea sorts worst first and refuses bad tables", {
  t <- data.frame(fsa = c("A", "B"), exposure = c(18, 28),
                  population = c(4e4, 6e4), brate = 0.008,
                  stringsAsFactors = FALSE)
  expect_identical(polbrdarea(t)$fsa, c("B", "A"))
  expect_s3_class(polbrdarea(t), "polbrdarea")
  expect_error(polbrdarea(t[, 1:2]), "missing columns")
  expect_error(polbrdarea(t[0, ]), "no rows")
  t2 <- t
  t2$exposure <- as.character(t2$exposure)
  expect_error(polbrdarea(t2), "must be numeric")
  t3 <- t
  t3$brate <- 800
  expect_error(polbrdarea(t3),
               paste("polbrdarea: the `basratecol` column (`brate`) of row 1",
                     "is a rate"), fixed = TRUE)
  expect_error(polbrdarea(t, areacol = 1), "`areacol` must name one column")
})

test_that("expconcidx: zero, pro-poor, pro-rich, and its refusals", {
  e <- expconcidx(data.frame(e = 20, i = c(1, 1, 2, 3)), "e", "i")
  expect_equal(e$concidx, 0, tolerance = testthat_tolerance())
  expect_match(e$interpretation, "evenly")
  expect_match(
    expconcidx(data.frame(e = 5:1, i = 1:5), "e", "i")$interpretation,
    "Pro-poor"
  )
  expect_match(
    expconcidx(data.frame(e = 1:5, i = 1:5), "e", "i")$interpretation,
    "Pro-rich"
  )
  expect_error(expconcidx(data.frame(e = 1:3), "e", "i"), "missing columns")
  expect_error(expconcidx(data.frame(e = 1, i = 1), "e", "i"), "at least 2")
  expect_error(expconcidx(data.frame(e = 1:3, i = c("a", "b", "c"),
                                     stringsAsFactors = FALSE), "e", "i"),
               "numeric")
  expect_error(expconcidx(data.frame(e = c(-1, 2), i = 1:2), "e", "i"),
               "non-negative")
  expect_error(expconcidx(data.frame(e = 1:2, i = 1), "e", "i"), "two values")
  expect_error(expconcidx(data.frame(e = 0, i = 1:2), "e", "i"), "zero")
  expect_error(expconcidx(data.frame(e = 1), c("e", "f"), "i"), "one column")
})

test_that("an exposure series joined to daily crime counts runs as it is", {
  # simulated: daily NO2 and assaults for three areas
  set.seed(1)
  d <- data.frame(area = rep(c("a", "b", "c"), each = 30), day = rep(1:30, 3),
                  stringsAsFactors = FALSE)
  d$no2 <- stats::rlnorm(90, log(20), 0.3)
  d$assaults <- stats::rpois(90, 3 + 0.05 * d$no2)
  rr <- crfno2(d$no2)$extra$rrperunit
  expect_length(rr, 90L)
  expect_true(all(rr >= 1))
  expect_identical(rr > 1, d$no2 > 10)
})

test_that("the air pollution results print", {
  c0 <- matrix(0, 5, 5)
  c0[3, 3] <- 1
  out <- c(
    capture.output(print(
      suppressWarnings(advdif2d(c0, 2, 0, 0, 0, 1, 1, 1, 1))
    )),
    capture.output(print(plmfld(100, 5, 50, x = c(100, 600), y = 0))),
    capture.output(print(lagrprt(5, 0, 0, 1, 0, 1, 1, 1, 2,
                                 grid = c(-5, 5, -5, 5, 2, 2)))),
    capture.output(print(crfpm25(c(12, 20)))),
    capture.output(print(attrfrac(0.8, 0.5, 0.7, 0.9))),
    capture.output(print(polbrdn(25, 1, 0.008, 1e6))),
    capture.output(print(polbrdarea(data.frame(
      fsa = "A", exposure = 20, population = 1, brate = 0.008,
      stringsAsFactors = FALSE
    )))),
    capture.output(print(expconcidx(data.frame(e = 5:1, i = 1:5), "e", "i")))
  )
  expect_true(any(grepl("sub-steps", out, fixed = TRUE)))
  expect_true(any(grepl("receptors; peak", out, fixed = TRUE)))
  expect_true(any(grepl("peak density", out, fixed = TRUE)))
  expect_true(any(grepl("their mean", out, fixed = TRUE)))
  expect_true(any(grepl("protects", out, fixed = TRUE)))
  expect_true(any(grepl("attributable", out, fixed = TRUE)))
  expect_true(any(grepl("by area: 1 areas", out, fixed = TRUE)))
  expect_true(any(grepl("Concentration index", out, fixed = TRUE)))
})

test_that("the stage functions refuse a logical or NULL number", {
  expect_error(attrfrac(TRUE, 1),
               "attrfrac: `rr` must be numeric, not logical.", fixed = TRUE)
  expect_error(attrfrac(NULL, 1), "attrfrac: `rr` must be numeric, got NULL.",
               fixed = TRUE)
})

test_that("plmfld() refuses in its own name", {
  expect_error(plmfld(q = NULL, u = 5, h = 50, x = seq(100, 300, 100),
                      y = c(-50, 50)), "^plmfld: `q` must be numeric")
})

test_that("an outcome or a class that is not text is refused in words", {
  expect_error(crfno2(20, y = sum), "^crfno2: Unknown outcome 'function'")
  expect_error(gplume(1, 5, 50, rbind(c(100, 0, 0)), stab = sum),
               "^gplume: `stab` must be one of the Pasquill classes")
  bad <- rawToChar(as.raw(c(0xff, 0xfe, 0x41)))
  expect_error(gplume(1, 5, 50, rbind(c(100, 0, 0)), stab = bad),
               "^gplume: `stab` holds text that is not valid")
})

test_that("advdif2d() refuses a source that is not numbers", {
  c0 <- matrix(0, 5, 5)
  c0[3, 3] <- 1
  for (s in list(sum, new.env())) {
    expect_error(advdif2d(c0, 0.7, -0.4, 0.2, 0.1, 1, 1, 0.5, 6, source = s),
                 "^advdif2d: `source` must be a finite numeric matrix\\.$")
  }
})

test_that("polbrdarea() and polbrdn() name the column, row and function", {
  t <- data.frame(fsa = c("A1", "B2"), exposure = c(NA, 18),
                  population = c(4e4, 6e4), brate = 0.008)
  expect_error(polbrdarea(t),
               paste("polbrdarea: the `expcol` column (`exposure`) of row 1",
                     "must be"), fixed = TRUE)
  expect_error(polbrdn(25, 1, 0.008, 1e6, poll = "PM2.5", y = "chiast"),
               "^polbrdn: Unknown outcome 'chiast'")
  expect_error(polbrdn(25, 1, 0.008, 1e6, poll = "NO2", y = "nope"),
               "^polbrdn: Unknown outcome 'nope'")
})
