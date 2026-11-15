test_that("estimates with intervals", {
  d <- did()
  f <- didfe(d$rate, d$town, d$year, d$bylaw, bacon = TRUE, loo = TRUE)
  p <- drw(plot(f))
  expect_identical(p$coef$est, f$coef$est)
  expect_identical(p$bacon, f$bacon)
  z <- stats::qnorm(0.975)
  expect_identical(p$loo$hi, f$loo$est + z * f$loo$se)
  expect_named(drw(plot(didfe(d$rate, d$town, d$year, d$bylaw))), "coef")
  pr <- data.frame(EndFiscalYear = rep(2020:2022, each = 2),
                   Gender = rep(c("Female", "Male"), 3),
                   NumberIndividuals_RestrictiveConfinement =
                     c(40, 310, 52, 335, 47, 360), stringsAsFactors = FALSE)
  cr <- cntreg(clprep(pr, "otis-c01-cnt"))
  pc <- drw(plot(cr))
  expect_identical(pc$est, cr$irr$irr)
  expect_identical(pc$group, cr$irr$fit)
  ds <- data.frame(age = rep(c("18-29", "30-49", "50+"), 2),
                   grp = rep(c("A", "B"), each = 3),
                   x = c(30, 22, 5, 41, 30, 12),
                   n = c(5000, 9000, 6000, 3000, 8000, 7000),
                   stringsAsFactors = FALSE)
  r <- dastd(ds$x, ds$n, ds$age, ds$grp, std = c(30, 45, 25), per = 1e4)
  expect_identical(drw(plot(r))$est, r$rates$rate)
  g <- grprr(c(a = 120, b = 95, c = 30), c(a = 6e4, b = 2.1e4, c = 1.5e4))
  expect_identical(drw(plot(g))$lo, g$rates$lo)
  gp <- grprr(c(1, 0, 0, 1, 1, 0, 1, 1), grp = c(1, 1, 1, 1, 2, 2, 2, 2),
              type = "prop")
  expect_identical(drw(plot(gp))$est, gp$rates$rate)
  si <- sirex(ds$x, ds$n, ds$age, ds$grp)
  expect_identical(drw(plot(si))$est, si$sir$sir)
  set.seed(1)
  sv <- data.frame(region = rep(c("north", "south"), each = 60),
                   block = rep(1:12, each = 10),
                   victim = stats::rbinom(120, 1, 0.3),
                   w = stats::runif(120, 50, 150), stringsAsFactors = FALSE)
  s <- svymx(sv$victim, sv$w, sv$region, sv$block, stat = "prop")
  expect_identical(drw(plot(s))$hi, s$est$hi)
  set.seed(6)
  md <- data.frame(id = sample(1:80, 150, TRUE),
                   yr = sample(2021:2022, 150, TRUE),
                   days = stats::rgeom(150, 0.08) + 1)
  m <- mndla(md$days, md$id, md$yr)
  expect_identical(drw(plot(m))$est, m$tab$share)
  td <- data.frame(grp = rep(c("a", "b"), 4), area = rep(1:4, each = 2),
                   stops = c(1200, 800, 900, 1100, 1500, 700, 1000, 1000),
                   srch = c(60, 70, 50, 95, 80, 66, 52, 90),
                   hits = c(24, 18, 19, 25, 30, 17, 21, 24),
                   stringsAsFactors = FALSE)
  th <- thrtst(td$stops, td$srch, td$hits, td$grp, area = td$area)
  expect_identical(drw(plot(th))$est, th$thr$thr)
  set.seed(10)
  pri <- stats::rbinom(300, 1, 0.4)
  tt <- ifelse(stats::rbinom(300, 1, 0.5) == 1,
               stats::rweibull(300, 1.3, 300), Inf)
  ce <- stats::runif(300, 365, 1500)
  sc <- spcure(pmin(tt, ce), as.numeric(tt <= ce), data.frame(prior = pri))
  ps <- drw(plot(sc))
  expect_identical(ps$est, sc$coef$est)
  expect_identical(ps$group, sc$coef$part)
  set.seed(13)
  tz <- factor(sample(c("fine", "custody"), 200, TRUE))
  ty <- stats::rbinom(200, 1, 0.3)
  mb <- mtrbnd(ty, tz, assume = "both")
  expect_identical(drw(plot(mb))$hi, mb$bnd$upper)
  pb <- pbbnd(c("1,240", "35", "<5", "x"), rnd = 5, suplim = 10)
  pp <- drw(plot(pb))
  expect_identical(pp$lo, pb$bnd$lower)
  expect_identical(pp$group, pb$bnd$status)
  dk <- dkfig(v = c(0.12, 0.08), r = 0.03, amax = 0.01, bmax = 0.2,
              thr = 0.15)
  expect_identical(drw(plot(dk))$hi, dk$bnd$hi)
  set.seed(1)
  bx <- data.frame(age = stats::rnorm(200, 40, 10),
                   prior = stats::rpois(200, 2))
  bd <- stats::rbinom(200, 1, stats::plogis((bx$age - 40) / 10))
  sb <- smdbal(bx, bd)
  ss <- drw(plot(sb))
  expect_identical(sort(ss$est), sort(abs(sb$tab$smd)))
  set.seed(9)
  j <- sample(paste0("J", 1:4), 400, TRUE)
  jm <- jmono(j, stats::rbinom(400, 1, 0.4))
  expect_identical(sort(drw(plot(jm))$est), sort(jm$judges$rate))
  tab <- matrix(c(20, 10, 80, 90), 2,
                dimnames = list(c("A", "B"), c("yes", "no")))
  rc <- relchk(tab, 2.25, ref = "B", event = "yes")
  pr2 <- drw(plot(rc))
  expect_identical(pr2$est, rc$reads$value)
  expect_identical(pr2$group == "reproduces the claim", rc$reads$match)
  dx <- dissix(cbind(a = c(90, 80, 10, 5, 50), b = c(10, 20, 85, 95, 50)))
  pd <- drw(plot(dx))
  expect_identical(pd$est[pd$group == "Gini"], dx$grp$gini)
})

test_that("spatial regression coefficients", {
  set.seed(1)
  p <- cbind(x = stats::runif(60), y = stats::runif(60))
  w <- spwt(p, k = 5)
  x <- stats::rnorm(60)
  f <- sarml(1 + 2 * x + stats::rnorm(60), x, w)
  s <- drw(plot(f))
  expect_identical(s$est, f$coef$est)
  expect_equal(s$hi - s$est, stats::qnorm(0.975) * f$coef$se,
               tolerance = 1e-12)
})

test_that("paths over time", {
  d <- did()
  one <- d[d$town %in% c("A", LETTERS[9:12]), ]
  sy <- synthc(one$rate, one$town, one$year, one$first)
  expect_identical(drw(plot(sy)), sy$path)
  sd <- sdidx(one$rate, one$town, one$year, one$first)
  expect_identical(drw(plot(sd)), sd$path)
  tr <- trndrr(c(402, 377, 390, 368, 331, 310, 305), 2017:2023)
  pt <- drw(plot(tr))
  expect_identical(pt$fitted, tr$fitted)
  sk <- stkflw(people = c(5200, 4800, 4300), days = c(310000, 318000, 330000),
               prd = 2021:2023)
  expect_identical(drw(plot(sk, main = "Stock"))$adp, sk$flow$adp)
  dk <- dktflw(recv = c(120, 130, 125, 140), disp = c(110, 135, 120, 150),
               pend = c(400, 395, 400, 390), win = 2, span = 91)
  expect_identical(drw(plot(dk))$clear, dk$flow$clear)
  jp <- jpipe(num = c(5, 7, 6, 9, 8, 10), den = c(20, 22, 19, 25, 24, 26),
              win = 3, lag = 1)
  expect_identical(drw(plot(jp))$rate, jp$rate)
  set.seed(1)
  hp <- cbind(x = c(stats::rnorm(150, 2, 0.3), stats::runif(150, 0, 10)),
              y = c(stats::rnorm(150, 2, 0.3), stats::runif(150, 0, 10)))
  hs <- hsper(hp, prd = rep(2021:2023, length.out = 300), cell = 1,
              top = 0.1)
  expect_identical(drw(plot(hs))$hotshr, hs$periods$hotshr)
  yd <- data.frame(year = rep(2019:2023, 2), svc = rep(c("A", "B"), each = 5),
                   n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33),
                   stringsAsFactors = FALSE)
  yy <- yoypx(yd$n, yd$year, yd$svc, rnd = 5)
  py <- drw(plot(yy))
  expect_identical(py$pct, yy$pctchg[is.finite(yy$pctchg)])
  cu <- cumrsk(18:22, c(120, 100, 80, 60, 40), rep(10000, 5))
  expect_identical(drw(plot(cu))$hi, cu$hi)
  ages <- 12:30
  am <- agemix(cbind(al = stats::dnorm(ages, 17, 2),
                     lcp = rep(0.08, 19)), c(0.8, 0.2), ages)
  pa <- drw(plot(am))
  expect_identical(pa$y[pa$series == "agg"], am$curve$agg)
  set.seed(11)
  y <- matrix(stats::rpois(20, 50), 5,
              dimnames = list(c("15-19", "20-24", "25-29", "30-34", "35-39"),
                              c(2000, 2005, 2010, 2015)))
  ap <- apcie(y, matrix(1e4, 5, 4))
  pp <- drw(plot(ap, main = "APC"))
  expect_identical(pp$est, ap$coef$est[ap$coef$term != "(Intercept)"])
})

test_that("curves", {
  km <- kmrec(c(3, 8, 8, 12, 20, 24), event = c(1, 1, 0, 1, 0, 1),
              grp = c("a", "a", "a", "b", "b", "b"))
  pk <- drw(plot(km))
  expect_identical(pk$surv[pk$time > 0], km$km$surv)
  lf <- lftab(c(0, 3, 6, 12, 24, Inf), events = c(40, 30, 20, 10, 5),
              cens = c(2, 2, 3, 3, 4))
  pl <- drw(plot(lf))
  expect_identical(pl$surv, c(1, lf$tab$surv[1:4]))
  set.seed(3)
  n <- 200
  tp <- rep(0:1, each = n / 2)
  ch <- c("suspect", "filler", "none")
  pick <- ifelse(tp == 1, sample(ch, n, TRUE, c(0.5, 0.15, 0.35)),
                 sample(ch, n, TRUE, c(0.1, 0.3, 0.6)))
  conf <- ifelse(pick == "none", NA, sample(1:5, n, TRUE))
  lr <- lnroc(tp, pick, conf, proc = rep(c("sim", "seq"), n / 2), nboot = 0)
  pr <- drw(plot(lr))
  expect_identical(sort(pr$hr[pr$hr > 0]), sort(lr$roc$hr[lr$roc$hr > 0]))
  set.seed(9)
  s <- rep(1:6, each = 5)
  cx <- stats::runif(6, 0, 5000)[s] + stats::rnorm(30, 0, 1500)
  cy <- stats::runif(6, 0, 5000)[s] + stats::rnorm(30, 0, 1500)
  day <- stats::runif(6, 0, 300)[s] + stats::rnorm(30, 0, 60)
  lk <- lnkpr(cbind(x = cx, y = cy), t = day, series = s)
  pk2 <- drw(plot(lk))
  l <- as.logical(lk$pairs$linked)
  expect_identical(utils::tail(pk2$tpr, 1), 1)
  expect_identical(pk2$fpr[2L], mean(lk$pairs$prob[!l] >= max(lk$pairs$prob)))
  lk$pairs$linked <- TRUE
  expect_error(drw(plot(lk)), "^plot: an ROC curve needs both")
  set.seed(2)
  cc <- cconc(stats::rpois(100, 2), nsim = 0)
  pc <- drw(plot(cc))
  expect_identical(pc$share, cc$lorenz$share)
  expect_identical(attr(pc, "gini"), cc$stat[["gini"]])
  set.seed(1)
  p <- cbind(x = stats::runif(80), y = stats::runif(80))
  ck <- crossk(p, marks = ifelse(p[, 1] < 0.5, "a", "b"), i = "a", j = "b",
               win = c(0, 1, 0, 1), r = c(0.05, 0.1, 0.2), sigma = 0.1,
               nsim = 0)
  expect_identical(drw(plot(ck)), ck$k)
  bx <- bandx(c("1-5", "6-15", "16-30", "31-60", "Over 60"),
              c(40, 25, 15, 12, 8), stat = mean)
  expect_identical(drw(plot(bx)), bx$sens)
  set.seed(7)
  z <- stats::rbinom(300, 1, 0.5)
  dd <- stats::rbinom(300, 1, ifelse(z == 1, 0.3, 0.5))
  jg <- jdgiv(z, dd, stats::rbinom(300, 1, 0.3))
  expect_identical(drw(plot(jg)), jg$sens)
  a <- c("2019" = 120, "2020" = 98, "2021" = 105, "2022" = 130, "2023" = 141)
  b <- c("2018" = 50, "2019" = 61, "2020" = 44, "2021" = 52, "2022" = 66,
         "2023" = 70)
  xc <- xcser(a, b)
  expect_identical(drw(plot(xc)), xc$ccf)
  hx <- hwksim(c(a0 = log(0.5), eta = 0.5, beta = 2), 100, seed = 1)
  hr <- hwkres(hwkfit(hx$t, 100))
  ph <- drw(plot(hr))
  expect_identical(ph$u, sort(hr$u))
  expect_equal(ph$hi - ph$q, rep(1.36 / sqrt(hr$n), hr$n), tolerance = 1e-12)
})

test_that("bars", {
  aw <- arwt(c(a = 120, b = 40, c = 60, d = 100),
             from = c("a", "b", "c", "d", "d"),
             to = c("N1", "N1", "N2", "N2", "N3"), area = c(1, 1, 1, 7, 3))
  expect_identical(drw(plot(aw))$value, aw$out$value)
  cs <- crsev(c(10, 2, 8, 3), rep(c("theft", "robbery"), 2),
              c(theft = 1, robbery = 50), unit = rep(c("A", "B"), each = 2),
              pop = rep(1000, 4), base = 1)
  expect_identical(drw(plot(cs))$value, cs$tab$idx)
  peaks <- data.frame(locus = rep(c("D3", "vWA"), c(4, 3)),
                      allele = c(14, 15, 16, 17, 16, 17, 18),
                      height = c(1210, 380, 1150, 420, 1400, 1960, 610),
                      stringsAsFactors = FALSE)
  freq <- data.frame(locus = rep(c("D3", "vWA"), c(5, 5)),
                     allele = c(13, 14, 15, 16, 17, 15, 16, 17, 18, 19),
                     freq = c(0.05, 0.1, 0.25, 0.25, 0.2,
                              0.1, 0.2, 0.3, 0.2, 0.1),
                     stringsAsFactors = FALSE)
  ref <- data.frame(id = rep(c("V", "S"), each = 2),
                    locus = rep(c("D3", "vWA"), 2),
                    a1 = c(14, 16, 15, 17), a2 = c(16, 17, 17, 18),
                    stringsAsFactors = FALSE)
  gm <- gmixlr(peaks, freq, ref, hp = c("V", "S"), nc = 2, hd = "V")
  expect_identical(drw(plot(gm))$value, gm$locus$log10lr)
  trio <- data.frame(id = c("AF", "M", "C"), fid = c(NA, NA, "AF"),
                     mid = c(NA, NA, "M"), stringsAsFactors = FALSE)
  alt <- data.frame(id = c("AF", "M", "C", "X"), fid = c(NA, NA, "X", NA),
                    mid = c(NA, NA, "M", NA), stringsAsFactors = FALSE)
  geno <- data.frame(id = rep(c("AF", "M", "C"), 2),
                     locus = rep(c("D3", "vWA"), each = 3),
                     a1 = c(14, 15, 15, 16, 17, 17),
                     a2 = c(16, 17, 16, 18, 17, 19), stringsAsFactors = FALSE)
  pf <- data.frame(locus = rep(c("D3", "vWA"), c(4, 4)),
                   allele = c(14, 15, 16, 17, 16, 17, 18, 19),
                   freq = c(0.1, 0.25, 0.25, 0.2, 0.2, 0.3, 0.2, 0.1),
                   stringsAsFactors = FALSE)
  pe <- pedlr(trio, geno, pf, alt)
  expect_identical(drw(plot(pe))$value, log10(pe$locus$lr))
  pref <- rbind(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1), c(3, 1, 2))
  tt <- matrix(c(2, 5, 8, 6, 3, 7, 9, 6, 2, 7, 8, 3), 3)
  hq <- hypq(c(0.4, 0.7, 0.3, 0.6), pref, mu = 1.5, tt = tt)
  expect_identical(drw(plot(hq))$value, hq$wl)
  ch <- c(rep(3, 15), rep(1, 5), rep(2, 6), rep(4, 4), rep(5, 7), rep(6, 3))
  ls <- lnsize(ch, k = 6, susp = 3)
  pl <- drw(plot(ls))
  expect_identical(pl$value, ls$count)
  expect_identical(pl$label[3L], "member 3 (suspect)")
  set.seed(3)
  st <- lostay(round(stats::rlnorm(200, 2.3, 1.1)) + 1,
               breaks = c(0, 5, 15, Inf))
  pst <- drw(plot(st))
  expect_identical(pst$value, st$bands$share)
  expect_identical(pst$label[3L], "15 or more")
  sq <- function(x0, id) {
    data.frame(x = x0 + c(0, 1, 1, 0), y = c(0, 0, 1, 1), id = id)
  }
  pp <- pinpl(cbind(x = c(0.5, 0.2, 2.5), y = c(0.5, 0.7, 0.5)),
              poly = rbind(sq(0, "A"), sq(1, "B"), sq(2, "C")),
              pop = c(1200, 800, 450))
  expect_identical(drw(plot(pp))$value, pp$regs$units)
  mm <- mmshr(c(365, 400, 300, 365, 730), 365, c("a", "a", "b", "b", "b"))
  pm <- drw(plot(mm))
  expect_identical(pm$above, mm$above)
  set.seed(8)
  rv <- rptvic(sample(1:200, 150, TRUE), pop = 500)
  prv <- drw(plot(rv))
  expect_identical(prv$expected, rv$tab$expected)
})

test_that("matrices, surfaces and points", {
  ct <- ctabx(c("none", "none", "knife", "gun", "none", "knife"),
              c("hands", "hands", "CEW", "firearm", "CEW", "hands"))
  pc <- drw(plot(ct))
  expect_identical(pc$value, as.vector(ct$res))
  set.seed(1)
  p <- cbind(x = stats::runif(100, 0, 2000), y = stats::runif(100, 0, 2000))
  kn <- knxnr(p, t = stats::runif(100, 0, 365), ds = c(200, 400),
              dt = c(14, 28), nsim = 19, seed = 2)
  pk <- drw(plot(kn))
  expect_identical(sort(pk$value), sort(kn$tab$ratio))
  q <- qdrtst(p, win = c(0, 2000, 0, 2000), nx = 3, ny = 2)
  expect_identical(drw(plot(q))$value, as.vector(q$counts))
  tm <- trnmat(c("N", "N", "E", "N", "W", "W", "E", "E"),
               c(1, 1, 1, 2, 2, 2, 3, 3), time = c(1, 2, 3, 1, 2, 3, 1, 2))
  expect_identical(drw(plot(tm, pal = "grey"))$value, as.vector(tm$prob))
  zc <- zcomp(c(5, 9, 1, 30, 22, 41), rep(c("a", "b", "c"), 2),
              rep(c("assault", "theft"), each = 3))
  expect_identical(drw(plot(zc))$value, as.vector(zc$z))
  kd <- kdsurf(p, win = c(0, 2000, 0, 2000), nx = 20)
  expect_identical(drw(plot(kd, pal = "grey"))$z, kd$z)
  ph <- prohs(p, t = sort(sample(1:120, 100, TRUE)), at = 100,
              win = c(0, 2000, 0, 2000), cell = 200, bw = 400, ahead = 14)
  pp <- drw(plot(ph))
  expect_identical(pp$z, ph$z)
  expect_identical(pp$hot, ph$hot)
  sb <- sdbsim(8, 10, a0 = 1 / 30, omega = 1 / 15, eta = 0.03, gam = 0.02,
               theta = 0.5, seed = 1)
  expect_identical(drw(plot(sb))$z, sb$b)
  bf <- bufas(p, site = cbind(c(500, 1500), c(600, 1500)),
              comp = cbind(c(1000, 1800), c(300, 1900)), r = 200)
  expect_identical(drw(plot(bf)), bf$pts)
  set.seed(1)
  a <- cumsum(stats::rnorm(600, sd = 0.5)) + sin(seq_len(600) / 9)
  cm <- cmpsc(a, a[31:600] + stats::rnorm(570, sd = 0.1))
  expect_identical(drw(plot(cm)), cm$segs)
  eb <- ebshr(c(40, 12, 9, 25, 7, 31, 5, 18))
  expect_identical(drw(plot(eb)), eb$est)
  set.seed(2)
  off <- sample(sprintf("o%02d", 1:6), 1500, TRUE)
  night <- stats::rbinom(1500, 1, 0.4)
  yb <- stats::rbinom(1500, 1, stats::plogis(-1 + night + (off == "o03")))
  ib <- ibench(yb, off, data.frame(night), minn = 50)
  pi <- drw(plot(ib))
  expect_identical(pi$rate, ib$tab$rate)
  set.seed(1)
  kp <- cbind(x = c(stats::runif(80, 0, 20), stats::rnorm(15, 5, 0.3)),
              y = c(stats::runif(80, 0, 20), stats::rnorm(15, 15, 0.3)))
  ks <- kscan(kp, t = c(sample(0:99, 80, TRUE), sample(40:44, 15, TRUE)),
              r = c(1, 3), tw = c(5, 20), top = 1, nsim = 19, seed = 1)
  expect_identical(drw(plot(ks)), ks$cl)
})

test_that("a result missing its parts is refused, with the plot's message", {
  s3 <- getNamespaceInfo("crimlab", "S3methods")
  cls <- s3[s3[, 1L] == "plot", 2L]
  expect_gt(length(cls), 50L)
  for (k in cls) {
    expect_error(drw(plot(structure(list(), class = k))),
                 "^plot: `x` is not a complete", info = k)
  }
  e <- evtst(did()$rate, did()$town, did()$year, did()$first, leads = 2,
             lags = 2)
  expect_error(drw(plot(e, main = 1)), "^plot: `main` must be one string")
  expect_error(drw(plot(e, col = "nocolour")), "^plot: `col` must be colours")
  expect_error(drw(plot(e, pal = "neon")), "^plot: the palette must be")
})
