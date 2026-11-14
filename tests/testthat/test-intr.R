# Interrupts: every long loop of the kernels polls R's interrupt check, and
# an interrupt there ends the call as R's own interrupt, with no result.

# TRUE when evaluating e is interrupted with every poll of the kernels made
# to read as an interrupt (C_clab_r_intrt); R's check is back afterwards.
stops <- function(e) {
  ns <- asNamespace("crimlab")
  .Call(ns$C_clab_r_intrt, TRUE)
  on.exit(.Call(ns$C_clab_r_intrt, FALSE))
  tryCatch(
    {
      force(e)
      FALSE
    },
    interrupt = function(i) TRUE
  )
}

test_that("each long loop stops at an interrupt", {
  set.seed(1)
  p <- cbind(x = runif(60), y = runif(60))
  v <- p[, 1] + rnorm(60, sd = 0.2)
  w <- spwt(p, k = 4)
  x <- matrix(rnorm(400), 100)
  yc <- x[, 1] + rnorm(100)
  yb <- rbinom(100, 1, plogis(x[, 1]))
  fold <- rep(1:2, 50)
  g <- expand.grid(year = 2015:2024, town = LETTERS[1:8])
  i <- match(g$town, LETTERS)
  g$first <- ifelse(i <= 4, 2019, NA)
  g$on <- as.integer(!is.na(g$first) & g$year >= g$first)
  g$rate <- 30 + i + 0.5 * (g$year - 2015) - 2 * g$on + sin(seq_len(80))
  h <- hwksim(c(log(0.5), 0.5, 2), 100, seed = 1)
  e <- c(1, 1, 2, 4, 4, 5, 3)
  f <- c(2, 3, 3, 5, 6, 6, 4)
  zy <- t(vapply(rep(1:2, each = 20), function(k) rpois(4, 2 * k),
                 integer(4)))
  zx <- cbind(1, 1:4)
  my <- array(rnorm(20 * 3 * 2), c(20, 3, 2))
  calls <- list(
    nnidx = quote(nnidx(p[, 1], p[, 2], win = c(0, 1, 0, 1), nsim = 9)),
    nnidm = quote(nnidx(
      p,
      win = c(0, 1, 0, 1),
      metric = "manhattan",
      nsim = 9
    )),
    hdbang = quote(hdbang(p, v = v)),
    crwalk = quote(crwalk(p, h = 3, nsim = 9, seed = 1)),
    ripenv = quote(ripenv(
      p[, 1],
      p[, 2],
      c(0, 1, 0, 1),
      c(0.1, 0.2),
      nsim = 9
    )),
    mantst = quote(mantst(p, t = v, nsim = 9)),
    fairdf = quote(fairdf(yb, rep(1:2, 50), yb, nboot = 9)),
    calaud = quote(calaud(
      rep(1:10, 10),
      runif(100),
      yb,
      rep(1:2, 50),
      nperm = 9
    )),
    bnft = quote(bnft(bdig(exp(runif(50, 0, 9)))$dig, nsim = 9)),
    skp = quote(skp(
      c(rnorm(5, 1.518, 4e-5), rnorm(5, 1.5195, 4e-5)),
      s2 = 16e-10,
      nsim = 99
    )),
    zgem = quote(zgem(
      zy,
      zx,
      par = cbind(c(0.5, 0), c(1.5, 0)),
      prob = c(0.5, 0.5)
    )),
    mvtem = quote(mvtem(
      my,
      cbind(1, 1:3),
      array(0, c(2, 2, 2)),
      array(diag(2), c(2, 2, 2)),
      prob = c(0.5, 0.5)
    )),
    gstr = quote(gstr(matrix(rpois(100, 2), 10), nsim = 9)),
    sptst = quote(sptst(rpois(20, 9), rpois(20, 9), nboot = 9)),
    gearc = quote(gearc(v, w, nsim = 9)),
    gorg = quote(gorg(exp(v), w, nsim = 9)),
    lmorn = quote(lmorn(v, w, nsim = 9)),
    cntreg = quote(cntreg(rpois(40, 5), rep(c("a", "b"), 20))),
    knxnr = quote(knxnr(p, t = v, ds = 0.3, dt = 0.5, nsim = 9)),
    kscan = quote(kscan(p, t = round(10 * v), r = 0.2, tw = 3, nsim = 9)),
    gmoran = quote(gmoran(v, w, nsim = 9)),
    bvmor = quote(bvmor(v, v, w, nsim = 9)),
    sarml = quote(sarml(v, p[, 2], w)),
    sarmls = quote(sarml(v, p[, 2], w, model = "sac")),
    crossk = quote(crossk(
      p,
      marks = rep(c("a", "b"), 30),
      i = "a",
      j = "b",
      win = c(0, 1, 0, 1),
      r = c(0.05, 0.1),
      sigma = 0.1,
      nsim = 9
    )),
    stpcp = quote(stpcp(c(100, 104, 98, 101, 97, 60, 63, 58), nperm = 9)),
    sdbsim = quote(sdbsim(6, 5, seed = 1)),
    sdbpde = quote(sdbsim(6, 5, mode = "pde", rho = matrix(0.07, 6, 6))),
    rfor = quote(rfor(x, yc, ntree = 5)),
    dmlx = quote(dmlx(
      yc,
      x[, 1],
      x[, -1],
      "plr",
      fold,
      lrn = "forest",
      ntree = 5
    )),
    enpath = quote(enpath(x, yc, nlam = 5)),
    enbin = quote(enpath(x, yb, "binomial", nlam = 5)),
    lgtreg = quote(lgtreg(x, yb)),
    synthc = quote(synthc(
      g$rate,
      g$town,
      g$year,
      ifelse(g$town == "A", 2019, NA)
    )),
    tfrp = quote(tfrp(
      c(
        80,
        50,
        0.5,
        0.05,
        0.8,
        0.9,
        0.1,
        0.2,
        0.4,
        0.6,
        0.05,
        0.1,
        0.5,
        0.7,
        1,
        1,
        0.5
      ),
      ymax = 5,
      nsim = 9
    )),
    stgdid = quote(stgdid(g$rate, g$town, g$year, g$first, reps = 9)),
    kdsurf = quote(kdsurf(p, nx = 5)),
    kdpts = quote(kdsurf(p, h = 0.2, at = "points")),
    kmpart = quote(kmpart(p[, 1], p[, 2], k = 2, seed = 1)),
    webpt = quote(webpt(p[, 1], p[, 2])),
    ntevc = quote(ntevc(e, f)),
    ntchg = quote(ntchg(e, f)),
    ntlvn = quote(ntlvn(e, f)),
    wcbx = quote(wcbx(g$rate, g$town, g$year, g$on)),
    didfe = quote(didfe(g$rate, g$town, g$year, g$on)),
    hwkfit = quote(hwkfit(h$t, 100)),
    etasf = quote(etasf(h$t, 100)),
    mvhwk = quote(mvhwk(h$t, rep(c("a", "b"), length.out = length(h$t)), 100)),
    cconc = quote(cconc(rpois(50, 2), nsim = 9)),
    agemix = quote(agemix(
      cbind(a = 1:3, b = 3:1),
      c(0.5, 0.5),
      1:3,
      panel = cbind(rnorm(9), rnorm(9)),
      nperm = 9
    ))
  )
  for (nm in names(calls)) {
    expect_true(stops(eval(calls[[nm]])), info = nm)
  }
})

test_that("results do not change and no state is left after an interrupt", {
  set.seed(2)
  p <- cbind(x = runif(80), y = runif(80))
  t <- p[, 1] + rnorm(80)
  a <- mantst(p, t = t, nsim = 99, seed = 3)
  expect_true(stops(mantst(p, t = t, nsim = 99, seed = 3)))
  expect_identical(mantst(p, t = t, nsim = 99, seed = 3), a)
})

# The child gets SIGINT from a shell started beside it once it is ready
# (pid file written); it runs in the foreground of its shell, since an
# asynchronous command starts with SIGINT ignored. Without the interrupt
# the call runs for minutes and returns "done".
test_that("Ctrl-C stops a long call in a running R session", {
  skip_on_cran()
  skip_on_os("windows")
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  pid <- file.path(d, "pid")
  out <- file.path(d, "out.rds")
  src <- file.path(d, "run.R")
  writeLines(
    c(
      sprintf(".libPaths(%s)", paste0(deparse(.libPaths()), collapse = "")),
      "library(crimlab)",
      "set.seed(1)",
      "p <- cbind(x = runif(300), y = runif(300))",
      sprintf("writeLines(\"ready\", %s)", deparse(pid)),
      "t0 <- Sys.time()",
      "r <- tryCatch({",
      "  mantst(p, t = p[, 1], nsim = 1e6, seed = 1)",
      "  \"done\"",
      "}, interrupt = function(i) \"stopped\")",
      "s <- as.numeric(Sys.time() - t0, units = \"secs\")",
      "a <- mantst(p, t = p[, 1], nsim = 99, seed = 1)",
      sprintf("saveRDS(list(r = r, s = s, a = a), %s)", deparse(out))
    ),
    src
  )
  rs <- file.path(R.home("bin"), "Rscript")
  sh <- sprintf(
    paste(
      "(while [ ! -s %s ]; do sleep 1; done; sleep 2; kill -INT $$) &",
      "exec %s --vanilla %s"
    ),
    shQuote(pid),
    shQuote(rs),
    shQuote(src)
  )
  system2("sh", c("-c", shQuote(sh)), stdout = FALSE, stderr = FALSE)
  r <- readRDS(out)
  expect_identical(r$r, "stopped")
  expect_lt(r$s, 30)
  set.seed(1)
  p <- cbind(x = runif(300), y = runif(300))
  expect_identical(r$a, mantst(p, t = p[, 1], nsim = 99, seed = 1))
})
