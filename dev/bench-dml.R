# Speed of crimlab's learners and cross-fitted estimators against the
# reference packages, on fixed seeded data:
#   Rscript dev/bench-dml.R [n...]       e.g. Rscript dev/bench-dml.R 5000
# BENCH_N and BENCH_P (comma lists) override the sizes; the default grid is
# n = 5000, 50000 and p = 20, 200. For each task the estimates (or the
# errors, for the forests) are compared first; then each side is timed
# (the median of BENCH_REPS runs, default 3; one run when a single run
# exceeds 60 s). Writes inst/bench/dml.csv and inst/bench/dml.md; nothing
# in them is typed by hand.
suppressPackageStartupMessages(library(crimlab))
has <- function(p) requireNamespace(p, quietly = TRUE)
num <- function(v, d) {
  if (nzchar(v)) as.numeric(strsplit(v, ",")[[1]]) else d
}
args <- commandArgs(TRUE)
ns <- if (length(args)) as.numeric(args) else
  num(Sys.getenv("BENCH_N"), c(5000, 50000))
ps <- num(Sys.getenv("BENCH_P"), c(20, 200))
reps <- num(Sys.getenv("BENCH_REPS"), 3)
if (has("lgr")) lgr::get_logger("mlr3")$set_threshold("warn")
# AIPW and tmle look their learners up on the search path
if (has("SuperLearner")) suppressPackageStartupMessages(library(SuperLearner))

# the 1-minute load average (Linux /proc, or sysctl on macOS); a timing
# whose load rose by more than 1.5 while it ran is taken again (twice at
# most), and every row records the load before and after its runs
load1 <- function() {
  if (file.exists("/proc/loadavg")) {
    return(as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1]))
  }
  v <- tryCatch(system("sysctl -n vm.loadavg", intern = TRUE),
                error = function(e) NA_character_)
  as.numeric(strsplit(gsub("[{}]", "", trimws(v)), " +")[[1]][2])
}
lstart <- load1()
cpu <- if (file.exists("/proc/cpuinfo")) {
  m <- grep("^model name|^Model", readLines("/proc/cpuinfo"), value = TRUE)
  if (length(m)) trimws(sub(".*:", "", m[1])) else NA_character_
} else {
  tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE),
           error = function(e) NA_character_)
}
# one thread everywhere: crimlab is single-threaded, and so are the
# counterparts here (ranger and grf told so, data.table and BLAS by
# environment)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")
if (has("data.table")) data.table::setDTthreads(1)
cat("crimlab from", find.package("crimlab"), "\n")
for (p in c("crimlab", "glmnet", "ranger", "grf", "DoubleML", "mlr3",
            "mlr3learners", "ddml", "AIPW", "tmle", "SuperLearner")) {
  if (has(p)) cat(p, as.character(packageVersion(p)), "\n")
}

lnow <- c(NA_real_, NA_real_)
tm <- function(f) {
  for (try in 1:3) {
    l0 <- load1()
    t1 <- system.time(r <- f())[["elapsed"]]
    t <- if (t1 > 60 || reps < 2) t1 else c(t1, vapply(seq_len(reps - 1), function(i) {
      system.time(f())[["elapsed"]]
    }, 0))
    l1 <- load1()
    if (!isTRUE(l1 - l0 > 1.5)) break
  }
  lnow <<- c(l0, l1)
  list(t = stats::median(t), r = r)
}

# the data: sparse linear nuisances, as the dossier's lasso benchmark
sim <- function(n, p) {
  set.seed(n + p)
  x <- matrix(rnorm(n * p), n, dimnames = list(NULL, paste0("x", 1:p)))
  b <- c(1, -0.5, 0.5, 0.25, rep(0, p - 4))
  dc <- drop(x %*% (b / 2)) + rnorm(n)
  d <- rbinom(n, 1, plogis(drop(x %*% (b / 2))))
  z <- rbinom(n, 1, 0.5)
  dz <- z * rbinom(n, 1, plogis(drop(x %*% (b / 2))))
  list(x = x, dc = dc, d = d, z = z, dz = dz, zc = x[, 5] + rnorm(n),
       yc = 0.5 * dc + drop(x %*% b) + rnorm(n),
       y = d + drop(x %*% b) + rnorm(n),
       yz = dz + drop(x %*% b) + rnorm(n),
       fold = sample(rep(1:5, length.out = n)))
}

rows <- list()
add <- function(task, n, p, ours, ref, refpkg, agree, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    task = task, n = n, p = p, secs = round(ours, 3),
    ref = refpkg, refsecs = if (is.na(ref)) NA else round(ref, 3),
    speedup = if (is.na(ref)) NA else round(ref / ours, 2),
    agree = agree, note = note, lbefore = lnow[1], lafter = lnow[2],
    stringsAsFactors = FALSE)
  cat(sprintf("%-10s n=%-6d p=%-4d crimlab %8.2fs  %-10s %8.2fs  agree=%s %s\n",
              task, n, p, ours, refpkg, ref, agree, note))
}
# estimates agree when they differ by less than half a standard error:
# the inner folds that choose the penalty differ between packages
near <- function(a, b, se) isTRUE(abs(a - b) < 0.5 * se)

for (n in ns) for (p in ps) {
  s <- sim(n, p)
  x <- s$x
  sm <- list(list(train_ids = lapply(1:5, function(k) which(s$fold != k)),
                  test_ids = lapply(1:5, function(k) which(s$fold == k))))
  # the elastic net: one cross-validated lasso path, the same folds
  a <- tm(function() enpath(x, s$yc, fold = s$fold))
  if (has("glmnet")) {
    g <- tm(function() glmnet::cv.glmnet(x, s$yc, foldid = s$fold))
    # the same penalties and the same coefficients on the full path, to
    # the default threshold (1e-7) of both; either may stop a penalty
    # before the other
    k <- seq_len(min(length(a$r$lambda), length(g$r$lambda)))
    gb <- as.matrix(g$r$glmnet.fit$beta)[, k]
    ok <- isTRUE(all.equal(a$r$lambda[k], g$r$lambda[k], tolerance = 1e-10)) &&
      max(abs(a$r$beta[, k] - gb)) < 1e-3 * max(abs(gb))
    add("enpath", n, p, a$t, g$t, "glmnet", ok,
        sprintf("cv.glmnet; %d and %d penalties fitted",
                length(a$r$lambda), length(g$r$lambda)))
  }
  # the forests: test error on fresh rows from the same model
  xt <- sim(2000, p)
  ft <- tm(function() rfor(x, s$yc, xt$x, ntree = 100))
  mse <- function(pr) mean((pr - xt$yc)^2)
  if (has("ranger")) {
    r <- tm(function() {
      ranger::ranger(x = x, y = s$yc, num.trees = 100, num.threads = 1,
                     seed = 1)
    })
    e1 <- mse(ft$r$pred)
    e2 <- mse(predict(r$r, xt$x, num.threads = 1)$predictions)
    add("rfor", n, p, ft$t, r$t, "ranger", e1 < 1.1 * e2,
        sprintf("test MSE %.3f vs %.3f; ranger single-threaded", e1, e2))
  }
  if (has("grf")) {
    r <- tm(function() {
      grf::regression_forest(x, s$yc, num.trees = 100, num.threads = 1,
                             seed = 1)
    })
    e2 <- mse(predict(r$r, xt$x)$predictions)
    add("rfor", n, p, ft$t, r$t, "grf", mse(ft$r$pred) < 1.1 * e2,
        sprintf("test MSE %.3f vs %.3f; grf honest, single-threaded",
                mse(ft$r$pred), e2))
  }
  # the cross-fitted estimators with lasso nuisances
  ours <- list(
    plr = tm(function() dmlx(s$yc, s$dc, x, "plr", s$fold)),
    pliv = tm(function() dmlx(s$yc, s$dc, x, "pliv", s$fold, z = s$zc)),
    irm = tm(function() dmlx(s$y, s$d, x, "irm", s$fold)),
    iivm = tm(function() dmlx(s$yz, s$dz, x, "iivm", s$fold, z = s$z)),
    aipw = tm(function() dmlx(s$y, s$d, x, "aipw", s$fold)),
    tmle = tm(function() dmlx(s$y, s$d, x, "tmle", s$fold))
  )
  if (has("DoubleML") && has("mlr3learners")) {
    lr <- mlr3::lrn("regr.cv_glmnet", s = "lambda.min")
    lc <- mlr3::lrn("classif.cv_glmnet", s = "lambda.min",
                    predict_type = "prob")
    dd <- function(y, d, z = NULL) {
      df <- data.table::data.table(x, y = y, d = d)
      if (!is.null(z)) df$z <- z
      DoubleML::DoubleMLData$new(df, y_col = "y", d_cols = "d",
                                 z_cols = if (is.null(z)) NULL else "z")
    }
    one <- function(nm, obj) {
      obj$set_sample_splitting(sm)
      r <- tm(function() obj$fit())
      o <- ours[[nm]]$r
      add(nm, n, p, ours[[nm]]$t, r$t, "DoubleML",
          near(o$coef, obj$coef, o$se))
    }
    one("plr", DoubleML::DoubleMLPLR$new(dd(s$yc, s$dc), lr$clone(),
                                        lr$clone()))
    one("pliv", DoubleML::DoubleMLPLIV$new(dd(s$yc, s$dc, s$zc), lr$clone(),
                                          lr$clone(), lr$clone()))
    one("irm", DoubleML::DoubleMLIRM$new(dd(s$y, s$d), lr$clone(),
                                        lc$clone()))
    # nobody with z = 0 is treated: DoubleML must be told (crimlab sees a
    # constant outcome and predicts it)
    one("iivm", DoubleML::DoubleMLIIVM$new(dd(s$yz, s$dz, s$z), lr$clone(),
                                          lc$clone(), lc$clone(),
                                          subgroups = list(
                                            always_takers = FALSE,
                                            never_takers = TRUE)))
  }
  if (has("ddml")) {
    sub <- lapply(1:5, function(k) which(s$fold == k))
    l <- list(what = ddml::mdl_glmnet, args = list(cv = TRUE))
    r <- tm(function() {
      ddml::ddml_plm(s$yc, s$dc, x, learners = l, sample_folds = 5,
                     subsamples = sub, silent = TRUE)
    })
    add("plr", n, p, ours$plr$t, r$t, "ddml",
        near(ours$plr$r$coef, r$r$coef[1], ours$plr$r$se))
  }
  if (has("AIPW") && has("SuperLearner") && has("glmnet")) {
    r <- tm(function() {
      AIPW::AIPW$new(Y = s$y, A = s$d, W = x, Q.SL.library = "SL.glmnet",
                     g.SL.library = "SL.glmnet", k_split = 5,
                     verbose = FALSE)$fit()$summary()
    })
    add("aipw", n, p, ours$aipw$t, r$t, "AIPW",
        near(ours$aipw$r$coef, r$r$estimates$RD[["Estimate"]],
             ours$aipw$r$se), "AIPW cross-fits SuperLearner's glmnet")
  }
  if (has("tmle") && has("SuperLearner") && has("glmnet")) {
    r <- tm(function() {
      tmle::tmle(s$y, s$d, as.data.frame(x), Q.SL.library = "SL.glmnet",
                 g.SL.library = "SL.glmnet", family = "gaussian")
    })
    add("tmle", n, p, ours$tmle$t, r$t, "tmle",
        near(ours$tmle$r$coef, r$r$estimates$ATE$psi, ours$tmle$r$se),
        "tmle cross-validates SuperLearner's glmnet, 10 folds")
  }
}
tab <- do.call(rbind, rows)
tab$faster <- !is.na(tab$speedup) & tab$speedup > 1
lend <- load1()
tab$machine <- Sys.info()[["nodename"]]
tab$cpu <- cpu
tab$r <- R.version.string
tab$lstart <- lstart
tab$lend <- lend
tab$date <- format(Sys.time(), "%Y-%m-%d %H:%M")

dir.create("inst/bench", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/bench/dml.csv", row.names = FALSE)
md <- c(
  "# crimlab speed against reference packages",
  "",
  sprintf("Generated by dev/bench-dml.R on %s: %s (%s, %s), %s; load average %.2f at the start and %.2f at the end; seconds of elapsed time, the median of %d runs.",
          tab$date[1], tab$machine[1], cpu, Sys.info()[["machine"]],
          R.version.string, lstart, lend, reps),
  "",
  "The maintainer's own packages (rmorie, rmoriebricklayer) are not in this table; every reference here is a third-party package. One thread on both sides (ranger and grf with num.threads = 1; OMP, BLAS and data.table threads set to 1). A machine that runs other sessions idles near a load of 2; a timing whose load rose by more than 1.5 while it ran was taken again, and each row's load before and after is in the CSV.",
  "",
  "| task | n | p | crimlab (s) | reference | reference (s) | speed-up | estimates agree | note |",
  "|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %d | %.3f | %s | %s | %s | %s | %s |", tab$task, tab$n,
          tab$p, tab$secs, tab$ref,
          format(tab$refsecs), format(tab$speedup), tab$agree, tab$note)
)
writeLines(md, "inst/bench/dml.md")
