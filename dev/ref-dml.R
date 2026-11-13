# Reference values for tests/testthat/test-dml.R:
#   Rscript dev/ref-dml.R      (from the package root)
# Needs DoubleML, mlr3, mlr3learners, data.table, AIPW, tmle, DRDID and
# rmoriebricklayer. On fixed simulated data and fixed folds it records,
# for each model: what DoubleML (R) estimates with least-squares and
# logistic-regression learners, the out-of-fold predictions it used, and
# its estimate and standard error; the AIPW package's and the tmle
# package's results from given nuisance predictions; DRDID's doubly robust
# difference-in-differences estimate from given predictions; and
# rmoriebricklayer's partialling-out and AIPW formulas on given
# predictions. Writes tests/testthat/ref-dml.rds.
suppressPackageStartupMessages({
  for (p in c("DoubleML", "mlr3", "mlr3learners", "data.table", "AIPW",
              "tmle", "DRDID", "rmoriebricklayer")) {
    stopifnot(requireNamespace(p, quietly = TRUE))
  }
})
lgr::get_logger("mlr3")$set_threshold("warn")
set.seed(20261112)
n <- 600
p <- 4
x <- matrix(rnorm(n * p), n, dimnames = list(NULL, paste0("x", 1:p)))
e <- plogis(0.8 * x[, 1] - 0.4 * x[, 2])
d <- rbinom(n, 1, e)
dc <- 0.7 * x[, 1] - 0.3 * x[, 3] + rnorm(n)
zc <- x[, 2] + rnorm(n)
z <- rbinom(n, 1, 0.5)
dz <- ifelse(z == 1, rbinom(n, 1, plogis(1 + x[, 1])), rbinom(n, 1, 0.1))
s <- rbinom(n, 1, plogis(0.5 + x[, 2] + 0.5 * d))
y <- d + x[, 1] + 0.5 * x[, 2] + rnorm(n)
yc <- 0.5 * dc + x[, 1] - x[, 3] + rnorm(n)
yiv <- 0.5 * (dc + zc) + x[, 3] + rnorm(n)
yz <- 2 * dz + x[, 1] + rnorm(n)
fold <- sample(rep(1:5, length.out = n))
folds3 <- cbind(fold, sample(rep(1:5, length.out = n)),
                sample(rep(1:5, length.out = n)))
smpl <- function(f) {
  lapply(seq_len(ncol(f)), function(r) {
    list(train_ids = lapply(1:5, function(k) which(f[, r] != k)),
         test_ids = lapply(1:5, function(k) which(f[, r] == k)))
  })
}
regr <- mlr3::lrn("regr.lm")
clsf <- mlr3::lrn("classif.log_reg")
dt <- function(...) data.table::data.table(x, ...)
# one DoubleML fit: estimate, standard error, predictions (n x k x reps)
run <- function(obj, nm, f = cbind(fold)) {
  obj$set_sample_splitting(smpl(f))
  obj$fit(store_predictions = TRUE)
  pr <- simplify2array(lapply(nm, function(k) {
    matrix(obj$predictions[[k]], n, ncol(f))
  }))
  list(coef = unname(obj$coef), se = unname(obj$se),
       pred = aperm(array(pr, c(n, ncol(f), length(nm))), c(1, 3, 2)))
}
dmd <- function(yy, dd, zz = NULL, ss = NULL, ...) {
  df <- dt(y = yy, d = dd)
  if (!is.null(zz)) df$z <- zz
  if (!is.null(ss)) df$s <- ss
  DoubleML::DoubleMLData$new(df, y_col = "y", d_cols = "d",
                             z_cols = if (is.null(zz)) NULL else "z",
                             s_col = if (is.null(ss)) NULL else "s", ...)
}
ref <- list(x = x, y = y, d = d, dc = dc, yc = yc, zc = zc, yiv = yiv,
            z = z, dz = dz, yz = yz, s = s, fold = fold, folds3 = folds3)
ref$plr <- run(DoubleML::DoubleMLPLR$new(dmd(yc, dc), regr$clone(),
                                         regr$clone()), c("ml_l", "ml_m"))
ref$plriv <- run(DoubleML::DoubleMLPLR$new(dmd(yc, dc), regr$clone(),
                                           regr$clone(), regr$clone(),
                                           score = "IV-type"),
                 c("ml_l", "ml_m", "ml_g"))
ref$plr3 <- run(DoubleML::DoubleMLPLR$new(dmd(yc, dc), regr$clone(),
                                          regr$clone(), n_rep = 3),
                c("ml_l", "ml_m"), folds3)
ref$pliv <- run(DoubleML::DoubleMLPLIV$new(dmd(yiv, dc, zc), regr$clone(),
                                           regr$clone(), regr$clone()),
                c("ml_l", "ml_m", "ml_r"))
ref$irm <- run(DoubleML::DoubleMLIRM$new(dmd(y, d), regr$clone(),
                                         clsf$clone()),
               c("ml_g0", "ml_g1", "ml_m"))
ref$atte <- run(DoubleML::DoubleMLIRM$new(dmd(y, d), regr$clone(),
                                          clsf$clone(), score = "ATTE"),
                c("ml_g0", "ml_g1", "ml_m"))
ref$iivm <- run(DoubleML::DoubleMLIIVM$new(dmd(yz, dz, z), regr$clone(),
                                           clsf$clone(), clsf$clone()),
                c("ml_g0", "ml_g1", "ml_m", "ml_r0", "ml_r1"))
ys <- ifelse(s == 1, y, 0)
ref$ssm <- run(DoubleML::DoubleMLSSM$new(dmd(ys, d, ss = s), regr$clone(),
                                         clsf$clone(), clsf$clone()),
               c("ml_pi", "ml_m", "ml_g_d0", "ml_g_d1"))
# clusters: DoubleML draws the cluster folds itself; they are recorded
cl <- rep(1:60, each = 10)
cdat <- DoubleML::DoubleMLClusterData$new(dt(y = yc, d = dc, cl = cl),
                                          y_col = "y", d_cols = "d",
                                          cluster_cols = "cl",
                                          x_cols = colnames(x))
cobj <- DoubleML::DoubleMLPLR$new(cdat, regr$clone(), regr$clone())
cobj$fit(store_predictions = TRUE)
cf <- integer(n)
for (k in 1:5) cf[cobj$smpls[[1]]$test_ids[[k]]] <- k
ref$clu <- list(coef = unname(cobj$coef), se = unname(cobj$se), cl = cl,
                fold = cf,
                pred = cbind(cobj$predictions$ml_l, cobj$predictions$ml_m))
# AIPW: its summary from the IRM predictions, filled in by hand
pr <- ref$irm$pred[, , 1]
a <- AIPW::AIPW$new(Y = y, A = d, W = x, Q.SL.library = "SL.glm",
                    g.SL.library = "SL.glm", k_split = 2, verbose = FALSE)
a$obs_est$mu0 <- pr[, 1]
a$obs_est$mu1 <- pr[, 2]
a$obs_est$mu <- ifelse(d == 1, pr[, 2], pr[, 1])
a$obs_est$raw_p_score <- pr[, 3]
a$summary(g.bound = 0.025)
ref$aipw <- list(coef = unname(a$estimates$RD[["Estimate"]]),
                 se = unname(a$estimates$RD[["SE"]]))
# tmle: the fluctuation of given predictions
tm <- tmle::tmle(Y = y, A = d, W = as.data.frame(x),
                 Q = cbind(Q0W = pr[, 1], Q1W = pr[, 2]), g1W = pr[, 3],
                 family = "gaussian")
ref$tmle <- list(coef = tm$estimates$ATE$psi,
                 se = sqrt(tm$estimates$ATE$var.psi))
# DRDID: the doubly robust DiD on a panel with in-sample nuisances (the
# estimate coincides with the normalised observational score; its
# standard error accounts for the in-sample fits, so it differs)
dy <- 0.5 * x[, 1] + d + rnorm(n)
y0 <- rnorm(n)
ps <- fitted(glm(d ~ x, family = binomial()))
od <- drop(cbind(1, x) %*% coef(lm(dy ~ x, subset = d == 0)))
ref$did <- list(dy = dy, pred = cbind(od, ps),
                coef = DRDID::drdid_panel(y0 + dy, y0, d, cbind(1, x),
                                          boot = FALSE)$ATT)
# rmoriebricklayer: the partialling-out estimate and the AIPW formulas on
# the DoubleML predictions
bl <- asNamespace("rmoriebricklayer")
pp <- ref$plr$pred[, , 1]
ref$blplr <- bl$.rmbl_dml_theta(yc - pp[, 1], dc - pp[, 2], yc, dc,
                                seq_len(n))$theta
# rmoriebricklayer's ATT takes the treated share of the whole sample,
# DoubleML's ATTE that of each test fold: they coincide when every fold
# has the same share, so the rows kept make both counts multiples of 5
# and the folds deal each group round robin
kt <- which(d == 1)[seq_len(5 * (sum(d) %/% 5))]
kc <- which(d == 0)[seq_len(5 * (sum(1 - d) %/% 5))]
keep <- sort(c(kt, kc))
sfold <- integer(n)
sfold[kt] <- rep(1:5, length.out = length(kt))
sfold[kc] <- rep(1:5, length.out = length(kc))
ref$keep <- keep
ref$sfold <- sfold[keep]
for (est in c("ATE", "ATT")) {
  k <- if (est == "ATT") keep else seq_len(n)
  sc <- bl$.rmbl_aipw_score(est, pr[k, 2], pr[k, 1], d[k], y[k],
                            pmin(pmax(pr[k, 3], 0.01), 0.99),
                            rep(1, length(k)))
  ref[[paste0("bl", est)]] <- list(
    coef = sc$th,
    se = bl$.rmbl_aipw_se(NULL, NULL, est, sc, d[k], rep(1, length(k)),
                          length(k))$se
  )
}
ref$made <- list(R = R.version.string,
                 pkgs = vapply(c("DoubleML", "AIPW", "tmle", "DRDID",
                                 "rmoriebricklayer"),
                               function(p) as.character(packageVersion(p)), ""))
saveRDS(ref, "tests/testthat/ref-dml.rds", version = 2)
