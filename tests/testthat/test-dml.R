# Cross-fitted estimators against DoubleML, AIPW, tmle, DRDID and
# rmoriebricklayer, with reference values from dev/ref-dml.R; the
# closed forms and the simulations are in src/t_dml.cpp.
ref <- readRDS(test_path("ref-dml.rds"))
x <- ref$x

# crimlab from the reference's own predictions: estimates and standard
# errors to 1e-10; then from its own least-squares and logistic learners:
# the same predictions to 1e-8, the same estimates to 1e-7
# (DoubleML leaves g1 unfitted for the ATTE, which does not read it: NA
# there, 0 here)
both <- function(r, ...) {
  pr <- r$pred
  na <- !is.finite(pr)
  pr[na] <- 0
  a <- dmlx(..., x = x, pred = pr)
  expect_equal(a$coef, r$coef, tolerance = 1e-10)
  expect_equal(a$se, r$se, tolerance = 1e-10)
  b <- dmlx(..., x = x, lrn = "glm")
  expect_lt(max(abs(b$pred - r$pred)[!na]), 1e-8)
  # DoubleML's logistic learner is glm() at its default deviance tolerance
  expect_equal(b$coef, r$coef, tolerance = 1e-7)
  expect_equal(b$se, r$se, tolerance = 1e-7)
  invisible(b)
}

test_that("dmlx agrees with DoubleML", {
  both(ref$plr, y = ref$yc, d = ref$dc, model = "plr", fold = ref$fold)
  both(ref$plriv, y = ref$yc, d = ref$dc, model = "plr", fold = ref$fold,
       score = "iv")
  both(ref$plr3, y = ref$yc, d = ref$dc, model = "plr", fold = ref$folds3)
  both(ref$pliv, y = ref$yiv, d = ref$dc, z = ref$zc, model = "pliv",
       fold = ref$fold)
  both(ref$irm, y = ref$y, d = ref$d, model = "irm", fold = ref$fold)
  both(ref$atte, y = ref$y, d = ref$d, model = "irm", fold = ref$fold,
       score = "atte")
  both(ref$iivm, y = ref$yz, d = ref$dz, z = ref$z, model = "iivm",
       fold = ref$fold)
  ys <- ifelse(ref$s == 1, ref$y, NA)
  both(ref$ssm, y = ys, d = ref$d, s = ref$s, model = "ssm",
       fold = ref$fold)
  # one-way clusters, DoubleML's own cluster folds
  c1 <- dmlx(ref$yc, ref$dc, x, "plr", ref$clu$fold, clu = ref$clu$cl,
             pred = ref$clu$pred)
  expect_equal(c1$coef, ref$clu$coef, tolerance = 1e-10)
  expect_equal(c1$se, ref$clu$se, tolerance = 1e-10)
  # two-way clusters, DoubleML's folds of each kind
  tw <- dmlx(ref$yc, ref$dc, x, "plr", ref$two$fold, clu = ref$two$cl,
             pred = ref$two$pred)
  expect_equal(tw$coef, ref$two$coef, tolerance = 1e-10)
  expect_equal(tw$se, ref$two$se, tolerance = 1e-10)
  tg <- dmlx(ref$yc, ref$dc, x, "plr", ref$two$fold, clu = ref$two$cl,
             lrn = "glm")
  expect_lt(max(abs(tg$pred[, , 1] - ref$two$pred)), 1e-8)
  expect_equal(tg$coef, ref$two$coef, tolerance = 1e-8)
})

test_that("dmlx agrees with AIPW, tmle and DRDID on the same nuisances", {
  pr <- ref$irm$pred[, , 1]
  a <- dmlx(ref$y, ref$d, x, "aipw", ref$fold, pred = pr)
  expect_equal(a$coef, ref$aipw$coef, tolerance = 1e-10)
  expect_equal(a$se, ref$aipw$se, tolerance = 1e-10)
  # tmle solves its fluctuation with glm(), to glm's deviance tolerance
  tm <- dmlx(ref$y, ref$d, x, "tmle", ref$fold, pred = pr)
  expect_equal(tm$coef, ref$tmle$coef, tolerance = 1e-7)
  expect_equal(tm$se, ref$tmle$se, tolerance = 1e-7)
  di <- dmlx(ref$did$dy, ref$d, x, "did", ref$fold, pred = ref$did$pred,
             trim = 1e-12)
  expect_equal(di$coef, ref$did$coef, tolerance = 1e-10)
})

test_that("dmlx agrees with rmoriebricklayer where the definitions meet", {
  pp <- ref$plr$pred[, , 1]
  expect_equal(dmlx(ref$yc, ref$dc, x, "plr", ref$fold, pred = pp)$coef,
               ref$blplr, tolerance = 1e-10)
  pr <- ref$irm$pred[, , 1]
  a <- dmlx(ref$y, ref$d, x, "irm", ref$fold, pred = pr, trim = 0.01)
  expect_equal(a$coef, ref$blATE$coef, tolerance = 1e-10)
  expect_equal(a$se, ref$blATE$se, tolerance = 1e-10)
  k <- ref$keep
  b <- dmlx(ref$y[k], ref$d[k], x[k, ], "irm", ref$sfold, score = "atte",
            pred = pr[k, ], trim = 0.01)
  expect_equal(b$coef, ref$blATT$coef, tolerance = 1e-10)
  expect_equal(b$se, ref$blATT$se, tolerance = 1e-10)
})

test_that("dmlx learns with the elastic net and repeats the folds", {
  r <- dmlx(ref$yc, ref$dc, x, "plr", ref$folds3, nlam = 20, nfin = 3)
  expect_identical(dim(r$pred), c(600L, 2L, 3L))
  expect_identical(dimnames(r$pred)[[2]], c("l", "m"))
  expect_named(r$reps, c("coef", "se"))
  expect_lt(abs(r$coef - 0.5), 4 * r$se)
  expect_identical(r$coef, median(r$reps$coef))
  # the same seed, the same inner folds and the same result
  expect_identical(dmlx(ref$yc, ref$dc, x, "plr", ref$folds3, nlam = 20,
                        nfin = 3), r)
  e <- dmlx(ref$y, ref$d, x, "didcs", ref$fold, t = ref$s, lrn = "ridge",
            lrnd = "enet", nlam = 10, nfin = 2)
  expect_identical(dimnames(e$pred)[[2]],
                   c("g00", "g01", "g10", "g11", "m"))
  expect_true(is.finite(e$coef))
  expect_gt(e$se, 0)
  ex <- dmlx(ref$y, ref$d, x, "did", ref$fold, score = "exp", lrn = "glm")
  expect_identical(dimnames(ex$pred)[[2]], c("g0", "g1"))
  ec <- dmlx(ref$y, ref$d, x, "didcs", ref$fold, t = ref$s, score = "exp",
             lrn = "glm")
  expect_identical(dimnames(ec$pred)[[2]], c("g00", "g01", "g10", "g11"))
  hf <- dmlx(ref$y, ref$d, x, "irm", ref$fold, lrn = "forest", ntree = 50,
             hon = TRUE)
  expect_lt(abs(hf$coef - 1), 4 * hf$se)
  fr <- dmlx(ref$y, ref$d, x, "irm", ref$fold, lrn = "forest", ntree = 50,
             mtry = 9)
  expect_lt(abs(fr$coef - 1), 4 * fr$se)
})

test_that("dmlx stops on bad input, naming the argument", {
  y <- ref$yc
  d <- ref$dc
  f <- ref$fold
  errs("dmlx", list(
    list("`model` must be one of", y, d, x, "xyz", f),
    list("`x` must be a numeric matrix", y, d, y, "plr", f),
    list("four rows", y[1:3], d[1:3], x[1:3, ], "plr", f[1:3]),
    list("`y` must have length", y[-1], d, x, "plr", f),
    list("`d` must be 0 or 1", y, d, x, "irm", f),
    list("`z` must be finite", y, d, x, "pliv", f),
    list("`z` must be 0 or 1", ref$yz, ref$dz, x, "iivm", f, z = d),
    list("`s` must be 0 or 1", y, ref$d, x, "ssm", f, s = d),
    list("`y` must be numbers", "a", ref$d, x, "ssm", f, s = ref$s),
    list("`y` must be finite where", NA * y, ref$d, x, "ssm", f, s = ref$s),
    list("`t` must be 0 or 1", y, ref$d, x, "didcs", f, t = d),
    list("`score` must be one of", y, d, x, "plr", f, score = "ate"),
    list("`score` must be one of", y, ref$d, x, "irm", f, score = "po"),
    list("`score` must be one of", y, ref$d, x, "did", f, score = "ate"),
    list("`score` must be NULL", y, ref$d, x, "aipw", f, score = "ate"),
    list("one row per row", y, d, x, "plr", f[-1]),
    list("`fold` must be finite", y, d, x, "plr", as.character(f)),
    list("same folds in every column", y, d, x, "plr",
         cbind(f, rep(1:2, 300))),
    list("`clu` is for models", y, ref$d, x, "irm", f, clu = f),
    list("constant within each cluster", y, d, x, "plr", f,
         clu = rep(1:2, 300)),
    list("two-way `clu` needs", y, d, x, "plr", f, clu = ref$two$cl),
    list("constant within each cluster", y, d, x, "plr",
         cbind(ref$two$fold[, 1], rep(1:3, 200)), clu = ref$two$cl),
    list("`alpha`", y, d, x, "plr", f, alpha = 2),
    list("`nlam`", y, d, x, "plr", f, nlam = 0),
    list("`nfin`", y, d, x, "plr", f, nfin = 1),
    list("`tol`", y, d, x, "plr", f, tol = 0),
    list("`scr`", y, d, x, "plr", f, scr = NA),
    list("`hon`", y, d, x, "plr", f, hon = "yes"),
    list("`seed`", y, d, x, "plr", f, seed = -1),
    list("`lrn` must be one of", y, d, x, "plr", f, lrn = "tree"),
    list("`ntree`", y, d, x, "plr", f, lrn = "forest", ntree = 0),
    list("`mtry`", y, d, x, "plr", f, lrn = "forest", mtry = 0),
    list("`lrnd` must be one of", y, d, x, "plr", f, lrnd = 1),
    list("`trim` must be at least 0", y, ref$d, x, "irm", f, trim = -1),
    list("`trim` must be below 0.5", y, ref$d, x, "irm", f, trim = 0.5),
    list("`pred` must be an array", y, d, x, "plr", f,
         pred = matrix(0, 600, 3)),
    list("`pred` must be an array", y, d, x, "plr", f, pred = 1:3),
    list("`pred` must be an array", y, d, x, "plr", f,
         pred = array(0, c(600, 2, 2))),
    list("`pred` must be finite", y, d, x, "plr", f,
         pred = matrix(NA_real_, 600, 2)),
    list("no estimate", y, d, x, "plr", f, pred = cbind(0, d))
  ))
})
