#' Simulation of the alert complexity estimators
#'
#' Draws person-year tables from the model of Ruhela (2026) and reports
#' the bias, spread and interval coverage of the estimators of
#' [cxvol()]. Each table has `n` person-years in `ncl` regional clusters
#' with normal intercepts (standard deviation `sdcl`), an age group (three
#' levels) and a gender, a treatment drawn from a logistic model in which
#' the covariates and the cluster raise the chance of a complex profile,
#' and a negative binomial count (dispersion `theta`) whose log mean is
#' `b0 + beta * cx` plus the covariates and the cluster. With `conf` the
#' person-year also has an unobserved trait that raises both the chance of
#' treatment and the count (by `conf` on the logit and on the log mean):
#' ignorability given the observed covariates then fails.
#'
#' The estimators: `crude`, the Poisson rate ratio of `cx` alone; `mixed`,
#' [cntmix()] on the full table with the covariates (log rate ratio,
#' Wald interval); `matched`, [cntmix()] on the frame matched as in
#' [cxvol()]. The target is `beta`, the log rate ratio given the covariates
#' and the cluster, which the model holds constant ([proof()] names the
#' log-link identity). Results are `bias` (mean estimate less `beta`),
#' `sd`, `rmse`, and `cover`, the share of 95% intervals holding `beta`
#' (`NA` for `crude`, which reports no interval). Each replicate is a
#' fresh table from `seed`; a fit that fails is left out and counted.
#'
#' @param reps The replicates: a whole number from 1 to 10000000.
#' @param n The person-years of each table: a whole number from 50 to
#'   10000000.
#' @param beta The log rate ratio of `cx`, from -20 to 20.
#' @param ncl,sdcl The clusters and the standard deviation of their
#'   intercepts.
#' @param theta The negative binomial dispersion.
#' @param b0 The intercept of the log mean, from -20 to 20.
#' @param conf 0, or the effect of the unobserved trait, from -20 to 20.
#' @param seed The seed of the whole run, set once before the first table is
#'   drawn: one whole number. The session's own random stream is left as it
#'   was.
#' @return An object of class `cxsim`, a list: `res`, one row per estimator
#'   (`est`, `bias`, `sd`, `rmse`, `cover`, `ok` the replicates fitted);
#'   `draws`, every replicate's estimates; and the settings.
#' @examples
#' cxsim(reps = 1, n = 100, ncl = 5, b0 = 0, seed = 1)
#' @export
cxsim <- function(
  reps = 200,
  n = 2000,
  beta = log(1.3),
  ncl = 24,
  sdcl = 1,
  theta = 2,
  b0 = -2.5,
  conf = 0,
  seed = 1
) {
  fn <- "cxsim"
  reps <- .cknum(reps, "reps", fn, min = 1, max = 1e7, whole = TRUE)
  n <- .cknum(n, "n", fn, min = 50, max = 1e7, whole = TRUE)
  ncl <- .cknum(ncl, "ncl", fn, min = 2, max = n / 10, whole = TRUE)
  # sizes a simulation can mean: the log rates stay far from overflow
  beta <- .cknum(beta, "beta", fn, min = -20, max = 20)
  sdcl <- .cknum(sdcl, "sdcl", fn, min = 0, max = 10)
  b0 <- .cknum(b0, "b0", fn, min = -20, max = 20)
  conf <- .cknum(conf, "conf", fn, min = -20, max = 20)
  theta <- .cknum(theta, "theta", fn, positive = TRUE)
  dr <- .clseed(.ckseed(seed, "seed", fn), do.call(
    rbind,
    lapply(seq_len(reps), function(r) {
      .csone(.csdraw(n, beta, ncl, sdcl, theta, b0, conf))
    })
  ))
  res <- do.call(
    rbind,
    lapply(c("crude", "mixed", "matched"), function(e) {
      v <- dr[dr$est == e & is.finite(dr$b), ]
      data.frame(
        est = e,
        bias = mean(v$b) - beta,
        sd = stats::sd(v$b),
        rmse = sqrt(mean((v$b - beta)^2)),
        cover = if (e == "crude") NA else mean(v$lo <= beta & beta <= v$hi),
        ok = nrow(v)
      )
    })
  )
  structure(
    list(res = res, draws = dr, beta = beta, n = n, reps = reps, conf = conf),
    class = "cxsim"
  )
}

#' @export
print.cxsim <- function(x, ...) {
  cat(
    "cxsim: ",
    x$reps,
    " tables of ",
    x$n,
    " person-years, log rate ratio ",
    format(x$beta, digits = 4),
    if (x$conf != 0) sprintf(", unobserved confounder %g", x$conf),
    "\n",
    sep = ""
  )
  print(x$res, row.names = FALSE, digits = 3)
  invisible(x)
}

#' One person-year table from the model.
#' @noRd
.csdraw <- function(n, beta, ncl, sdcl, theta, b0, conf) {
  cl <- sample.int(ncl, n, TRUE)
  u <- stats::rnorm(ncl, 0, sdcl)[cl]
  age <- sample(c("18 to 24", "25 to 49", "50+"), n, TRUE, c(0.15, 0.74, 0.11))
  male <- stats::rbinom(n, 1, 0.9)
  z <- stats::rnorm(n)
  lp <- -2 + 0.3 * (age == "18 to 24") + 0.2 * male + 0.5 * u + conf * z
  cx <- stats::rbinom(n, 1, stats::plogis(lp))
  mu <- exp(b0 + beta * cx + 0.2 * (age == "50+") - 0.1 * male + u + conf * z)
  data.frame(
    vm = stats::rnbinom(n, size = theta, mu = mu),
    cx = cx,
    age = age,
    gender = ifelse(male == 1, "Male", "Female"),
    rc = cl,
    year = 1
  )
}

#' The three estimates of one table.
#' @noRd
.csone <- function(d) {
  # a fit that fails is counted; its warnings (glm.fit()'s on a draw it
  # cannot fit) are not the caller's
  ft <- function(f) tryCatch(suppressWarnings(f), error = function(e) NULL)
  cr <- ft(stats::glm.fit(cbind(1, d$cx), d$vm, family = stats::poisson()))
  mx <- ft(cntmix(d$vm, d[c("cx", "age", "gender")], d$rc, fam = "nb"))
  mt <- ft(cxvol(d, yterm = FALSE, fam = "nb")$fits[[1L]])
  get1 <- function(f) {
    if (is.null(f)) {
      return(c(NA, NA, NA))
    }
    k <- f$coef[f$coef$term == "cx", ]
    c(k$est, log(k$lo), log(k$hi))
  }
  m <- rbind(
    c(if (is.null(cr)) NA else cr$coefficients[2L], NA, NA),
    get1(mx),
    get1(mt)
  )
  data.frame(
    est = c("crude", "mixed", "matched"),
    b = m[, 1L],
    lo = m[, 2L],
    hi = m[, 3L]
  )
}
