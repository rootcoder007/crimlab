#' Alert complexity and placement volatility, pooled or by fiscal year
#'
#' The analysis of Ruhela (2026) on any release of Ontario's restrictive
#' confinement records: from the person-year table ([otpyr()]) it gives,
#' for the years pooled or for each fiscal year, the share of person-years
#' with complex alert profiles (`cx = 1`), mean placement volatility by
#' group, the crude Poisson rate ratio (the ratio of the two means,
#' infinite when every count without complex alerts is 0), the one-to-one
#' nearest-neighbour propensity matching without replacement of each
#' complex person-year to a control on age group, gender and fiscal year
#' ([cbatt()], in 'MatchIt''s default order), the count model on the
#' matched frame with a random intercept per regional cluster
#' ([cntmix()]), and the bounds on the average effect on the matched
#' frame with no assumption and under monotone treatment response and
#' selection ([mtrbnd()]).
#'
#' A single fiscal year cannot identify a year term: `by = "pool"` with one
#' year and `yterm = TRUE` is refused; `by = "year"` drops the term in each
#' year and says so. A sample with no complex person-year, or none without,
#' is refused, naming it.
#'
#' @param x A person-year table from [otpyr()], or what [otpyr()] takes (a
#'   data frame of the release or a table key).
#' @param by `"pool"` for the years together, `"year"` for each year.
#' @param yterm `TRUE` to adjust for the fiscal year in the pooled
#'   propensity and count models.
#' @param fam The family of [cntmix()].
#' @param ... Passed to [otpyr()] when `x` is not a person-year table.
#' @return An object of class `cxvol`, a list: `res`, one row per sample
#'   (`sample`, `n`, `n1`, `share`, `vm`, `vm1`, `vm0`, `crude`, `nm`,
#'   `vm1m`, `vm0m`, `irr`, `lo`, `hi`, `est`, `se`, `p`, `ncl`, `fam`,
#'   `mlo`, `mhi`, `bhi`: the Manski bounds and the upper bound under both
#'   monotonicity assumptions, whose lower bound is 0; `ymax`, the largest
#'   count on the matched frame, the range used by the bounds; `note`, why
#'   the rate ratio has no estimate when it has none, such as a group whose
#'   counts on the matched frame are all 0, else `""`); `fits`, the
#'   [cntmix()] fits; `by`; `yterm`.
#' @inherit plcvol references
#' @examples
#' set.seed(4)
#' n <- 600
#' p <- data.frame(year = sample(2025:2026, n, TRUE),
#'                 gender = sample(c("Male", "Female"), n, TRUE, c(9, 1)),
#'                 age = sample(c("18 to 24", "25 to 49", "50+"), n, TRUE),
#'                 rc = sample(c("A", "B", "A + B", "C"), n, TRUE),
#'                 cx = rbinom(n, 1, 0.2))
#' p$vm <- rpois(n, exp(-2 + 0.3 * p$cx + (p$rc == "A + B")))
#' cxvol(p, fam = "poisson")
#' @export
cxvol <- function(
  x,
  by = c("pool", "year"),
  yterm = TRUE,
  fam = c("nb", "poisson"),
  ...
) {
  fn <- "cxvol"
  by <- .ckmarg("by", fn)
  fam <- .ckmarg("fam", fn)
  yterm <- .ckflag(yterm, "yterm", fn)
  p <- .cxpy(x, fn, ...)
  s <- .cxsplit(p, by, yterm, fn)
  out <- lapply(names(s), function(nm) {
    .cxone(s[[nm]]$d, s[[nm]]$yt, fam, nm, fn)
  })
  structure(
    list(
      res = do.call(rbind, lapply(out, `[[`, "row")),
      fits = stats::setNames(lapply(out, `[[`, "fit"), names(s)),
      by = by,
      yterm = yterm
    ),
    class = "cxvol"
  )
}

#' @export
print.cxvol <- function(x, ...) {
  cat(
    "cxvol: alert complexity (cx) and placement volatility (vm), ",
    if (x$by == "pool") "years pooled" else "each fiscal year (no year term)",
    if (x$by == "pool" && !x$yterm) ", no year term",
    "\n",
    sep = ""
  )
  r <- x$res
  print(
    data.frame(
      sample = r$sample,
      n = r$n,
      cx = r$share,
      vm1 = r$vm1,
      vm0 = r$vm0,
      crude = r$crude,
      matched = r$nm,
      irr = r$irr,
      lo = r$lo,
      hi = r$hi,
      p = signif(r$p, 3),
      bound = r$bhi
    ),
    row.names = FALSE,
    digits = 4
  )
  cat(
    "irr: rate ratio of the matched count model (",
    r$fam[1L],
    "); bound: upper end of [0, bound] under monotone response and\n",
    "selection (below 0, the data contradict the two assumptions)\n",
    sep = ""
  )
  nt <- nzchar(r$note)
  if (any(nt)) {
    cat(paste0("note (", r$sample[nt], "): ", r$note[nt], "\n"), sep = "")
  }
  invisible(x)
}

#' Double machine learning of the effect of alert complexity
#'
#' The interactive regression model of Chernozhukov et al. (2018) for the
#' effect of complex alert profiles (`cx`) on placement volatility (`vm`)
#' on the person-year table, by [dmlx()]: random forests (or another
#' learner of [dmlx()]) for the outcome regressions and the propensity,
#' `nfold`-fold cross-fitting, propensities kept within `trim` of 0 and 1,
#' covariates gender, age group, regional cluster and (pooled, with
#' `yterm`) fiscal year. Both the average effect (ATE) and the effect on
#' the treated (ATTE) are given, from one set of cross-fitted nuisances.
#' The same refusals as [cxvol()] apply.
#'
#' @inheritParams cxvol
#' @param lrn The learner of [dmlx()], for both nuisances.
#' @param ntree The trees of a forest.
#' @param nfold The folds of the cross-fitting: a whole number from 2.
#' @param trim Propensities are kept within `trim` of 0 and 1.
#' @param seed The seed of the folds and of the learners.
#' @return An object of class `cxdml`: a data frame with one row per sample
#'   and estimand (`sample`, `n`, `estimand`, `est`, `se`, `lo`, `hi`).
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E.,
#'   Hansen, C., Newey, W. and Robins, J. (2018). Double/debiased machine
#'   learning for treatment and structural parameters. *The Econometrics
#'   Journal*, 21(1), C1-C68. \doi{10.1111/ectj.12097}
#' @examples
#' set.seed(5)
#' n <- 400
#' p <- data.frame(year = sample(2025:2026, n, TRUE),
#'                 gender = sample(c("Male", "Female"), n, TRUE),
#'                 age = sample(c("18 to 24", "25 to 49"), n, TRUE),
#'                 rc = sample(c("A", "B"), n, TRUE), cx = rbinom(n, 1, 0.3))
#' p$vm <- rpois(n, exp(-1.5 + 0.4 * p$cx))
#' cxdml(p, lrn = "glm")
#' @export
cxdml <- function(
  x,
  by = c("pool", "year"),
  yterm = TRUE,
  lrn = "forest",
  ntree = 500,
  nfold = 3,
  trim = 0.01,
  seed = 1,
  ...
) {
  fn <- "cxdml"
  by <- .ckmarg("by", fn)
  yterm <- .ckflag(yterm, "yterm", fn)
  nfold <- .cknum(nfold, "nfold", fn, min = 2, max = 1e6, whole = TRUE)
  seed <- .ckseed(seed, "seed", fn)
  s <- .cxsplit(.cxpy(x, fn, ...), by, yterm, fn)
  out <- do.call(
    rbind,
    lapply(names(s), function(nm) {
      .cxdone(s[[nm]]$d, s[[nm]]$yt, nm, lrn, ntree, nfold, trim, seed, fn)
    })
  )
  structure(out, class = c("cxdml", "data.frame"))
}

#' @export
print.cxdml <- function(x, ...) {
  cat("cxdml: interactive model, effect of cx on vm (additive scale)\n")
  print(as.data.frame(unclass(x)), row.names = FALSE, digits = 4)
  invisible(x)
}

#' Placement volatility events by fiscal year
#'
#' The year-over-year table of placement volatility (Ruhela 2026,
#' Appendix E, Table E15) from the person-year table: for each fiscal year
#' the regional changes (`vm` summed), placements, person-years, the rates
#' of changes per 1,000 placements and per 1,000 person-years with exact
#' Poisson intervals, each year's share of all changes, and the percent
#' changes from the previous year of the count and of both rates. The rate
#' intervals are `qgamma(a/2, n)` to `qgamma(1 - a/2, n + 1)` over the
#' exposure, so a count of 0 has lower limit 0.
#'
#' The two rates can disagree in direction: placements per person-year
#' move too, and the placement denominator absorbs their growth.
#'
#' @inheritParams cxvol
#' @param level The confidence level of the gamma intervals of both rates:
#'   one number strictly between 0 and 1.
#' @return A data frame of class `cxyoy`, one row per fiscal year: `year`,
#'   `vm`, `np`, `py`, `rnp`, `rnplo`, `rnphi` (per 1,000 placements),
#'   `rpy`, `rpylo`, `rpyhi` (per 1,000 person-years), `share`, `dvm`,
#'   `drnp`, `drpy` (percent changes from the previous year).
#' @inherit plcvol references
#' @examples
#' p <- data.frame(year = rep(2024:2026, c(50, 60, 70)), vm = rpois(180, 0.2),
#'                 np = rpois(180, 20) + 1)
#' cxyoy(p)
#' @export
cxyoy <- function(x, level = 0.95, ...) {
  fn <- "cxyoy"
  level <- .xpopen(level, "level", fn)
  p <- if (is.data.frame(x) && all(c("year", "vm", "np") %in% names(x))) {
    x
  } else {
    .clas(fn, otpyr(x, ...))
  }
  p <- .cyck(p, fn)
  y <- factor(p$year, sort(unique(p$year)))
  vm <- as.vector(tapply(p$vm, y, sum))
  np <- as.vector(tapply(p$np, y, sum))
  py <- as.vector(table(y))
  a <- (1 - level) / 2
  lo <- ifelse(vm == 0, 0, stats::qgamma(a, vm))
  hi <- stats::qgamma(1 - a, vm + 1)
  pc <- function(v) c(NA, 100 * (v[-1L] / v[-length(v)] - 1))
  out <- data.frame(
    year = levels(y),
    vm = vm,
    np = np,
    py = py,
    rnp = 1000 * vm / np,
    rnplo = 1000 * lo / np,
    rnphi = 1000 * hi / np,
    rpy = 1000 * vm / py,
    rpylo = 1000 * lo / py,
    rpyhi = 1000 * hi / py,
    share = vm / sum(vm),
    dvm = pc(vm),
    drnp = pc(vm / np),
    drpy = pc(vm / py)
  )
  structure(out, class = c("cxyoy", "data.frame"))
}

#' @export
print.cxyoy <- function(x, ...) {
  cat("cxyoy: regional changes per fiscal year\n")
  print(as.data.frame(unclass(x)), row.names = FALSE, digits = 4)
  invisible(x)
}

#' cxyoy()'s table checked: a person-year or more, every year given, the
#' counts of changes and of placements whole numbers from 0, and each
#' year with a placement (the denominator of its rate).
#' @noRd
.cyck <- function(p, fn) {
  if (!nrow(p) || anyNA(p$year)) {
    .ckstop(fn, "`x` must have a person-year or more, each with its `year`")
  }
  p$vm <- .cknum(p$vm, "vm", fn, scalar = FALSE, min = 0, whole = TRUE)
  p$np <- .cknum(p$np, "np", fn, scalar = FALSE, min = 0, whole = TRUE)
  np <- tapply(p$np, p$year, sum)
  if (any(np == 0)) {
    .ckstop(fn, "year %s has no placement: no rate per placement",
            names(np)[np == 0][1L])
  }
  p
}

#' The person-year table, built when x is not one; its columns checked.
#' @noRd
.cxpy <- function(x, fn, ...) {
  p <- if (
    inherits(x, "otpyr") ||
      (is.data.frame(x) &&
        all(
          c("year", "vm", "cx", "gender", "age", "rc") %in% names(x)
        ))
  ) {
    x
  } else {
    .clas(fn, otpyr(x, ...))
  }
  need <- c("year", "vm", "cx", "gender", "age", "rc")
  if (!all(need %in% names(p)) || anyNA(p[need])) {
    .ckstop(
      fn,
      "`x` must be a person-year table with %s, none missing",
      toString(need)
    )
  }
  if (is.logical(p$cx)) {
    p$cx <- as.numeric(p$cx)
  }
  if (!is.numeric(p$cx) || !all(p$cx %in% 0:1)) {
    .ckstop(fn, "`cx` must be 0 or 1 (or TRUE or FALSE) in every person-year")
  }
  p$vm <- .cknum(p$vm, "vm", fn, scalar = FALSE, min = 0, whole = TRUE)
  p
}

#' The samples: all years, or one per year; whether each takes a year term.
#' @noRd
.cxsplit <- function(p, by, yterm, fn) {
  yrs <- sort(unique(p$year))
  if (by == "pool") {
    if (yterm && length(yrs) < 2L) {
      .ckstop(
        fn,
        paste0(
          "one fiscal year (%s) cannot identify a year term: ",
          "set yterm = FALSE, or by = \"year\""
        ),
        yrs
      )
    }
    nm <- if (length(yrs) > 1L) {
      paste(range(yrs), collapse = "-")
    } else {
      as.character(yrs)
    }
    return(stats::setNames(
      list(list(d = p, yt = yterm && length(yrs) > 1L)),
      nm
    ))
  }
  stats::setNames(
    lapply(yrs, function(y) list(d = p[p$year == y, ], yt = FALSE)),
    as.character(yrs)
  )
}

#' One sample: descriptives, crude ratio, matching, count model, bounds.
#' @noRd
.cxone <- function(d, yt, fam, nm, fn) {
  # cx is 0 or 1 in every row (.cxpy()), never missing
  if (!any(d$cx == 1) || !any(d$cx == 0)) {
    .ckstop(
      fn,
      "sample %s has no person-year with cx = %d: no contrast can be estimated",
      nm,
      if (any(d$cx == 1)) 0L else 1L
    )
  }
  ds <- .cxdes(d[c("age", "gender", if (yt) "year")], nm, fn)
  cv <- ds$cv
  X <- ds$X
  m <- .clas(fn, cbatt(d$vm, d$cx, X, est = "psm", tgt = "att", reps = 0))
  rows <- c(m$match$focal, m$match$match)
  md <- d[rows, ]
  fit <- .clas(
    fn,
    cntmix(
      md$vm,
      data.frame(cx = md$cx, cv[rows, , drop = FALSE]),
      md$rc,
      fam = fam
    )
  )
  list(fit = fit, row = .cxrow(d, md, fit, nm, fn))
}

#' The row of one sample: descriptives, crude ratio, model and bounds.
#' @noRd
.cxrow <- function(d, md, fit, nm, fn) {
  z <- factor(md$cx, 0:1)
  b0 <- .clas(fn, mtrbnd(md$vm, z, 0, max(md$vm), "none"))$bnd
  b2 <- .clas(fn, mtrbnd(md$vm, z, 0, max(md$vm), "both"))$bnd
  t1 <- fit$coef[fit$coef$term == "cx", ]
  data.frame(
    sample = nm,
    n = nrow(d),
    n1 = sum(d$cx),
    share = mean(d$cx),
    vm = mean(d$vm),
    vm1 = mean(d$vm[d$cx == 1]),
    vm0 = mean(d$vm[d$cx == 0]),
    crude = mean(d$vm[d$cx == 1]) / mean(d$vm[d$cx == 0]),
    nm = nrow(md),
    vm1m = mean(md$vm[md$cx == 1]),
    vm0m = mean(md$vm[md$cx == 0]),
    t1[c("irr", "lo", "hi", "est", "se", "p")],
    ncl = fit$ngrp,
    fam = fit$fam,
    mlo = b0$atelo[2L],
    mhi = b0$atehi[2L],
    bhi = b2$atehi[2L],
    ymax = max(md$vm),
    note = paste(grep("(^|, )cx(, |: )", fit$notes, value = TRUE),
                 collapse = "; "),
    row.names = NULL
  )
}

#' The covariates as factors (cv) and as indicator columns (X). A covariate
#' with one value in the sample (one gender in a year, say) tells nothing
#' apart and is left out of both; model.matrix() would stop on it.
#' @noRd
.cxdes <- function(cv, nm, fn) {
  cv[] <- lapply(cv, function(v) factor(as.character(v)))
  cv <- cv[vapply(cv, nlevels, 1L) > 1L]
  if (!length(cv)) {
    .ckstop(
      fn,
      "sample %s has one value of every covariate: nothing to adjust for",
      nm
    )
  }
  X <- .dsgn(cv, fn)[, -1L, drop = FALSE]
  list(cv = cv, X = X[, apply(X, 2L, stats::var) > 0, drop = FALSE])
}

#' One sample's ATE and ATTE, the nuisances learned once.
#' @noRd
.cxdone <- function(d, yt, nm, lrn, ntree, nfold, trim, seed, fn) {
  if (sum(d$cx == 1) < nfold || sum(d$cx == 0) < nfold) {
    .ckstop(
      fn,
      "sample %s has too few person-years in a group for %d folds",
      nm,
      nfold
    )
  }
  X <- .cxdes(d[c("gender", "age", "rc", if (yt) "year")], nm, fn)$X
  fold <- .clseed(seed, sample(rep_len(seq_len(nfold), nrow(d))))
  a <- .clas(
    fn,
    dmlx(
      d$vm,
      d$cx,
      X,
      "irm",
      fold,
      score = "ate",
      lrn = lrn,
      lrnd = lrn,
      ntree = ntree,
      trim = trim,
      seed = seed
    )
  )
  b <- .clas(
    fn,
    dmlx(d$vm, d$cx, X, "irm", fold, score = "atte", pred = a$pred, trim = trim)
  )
  z <- stats::qnorm(0.975)
  data.frame(
    sample = nm,
    n = nrow(d),
    estimand = c("ATE", "ATTE"),
    est = c(a$coef, b$coef),
    se = c(a$se, b$se),
    lo = c(a$coef, b$coef) - z * c(a$se, b$se),
    hi = c(a$coef, b$coef) + z * c(a$se, b$se)
  )
}
