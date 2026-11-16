#' Generalisability weights for a trial
#'
#' The weights that make the units of a randomised trial stand for a
#' target population, and how far the trial is from it. The data hold the
#' trial's units (`s = 1`) and a sample of the target (`s = 0`) with the
#' same covariates. The sampling score `ps`, the chance of being in the
#' trial given the covariates, is a logistic regression of `s` on `x`,
#' trimmed to `[trim, 1 - trim]`.
#'
#' With `target = "generalize"` the target is the population the trial
#' was drawn from, all the rows, and a trial unit's weight is `1 / ps`;
#' with `"transport"` it is the units outside the trial, and the weight
#' is the odds `(1 - ps) / ps` (Stuart et al. 2011; Westreich et al.
#' 2017; Tipton and Hartman 2023, eq. 3.11-3.12). Rows outside the trial
#' have weight 0.
#'
#' The diagnostics: `smd`, the difference of the mean logits of the
#' scores between the trial and the target over their standard deviation
#' in the target (Stuart et al. 2011); `dmean`, the difference of the mean
#' scores; `gindex`, Tipton's (2014) generalisability index, the sum over
#' `bins` equal-width bins of the scores of the square root of the trial's
#' share times the target's, 1 when they are spread alike and 0 when they
#' never meet (above 0.9 is very high, below 0.5 low); `ess`, Kish's
#' effective number of trial units under the weights. `bal` gives each
#' covariate's mean in the target and in the trial, unweighted and
#' weighted, with the standardised differences from the target (over its
#' standard deviation in the target).
#'
#' @param data A data frame, one row per unit of the trial or the target.
#' @param s The name of the column marking the trial's units: 1 in the
#'   trial, 0 outside it (or `TRUE`, `FALSE`), both present.
#' @param x The names of the covariate columns (numbers) of the sampling
#'   model; for `print()`, an object of class `gensw`.
#' @param target `"generalize"` or `"transport"`.
#' @param trim The bound of the scores, from 0 to below 0.5.
#' @param bins The bins of Tipton's index: a whole number from 1 to
#'   10000.
#' @return An object of class `gensw`, a list: `ps` and `w`, one per row;
#'   `coef`, the sampling model's coefficients; `diag`, a named vector
#'   (`ntrial`, `ntarget`, `smd`, `dmean`, `gindex`, `ess`); `bal`, a data
#'   frame (`var`, `target`, `trial`, `wtrial`, `smd`, `wsmd`); and
#'   `target`.
#' @references
#'   Stuart, E. A., Cole, S. R., Bradshaw, C. P. and Leaf, P. J. (2011).
#'   The use of propensity scores to assess the generalizability of
#'   results from randomized trials. *Journal of the Royal Statistical
#'   Society, Series A*, 174(2), 369-386.
#'   \doi{10.1111/j.1467-985X.2010.00673.x}
#'
#'   Tipton, E. (2014). How generalizable is your experiment? An index for
#'   comparing experimental samples and populations. *Journal of
#'   Educational and Behavioral Statistics*, 39(6), 478-501.
#'   \doi{10.3102/1076998614558486}
#'
#'   Westreich, D., Edwards, J. K., Lesko, C. R., Stuart, E. and Cole, S.
#'   R. (2017). Transportability of trial results using inverse odds of
#'   sampling weights. *American Journal of Epidemiology*, 186(8),
#'   1010-1014. \doi{10.1093/aje/kwx164}
#'
#'   Tipton, E. and Hartman, E. (2023). Generalizability and
#'   transportability. In Zubizarreta, J. R., Stuart, E. A., Small, D. S.
#'   and Rosenbaum, P. R. (eds), *Handbook of Matching and Weighting
#'   Adjustments for Causal Inference*, ch. 3. Chapman and Hall/CRC.
#'   \doi{10.1201/9781003102670}
#' @examples
#' set.seed(8)
#' n <- 600
#' d <- data.frame(age = rnorm(n, 35, 10), prior = rpois(n, 1.5))
#' d$s <- rbinom(n, 1, plogis(-1 + 0.04 * (d$age - 35) + 0.3 * d$prior))
#' g <- gensw(d, "s", c("age", "prior"), target = "transport")
#' g
#' g$bal
#' @export
gensw <- function(data, s, x, target = c("generalize", "transport"),
                  trim = 0, bins = 10) {
  f <- "gensw"
  dd <- .axdat(data, f)
  .axcols(dd, s, "s", f, n = 1L)
  .axcols(dd, x, "x", f)
  if (s %in% x) .ckstop(f, "`s` and `x` must name different columns")
  target <- .ckmarg("target", f)
  sv <- .axmat(dd, s, f)[, 1L]
  xm <- .axmat(dd, x, f)
  if (!all(sv %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", s)
  if (all(sv == 1) || all(sv == 0)) {
    .ckstop(f, "column `%s` must mark units in and out of the trial", s)
  }
  trim <- .cknum(trim, "trim", f, min = 0)
  if (trim >= 0.5) .ckstop(f, "`trim` must be below 0.5")
  bins <- .cknum(bins, "bins", f, min = 1, max = 1e4, whole = TRUE)
  k <- .Call(C_clab_r_gensw, xm, as.integer(sv), target == "transport",
             trim, as.integer(bins))
  st <- k$st
  b <- matrix(k$bal, length(x))
  cf <- k$coef
  names(cf) <- c("(Intercept)", x)
  structure(
    list(ps = k$ps, w = k$w, coef = cf,
         diag = c(ntrial = st[1L], ntarget = st[2L], smd = st[3L],
                  dmean = st[4L], gindex = st[5L], ess = st[6L]),
         bal = data.frame(var = x, target = b[, 1L], trial = b[, 2L],
                          wtrial = b[, 3L], smd = b[, 4L], wsmd = b[, 5L],
                          stringsAsFactors = FALSE),
         target = target),
    class = "gensw"
  )
}

#' @rdname gensw
#' @param ... Not used by `print()`.
#' @export
print.gensw <- function(x, ...) {
  g <- x[["diag"]]
  cat("gensw: ", g[["ntrial"]], " trial units to ", x[["target"]], " to ",
      g[["ntarget"]], "; logit SMD ", format(g[["smd"]], digits = 3),
      ", Tipton's index ", format(g[["gindex"]], digits = 3),
      ", effective size ", format(g[["ess"]], digits = 4), "\n", sep = "")
  print(x[["bal"]], row.names = FALSE, digits = 4)
  invisible(x)
}
