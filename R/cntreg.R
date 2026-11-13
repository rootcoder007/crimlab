#' Count regression: rate ratios between groups
#'
#' Fits a log-linear model of a count on a group, with adjustment
#' columns, an exposure and clusters, by Poisson and negative binomial
#' maximum likelihood, and reports the incidence rate ratio (IRR) of each
#' group against the reference with its interval. The same model serves
#' counts of people, placements, days or deaths in the cells of a
#' published table and counts per person or per area.
#'
#' The Poisson fit takes the steps of [stats::glm()] (iteratively
#' reweighted least squares from `y + 0.1`, stopping when the relative
#' change of the deviance is below 1e-8), and the negative binomial fit
#' (variance `mu + mu^2 / theta`) those of `MASS::glm.nb()`: `theta` by
#' maximum likelihood and the coefficients in turn. A column of the design
#' that the columns before it explain (an adjustment level the group
#' determines, say) is aliased and dropped, as R does: its coefficients are
#' `NA`.
#'
#' With `clus`, each fit also gets cluster-robust (sandwich) standard
#' errors, rows `"poiscl"` and `"nbcl"`: `(X'WX)^-1 (sum_g U_g U_g')
#' (X'WX)^-1` over the clusters' scores `U_g`, times `G / (G - 1)` for `G`
#' clusters when `cadj` (the default of `sandwich::vcovCL()`; `cadj =
#' FALSE` gives the independence GEE variance of 'geepack'). With `gee =
#' TRUE` the Poisson model is also fitted by generalised estimating
#' equations with an exchangeable working correlation (row `"gee"`), as
#' `geepack::geeglm(corstr = "exchangeable")` estimates it.
#'
#' Intervals are Wald intervals on the log scale, `exp(est +- z se)`, with
#' `z` the normal quantile for `level`. `sep` flags a ratio whose group or
#' reference has no counts at all: the estimate then runs off towards 0 or
#' infinity and its interval means nothing.
#'
#' Counts need not be whole numbers (a median duration, say): the
#' estimating equations are those of the Poisson pseudo-likelihood, and the
#' log-likelihood uses `lgamma(y + 1)`, where [stats::glm()] reports an
#' infinite AIC.
#'
#' @param y The counts: finite numbers from 0, not all 0; or a `clpre` list
#'   from [clprep()], whose elements fill the arguments by name (arguments
#'   given here override them).
#' @param grp The group of each count: labels (a factor by its levels, text
#'   or numbers by their sorted values), at least two groups.
#' @param adj `NULL` or adjustment columns: a data frame or list of columns,
#'   a numeric matrix, or one vector, one value per count. Factor, text and
#'   logical columns enter as indicators of each level but the first; give a
#'   number such as a year as a factor for one effect per value.
#' @param expo `NULL` or the exposure of each count (person-days,
#'   population): positive numbers; its logarithm is the offset.
#' @param clus `NULL` or the cluster of each count (institution, region,
#'   year): labels, at least two clusters.
#' @param ref The reference group: a label or a position; the first group
#'   by default.
#' @param gee `TRUE` to fit the exchangeable GEE as well (needs `clus`).
#' @param cadj `TRUE` for the `G / (G - 1)` factor on the cluster-robust
#'   variances.
#' @param level The confidence level of the intervals: one number strictly
#'   between 0 and 1.
#' @return An object of class `cntreg`, a list: `irr`, a data frame with one
#'   row per fit and group (`fit`, `term`, `irr`, `lo`, `hi`, `p`, `sep`);
#'   `coef`, every coefficient of every fit on the log scale (`fit`,
#'   `term`, `est`, `se`, `z`, `p`); `fit`, one row per model (`fit`,
#'   `loglik`, `aic`, `theta`, `setheta`, `alpha`, `iter`, `conv`, `rank`);
#'   `ref`, the reference group; `n`, the number of counts; `ncl`, the
#'   number of clusters (0 without `clus`); and `level`.
#' @references
#'   Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics
#'   with S*, 4th ed. Springer. \doi{10.1007/978-0-387-21706-2}
#'
#'   Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
#'   generalized linear models. *Biometrika*, 73(1), 13-22.
#'   \doi{10.1093/biomet/73.1.13}
#'
#'   Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
#'   cluster-robust inference. *Journal of Human Resources*, 50(2),
#'   317-372. \doi{10.3368/jhr.50.2.317}
#' @examples
#' # placements by group in eight institutions over three years
#' d <- data.frame(
#'   n = c(12, 30, 8, 25, 14, 41, 9, 22, 11, 35, 7, 28),
#'   female = rep(c("no", "yes"), 6),
#'   year = factor(rep(c(2021, 2022, 2023), each = 4)),
#'   inst = rep(c("A", "B", "C", "D"), each = 2, length.out = 12)
#' )
#' cntreg(d$n, d$female, adj = d["year"], clus = d$inst)
#' # rates: a count per person-days of exposure
#' cntreg(c(5, 12), c("a", "b"), expo = c(1000, 1500))
#' @export
cntreg <- function(
  y,
  grp = NULL,
  adj = NULL,
  expo = NULL,
  clus = NULL,
  ref = NULL,
  gee = FALSE,
  cadj = TRUE,
  level = 0.95
) {
  if (inherits(y, "clpre")) {
    return(.prep(cntreg, y, match.call(), parent.frame()))
  }
  .Call(C_clab_r_cntreg, y, grp, adj, expo, clus, ref, gee, cadj, level)
}

#' @export
print.cntreg <- function(x, ...) {
  cat(
    "cntreg: ",
    x[["n"]],
    " counts",
    if (x[["ncl"]] > 0L) {
      paste0(", ", x[["ncl"]], " clusters")
    },
    "; reference ",
    x[["ref"]],
    "\n",
    sep = ""
  )
  r <- x[["irr"]]
  cat("rate ratios (", format(100 * x[["level"]]), "% intervals):\n", sep = "")
  print(r, row.names = FALSE)
  f <- x[["fit"]]
  nb <- f[f[["fit"]] == "nb", ]
  if (nrow(nb) == 1L && !is.na(nb[["theta"]])) {
    cat(
      "negative binomial theta ",
      format(nb[["theta"]]),
      "; AIC Poisson ",
      format(f[["aic"]][1L]),
      ", negative binomial ",
      format(nb[["aic"]]),
      "\n",
      sep = ""
    )
  }
  if (any(r[["sep"]], na.rm = TRUE)) {
    cat(
      "a group or the reference has no counts: those ratios are not ",
      "estimable\n",
      sep = ""
    )
  }
  invisible(x)
}

# The generic f called on a clpre list p from clprep(): p's elements fill
# f's arguments by name, and the arguments of the original call cl (but the
# first, p itself), evaluated in env, override them.
.prep <- function(f, p, cl, env) {
  a <- lapply(as.list(cl)[-(1:2)], eval, envir = env)
  b <- unclass(p)
  attributes(b) <- list(names = names(b))
  b[names(a)] <- a
  do.call(f, b)
}
