#' Falsification with negative-control outcomes
#'
#' Estimates the "effect" of the treatment on outcomes it cannot have
#' affected (negative controls: an outcome measured before treatment, or
#' one with no plausible path from it), by the same regression as the
#' outcome of interest. An association with a negative control points to
#' confounding or another bias that the adjustment did not remove
#' (Lipsitch, Tchetgen Tchetgen and Cohen 2010; Shi, Miao and Tchetgen
#' Tchetgen 2020).
#'
#' Each outcome is regressed on an intercept, the covariates `x` and the
#' treatment `d` by least squares; `est` is the coefficient of `d` with a
#' classical (`se = "iid"`) or HC1 (`"hc1"`) standard error, Student's t
#' on the residual degrees of freedom and its interval at `level`. With
#' `nperm > 0` each outcome also has the permutation p-value of Kennedy
#' (1995): the residuals of the outcome on the covariates are permuted
#' against those of `d`, one permutation of the rows for every outcome,
#' `pperm = (1 + k) / (1 + nperm)`. The controls' p-values (the
#' permutation ones when drawn) are adjusted for their number by `adj`
#' (the methods of `p.adjust()`); `nflag` counts the controls whose
#' adjusted p-value is at most `alpha`.
#'
#' @inheritParams plntfx
#' @param y `NULL` or the name of the outcome of interest (numbers).
#' @param nc The names of the negative-control outcome columns (numbers).
#' @param x `NULL` or the names of the covariate columns (numbers); for
#'   `print()`, an object of class `falsctl`.
#' @param nperm The permutations: a whole number from 0 (none) to
#'   1000000.
#' @param adj The adjustment of the controls' p-values: `"holm"`,
#'   `"bonferroni"`, `"hochberg"`, `"BH"`, `"BY"` or `"none"`.
#' @param alpha The level at which a control is flagged: one number
#'   strictly between 0 and 1.
#' @param seed `NULL` (drawn from R's generator, so [set.seed()] fixes it)
#'   or a whole number from 0 to 2147483647; not used without
#'   permutations.
#' @return An object of class `falsctl`, a list: `tab`, a data frame with
#'   one row per outcome (`outcome`, `role`, `primary` or `control`, `est`,
#'   `se`, `stat`, `p`, `lo`, `hi`, `pperm`, `padj`; `NA` where not
#'   computed); `nflag`; `ncontrol`; and `se`, `adj`, `alpha`, `level`,
#'   `nperm`, `seed` (`NA` without permutations).
#' @references
#'   Lipsitch, M., Tchetgen Tchetgen, E. and Cohen, T. (2010). Negative
#'   controls: a tool for detecting confounding and bias in observational
#'   studies. *Epidemiology*, 21(3), 383-388.
#'   \doi{10.1097/EDE.0b013e3181d61eeb}
#'
#'   Shi, X., Miao, W. and Tchetgen Tchetgen, E. (2020). A selective review
#'   of negative control methods in epidemiology. *Current Epidemiology
#'   Reports*, 7(4), 190-202. \doi{10.1007/s40471-020-00243-4}
#'
#'   Kennedy, P. E. (1995). Randomization tests in econometrics. *Journal
#'   of Business and Economic Statistics*, 13(1), 85-94.
#'   \doi{10.1080/07350015.1995.10524581}
#' @examples
#' set.seed(5)
#' n <- 300
#' d <- data.frame(age = runif(n, 18, 60), u = rnorm(n))
#' d$t <- rbinom(n, 1, plogis(0.8 * d$u))
#' d$arrest <- 0.4 * d$t + 0.5 * d$u + rnorm(n)
#' d$prior <- 0.5 * d$u + rnorm(n)     # before treatment: a control
#' d$height <- rnorm(n)                # no path from treatment
#' falsctl(d, "arrest", "t", c("prior", "height"), x = "age")
#' @export
falsctl <- function(
  data,
  y = NULL,
  d,
  nc,
  x = NULL,
  se = c("hc1", "iid"),
  nperm = 0,
  adj = c("holm", "bonferroni", "hochberg", "BH", "BY", "none"),
  alpha = 0.05,
  level = 0.95,
  seed = NULL
) {
  f <- "falsctl"
  dd <- .axdat(data, f)
  .axcols(dd, y, "y", f, n = 1L, null = TRUE)
  .axcols(dd, d, "d", f, n = 1L)
  .axcols(dd, nc, "nc", f)
  .axcols(dd, x, "x", f, null = TRUE)
  if (any(nc %in% c(y, d, x)) || d %in% c(y, x) ||
        (!is.null(y) && y %in% x)) {
    .ckstop(f, "`y`, `d`, `nc` and `x` must name different columns")
  }
  se <- .ckmarg("se", f)
  adj <- .ckmarg("adj", f)
  ym <- .axmat(dd, c(y, nc), f)
  dv <- .axmat(dd, d, f)[, 1L]
  xm <- .axmat(dd, x, f)
  if (nrow(dd) <= length(x) + 2L) {
    .ckstop(f, "`data` needs more rows than covariates plus two")
  }
  if (all(dv == dv[1L])) .ckstop(f, "column `%s` must vary", d)
  nperm <- .cknum(nperm, "nperm", f, min = 0, max = 1e6, whole = TRUE)
  alpha <- .xpopen(alpha, "alpha", f)
  level <- .xpopen(level, "level", f)
  seed <- if (nperm > 0) .xpseed(seed, f) else NA_real_
  ctl <- c(if (!is.null(y)) 0L, rep(1L, length(nc)))
  o <- .Call(C_clab_r_falsctl, xm, dv, ym, ctl, se == "hc1", level,
             match(adj, c("none", "bonferroni", "holm", "hochberg", "BH",
                          "BY")) - 1L, as.integer(nperm), seed)
  m <- ncol(ym)
  col <- function(j) o[(j - 1L) * m + seq_len(m)]
  tab <- data.frame(outcome = c(y, nc),
                    role = ifelse(ctl == 1L, "control", "primary"),
                    est = col(1L), se = col(2L), stat = col(3L), p = col(4L),
                    lo = col(5L), hi = col(6L), pperm = col(7L),
                    padj = col(8L), stringsAsFactors = FALSE)
  tab$pperm[is.nan(tab$pperm)] <- NA_real_
  tab$padj[is.nan(tab$padj)] <- NA_real_
  structure(
    list(tab = tab, nflag = sum(tab$padj <= alpha, na.rm = TRUE),
         ncontrol = length(nc), se = se, adj = adj, alpha = alpha,
         level = level, nperm = nperm, seed = seed),
    class = "falsctl"
  )
}

#' @rdname falsctl
#' @inheritParams cemtch
#' @export
print.falsctl <- function(x, ...) {
  cat("falsctl: ", x[["nflag"]], " of ", x[["ncontrol"]],
      " negative control(s) associated with the treatment at alpha ",
      format(x[["alpha"]]), " (", x[["adj"]], ")\n", sep = "")
  print(x[["tab"]], row.names = FALSE, digits = 4)
  invisible(x)
}
