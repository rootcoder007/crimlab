#' Coarsened exact matching
#'
#' The matching of Iacus, King and Porro (2012): each covariate is
#' coarsened into bins, the rows are grouped by their bins on every
#' covariate at once, and only the groups (strata) holding both treated
#' and control rows are kept, so that a later regression compares like
#' with like without extrapolating.
#'
#' A numeric covariate is cut into equal-width bins over its range,
#' `nclass.Sturges()` of them unless `cut` gives the number of bins or the
#' breaks (bins closed on the left, the last closed on both sides); a
#' character, factor or logical covariate is matched exactly on its
#' values. A kept stratum `s` with `m_T^s` treated and `m_C^s` control
#' rows gives its treated rows weight 1 and its controls weight
#' `(m_T^s / m_C^s) (m_C / m_T)`, `m_T` and `m_C` the matched totals:
#' with these weights the controls stand for the treated rows'
#' distribution, and a weighted regression of the outcome on the
#' treatment (with `w` as weights, on the rows of `keep`) estimates the
#' effect on the treated (Iacus, King and Porro 2012, sec. 3; Blackwell et
#' al. 2009, p. 529-531).
#'
#' `l1` is the multivariate imbalance `L1` of Iacus, King and Porro
#' (2011): half the sum over the cells of a coarsening of the absolute
#' differences between the treated and the control shares of each cell, 0
#' for identical histograms and 1 for disjoint ones, before matching and
#' after (with the weights). On the matching coarsening it is 0 after
#' matching by construction, so it is measured on the coarsening of
#' `lcut`, by default Scott's rule (`nclass.scott()` equal-width bins) for
#' each numeric covariate.
#'
#' Sampson and Winter (2018) define children with a mean blood lead level
#' of 6 micrograms per decilitre or more as treated, match them on
#' categorical traits of the child, the caregiver and the neighbourhood
#' (Blackwell et al. 2009), report `L1` falling "from .71 to nearly 0"
#' (p. 283), and then fit ordinary least squares on the matched sample;
#' they assess mediators (impulsivity, anxiety) by adding them to that
#' regression (p. 289).
#'
#' Talayero et al. (2023) review seventeen studies of measured lead (blood,
#' bone or dentine) and later delinquent or criminal behaviour, rating
#' each for confounding, and find "limited data at the individual level"
#' (abstract).
#'
#' @param data A data frame, one row per unit.
#' @param d The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`).
#' @param x The names of the covariate columns to match on.
#' @param cut `NULL` or a named list: for a numeric covariate, the number
#'   of bins (a whole number from 1) or the breaks (increasing, covering
#'   its range).
#' @param lcut The same for the coarsening on which `l1` is measured.
#' @return An object of class `cemtch`, a list: `keep` (logical, the rows
#'   matched), `strata` (the stratum of each row, from 1; `NA` when not
#'   matched), `w` (the weights, 0 when not matched), `tab` (rows by group,
#'   `all`, `matched` and `unmatched`), `l1` (`before`, `after`), `nstrata`
#'   (the strata kept and in all) and `cut` (the breaks used).
#' @references
#'   Iacus, S. M., King, G. and Porro, G. (2012). Causal inference without
#'   balance checking: coarsened exact matching. *Political Analysis*, 20,
#'   1-24. \doi{10.1093/pan/mpr013}
#'
#'   Iacus, S. M., King, G. and Porro, G. (2011). Multivariate matching
#'   methods that are monotonic imbalance bounding. *Journal of the
#'   American Statistical Association*, 106, 345-361.
#'   \doi{10.1198/jasa.2011.tm09599}
#'
#'   Blackwell, M., Iacus, S., King, G. and Porro, G. (2009). cem: coarsened
#'   exact matching in Stata. *The Stata Journal*, 9, 524-546.
#'
#'   Sampson, R. J. and Winter, A. S. (2018). Poisoned development:
#'   assessing childhood lead exposure as a cause of crime in a birth
#'   cohort followed through adolescence. *Criminology*, 56(2), 269-301.
#'   \doi{10.1111/1745-9125.12171}
#'
#'   Talayero, M. J., Robbins, C. R., Smith, E. R. and Santos-Burgoa, C.
#'   (2023). The association between lead exposure and crime: a systematic
#'   review. *PLOS Global Public Health*, 3(8), e0002177.
#'   \doi{10.1371/journal.pgph.0002177}
#' @examples
#' set.seed(2)
#' n <- 300
#' d <- data.frame(pov = rbinom(n, 1, 0.4), age = runif(n, 14, 18))
#' d$lead <- rbinom(n, 1, plogis(-1 + 1.5 * d$pov))
#' d$y <- 0.5 * d$lead + d$pov + 0.2 * d$age + rnorm(n)
#' m <- cemtch(d, "lead", c("pov", "age"), cut = list(age = 4))
#' m
#' stats::lm(y ~ lead, d, weights = m$w, subset = m$keep)
#' @export
cemtch <- function(data, d, x, cut = NULL, lcut = NULL) {
  f <- "cemtch"
  dd <- .axdat(data, f)
  if (nrow(dd) < 2L) .ckstop(f, "`data` needs two rows at least")
  .axcols(dd, d, "d", f, n = 1L)
  .axcols(dd, x, "x", f)
  tr <- .axmat(dd, d, f)[, 1L]
  if (!all(tr %in% c(0, 1))) .ckstop(f, "column `%s` must hold 0 or 1", d)
  a <- .cemcode(dd, x, cut, "cut", grDevices::nclass.Sturges, f)
  b <- .cemcode(dd, x, lcut, "lcut", grDevices::nclass.scott, f)
  k <- .Call(C_clab_r_cemtch, a$code, as.integer(tr), b$code)
  st <- k$st
  s <- k$strata
  s[s == 0L] <- NA_integer_
  keep <- !is.na(s)
  g <- factor(ifelse(tr == 1, "treated", "control"), c("control", "treated"))
  tab <- rbind(all = table(g), matched = table(g[keep]),
               unmatched = table(g[!keep]))
  structure(
    list(keep = keep, strata = s, w = k$w, tab = tab,
         l1 = c(before = st[1L], after = st[2L]),
         nstrata = c(kept = as.integer(st[3L]), all = as.integer(st[6L])),
         cut = a$cut),
    class = "cemtch"
  )
}

#' @rdname cemtch
#' @param ... Not used by `print()`.
#' @export
print.cemtch <- function(x, ...) {
  cat("cemtch: coarsened exact matching, ", x[["nstrata"]][["kept"]],
      " of ", x[["nstrata"]][["all"]], " strata kept\n", sep = "")
  print(x[["tab"]])
  cat("L1 imbalance: ", format(x[["l1"]][["before"]], digits = 3),
      " before, ", format(x[["l1"]][["after"]], digits = 3), " after\n",
      sep = "")
  invisible(x)
}

# The covariates of x as integer codes (n x k) and the breaks of the
# numeric ones: cut[[v]] a number of bins or the breaks, else rule(v)
# bins.
.cemcode <- function(d, x, cut, what, rule, f) {
  if (!is.null(cut) && (!is.list(cut) || is.null(names(cut)) ||
                          !all(names(cut) %in% x))) {
    .ckstop(f, "`%s` must be a list named by columns of `x`", what)
  }
  br <- list()
  m <- matrix(0L, nrow(d), length(x))
  for (j in seq_along(x)) {
    z <- d[[x[j]]]
    if (!is.numeric(z)) {
      m[, j] <- match(as.character(z), unique(as.character(z)))
      next
    }
    if (!all(is.finite(z))) .ckstop(f, "column `%s` must be finite", x[j])
    br[[x[j]]] <- .cembr(z, cut[[x[j]]], rule, x[j], what, f)
    m[, j] <- findInterval(z, br[[x[j]]], rightmost.closed = TRUE,
                           all.inside = TRUE)
  }
  list(code = m, cut = br)
}

.cembr <- function(z, k, rule, v, what, f) {
  r <- range(z)
  if (is.null(k)) k <- rule(z)
  if (length(k) == 1L) {
    if (!is.numeric(k) || !is.finite(k) || k < 1 || k != round(k)) {
      .ckstop(f, "`%s$%s` must be a whole number of bins from 1", what, v)
    }
    if (r[1L] == r[2L]) return(r[1L] + c(-0.5, 0.5))
    return(seq(r[1L], r[2L], length.out = k + 1))
  }
  if (!is.numeric(k) || anyNA(k) || is.unsorted(k, strictly = TRUE) ||
        k[1L] > r[1L] || k[length(k)] < r[2L]) {
    .ckstop(f, "`%s$%s` must be increasing breaks covering the column", what,
            v)
  }
  k
}
