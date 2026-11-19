#' Empirical orthogonal functions of a space-time matrix
#'
#' `eofst()` decomposes a matrix of counts or rates, one row per period
#' and one column per place, into empirical orthogonal functions (EOFs):
#' spatial patterns (the columns of `eofs`) whose weights over time (the
#' principal component series, the columns of `pcs`) reproduce the
#' centred matrix, ordered by the share of its variance each explains. A
#' leading pattern that loads on a cluster of places, with a series that
#' rises, is a hot spot that grows; a second pattern of opposite signs is
#' a shift of crime from one set of places to another.
#'
#' @details Each column is centred on its mean over the periods and the
#'   centred matrix `A` is decomposed as `A = U D V'` (Lorenz 1956;
#'   Hannachi, Jolliffe and Stephenson 2007). The EOFs are the columns of
#'   `V`, the series the columns of `A V = U D`, and `explained` is
#'   `d^2 / sum(d^2)`; `sdev` is `d / sqrt(n - 1)`, the standard deviation
#'   of each series, as [stats::prcomp()] gives it. An EOF is defined up to
#'   its sign: each is signed so that its largest loading in absolute
#'   value is positive. The compiled kernel finds `V` and `d` from the
#'   eigendecomposition (cyclic Jacobi) of the smaller of `A'A` and `A A'`;
#'   components with `d^2` below `1e-10` of the largest are dropped, so a
#'   matrix of `n` periods has at most `n - 1` components.
#'
#' @param z A numeric matrix (or a data frame of numeric columns) with no
#'   missing values: one row per period (at least two), one column per
#'   place.
#' @param k `NULL` (every component) or the most components to keep: one
#'   whole number from 1.
#' @param x For `print()`, an `eofst` result.
#' @inheritParams cemtch
#' @return An object of class `eofst`, a list: `eofs` (places by
#'   components), `pcs` (periods by components), `explained` (the share of
#'   the variance of each component), `sdev` and `center` (the column
#'   means taken out).
#' @references Lorenz, E. N. (1956). *Empirical Orthogonal Functions and
#'   Statistical Weather Prediction*. Statistical Forecasting Project,
#'   Scientific Report 1. Massachusetts Institute of Technology.
#'
#'   Hannachi, A., Jolliffe, I. T. and Stephenson, D. B. (2007). Empirical
#'   orthogonal functions and related techniques in atmospheric science: a
#'   review. *International Journal of Climatology*, 27(9), 1119-1152.
#'   \doi{10.1002/joc.1499}
#' @examples
#' set.seed(3)
#' yr <- 1:10
#' z <- outer(yr, c(1, 1, 1, 0, 0, 0)) + matrix(rnorm(60), 10)
#' colnames(z) <- paste0("ward", 1:6)
#' e <- eofst(z)
#' e
#' round(e$eofs[, 1], 2)
#' @export
eofst <- function(z, k = NULL) {
  f <- "eofst"
  z <- .eozmat(z, f)
  if (!is.null(k) && !.pnum1(k, 1, Inf, TRUE)) {
    .ckstop(f, "`k` must be NULL or one whole number from 1")
  }
  storage.mode(z) <- "double"
  kk <- if (is.null(k)) 0L else as.integer(min(k, dim(z)))
  r <- .Call(C_clab_r_eofst, z, kk)
  nm <- paste0("EOF", seq_along(r$d2))
  dimnames(r$eofs) <- list(colnames(z), nm)
  dimnames(r$pcs) <- list(rownames(z), nm)
  ce <- r$center
  names(ce) <- colnames(z)
  structure(
    list(eofs = r$eofs, pcs = r$pcs, explained = r$d2 / r$tot,
         sdev = sqrt(r$d2) / sqrt(nrow(z) - 1), center = ce),
    class = "eofst"
  )
}

#' @rdname eofst
#' @export
print.eofst <- function(x, ...) {
  cat("eofst: ", length(x[["sdev"]]), " components of ", nrow(x[["pcs"]]),
      " periods by ", nrow(x[["eofs"]]), " places\n", sep = "")
  e <- x[["explained"]]
  print(data.frame(component = seq_along(e), sdev = signif(x[["sdev"]], 4),
                   explained = signif(e, 4),
                   cumulative = signif(cumsum(e), 4)), row.names = FALSE)
  invisible(x)
}

#' The data checked: a finite numeric matrix (a data frame of numbers
#' made one) of at least two rows.
#' @noRd
.eozmat <- function(z, f) {
  if (is.data.frame(z) && length(z) && all(vapply(z, is.numeric, NA))) {
    z <- as.matrix(z)
  }
  if (!is.matrix(z) || !is.numeric(z) || nrow(z) < 2L || ncol(z) < 1L) {
    .ckstop(f, paste("`z` must be a numeric matrix of at least two rows:",
                     "periods by places"))
  }
  if (!all(is.finite(z))) .ckstop(f, "`z` must be finite")
  z
}
