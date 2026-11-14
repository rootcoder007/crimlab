#' Internal benchmarking of officers
#'
#' Asks whether an officer stops more people of a group (or searches,
#' arrests, uses force more often) than other officers do in the same
#' circumstances. Officers work different places and hours, so each
#' officer's stops are compared not with all other stops but with the
#' other officers' stops weighted to look like this officer's: the
#' propensity `e(x)` that a stop with circumstances `x` (place, hour, day,
#' reason) is this officer's comes from a logistic regression on `x`
#' (over all stops), and each other officer's stop gets the weight `e / (1
#' - e)`, the weights of the effect on the treated (Ridgeway and MacDonald
#' 2009). The benchmark is the weighted share of the outcome among those
#' stops, and the z statistic compares the officer's share with it, with
#' the variance of a proportion for the officer and, for the benchmark,
#' that of a proportion over its effective number of stops, `(sum w)^2 /
#' sum w^2`. With many officers some will stand out by chance: `q` is the
#' false discovery rate of Benjamini and Hochberg (1995) over the officers
#' tested, and `flag` marks those with `q` below `fdr`. (The paper adds a
#' regression adjustment to the weights, a doubly robust estimate, and
#' estimates the local false discovery rate; the plain weights and the
#' step-up rate are used here.)
#'
#' The comparison is only as good as `x`: circumstances that differ between
#' officers and are not in `x` are left in the difference. An officer with
#' fewer than `minn` stops, or whose propensity model does not converge
#' (circumstances that separate this officer's stops from all others), is
#' not tested (`NA`).
#'
#' @param y The outcome of each stop: 0 or 1 (or `FALSE` or `TRUE`), such
#'   as whether the person stopped belongs to the group of interest.
#' @param off The officer of each stop: labels, two officers at least.
#' @param x `NULL` or the circumstances of each stop: a numeric matrix, a
#'   data frame or list of numeric columns, or one numeric vector, one row
#'   per stop (code a place or an hour as indicator columns).
#' @param minn The fewest stops of an officer to test: one whole number
#'   from 1.
#' @param fdr The false discovery rate at which an officer is flagged: one
#'   number strictly between 0 and 1.
#' @return An object of class `ibench`, a list: `tab`, one row per officer
#'   (`off`, `n` stops, `rate`, the officer's share of the outcome, `bench`,
#'   the benchmark, `ess`, its effective number of stops, `z`, `p`, `q` and
#'   `flag`); `fdr`; `minn`; and `nx`, the number of columns of `x`.
#' @references
#'   Ridgeway, G. and MacDonald, J. M. (2009). Doubly robust internal
#'   benchmarking and false discovery rates for detecting racial bias in
#'   police stops. *Journal of the American Statistical Association*,
#'   104(486), 661-668. \doi{10.1198/jasa.2009.0034}
#'
#'   Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
#'   rate: a practical and powerful approach to multiple testing. *Journal
#'   of the Royal Statistical Society, Series B*, 57(1), 289-300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' @examples
#' set.seed(2)
#' n <- 3000
#' off <- sample(sprintf("o%02d", 1:12), n, TRUE)
#' night <- rbinom(n, 1, ifelse(off %in% c("o01", "o02"), 0.7, 0.3))
#' y <- rbinom(n, 1, plogis(-1 + night + 0.8 * (off == "o03")))
#' ibench(y, off, data.frame(night), minn = 50)
#' @export
ibench <- function(y, off, x = NULL, minn = 20, fdr = 0.05) {
  .Call(C_clab_r_ibench, y, off, x, minn, fdr)
}

#' @export
print.ibench <- function(x, ...) {
  t <- x[["tab"]]
  cat(
    "ibench: ",
    nrow(t),
    " officers, ",
    sum(!is.na(t[["z"]])),
    " tested; flagged at a false discovery rate of ",
    format(x[["fdr"]]),
    ": ",
    sum(t[["flag"]], na.rm = TRUE),
    "\n",
    sep = ""
  )
  print(t, row.names = FALSE)
  invisible(x)
}
