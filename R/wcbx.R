#' Inference for a difference in differences with few clusters
#'
#' Tests of the two-way fixed-effects effect of [didfe()] that stay
#' accurate when a policy was assigned to a handful of clusters (districts,
#' states), where the cluster-robust t test over-rejects: the restricted
#' wild cluster bootstrap, and the permutation of the adoption dates across
#' clusters; and the power of a design simulated on real untreated data.
#'
#' `meth = "wild"` is the restricted wild cluster bootstrap-t (WCR) of
#' Cameron, Gelbach and Miller (2008): the model is fitted with the null
#' `b = b0` imposed, each bootstrap outcome is that fit plus its residuals
#' multiplied cluster by cluster by a random sign (`wtype = "rad"`) or one
#' of Webb's six values (`"webb"`), and the CR1 t statistic of [didfe()] is
#' recomputed. The p-value is the share of draws with `|t*|` at least `|t|`
#' (`ties = "count"`: the sign patterns that reproduce the data are among
#' the draws, as the observed assignment is in a permutation test) or above
#' it (`"strict"`, as 'fwildclusterboot' counts). With Rademacher weights
#' and `2^G <= reps` every sign pattern is enumerated. The interval collects
#' the `b0` not rejected at `1 - level`, by bisection with the same draws;
#' with Rademacher weights the smallest attainable p-value is `2 / 2^G`, so
#' with four clusters nothing is rejected at 5% and the interval is the
#' whole line: the honest answer.
#'
#' `meth = "perm"` shuffles the adoption period of each cluster (all its
#' units adopt together; never is a period too) among the clusters,
#' keeping the number adopting at each date, and recomputes the t statistic
#' (`stat = "t"`) or the estimate (`"est"`). When the distinct assignments
#' number at most `reps` they are all enumerated and `p` is the exact share
#' at least as extreme; otherwise `p = (b + 1) / (R + 1)` over `R` random
#' ones. The test is Fisher's: valid in any sample when the dates were
#' assigned at random. With several outcomes (the columns of `y`), `pstep`
#' holds the Westfall-Young stepdown p-values of [wystep()] over the same
#' assignments.
#'
#' `meth = "power"` draws `ntrt` clusters, gives each an adoption period
#' drawn from `adopt`, adds each effect of `eff` to their outcomes from
#' adoption on, and counts the CR1 t tests that reject at `alpha`; on real
#' untreated data the serial correlation and clustering are carried, as
#' formula-based power calculations do not (Burlig, Preonas and Woerman
#' 2020). `size` is the rejection rate at effect 0 (it should be near
#' `alpha`); `mde`, where the power curve first crosses `power`, by linear
#' interpolation.
#'
#' @inheritParams didfe
#' @param y The outcome: numbers, one per row; for `meth = "perm"` also a
#'   matrix or data frame of several outcomes, one column each.
#' @param trt The treatment of each row, 0 or 1 (`"wild"`, `"perm"`; for
#'   `"perm"` it must stay on, from one period within each cluster); `NULL`
#'   for `"power"`.
#' @param clus `NULL` (cluster by unit) or the cluster of each row.
#' @param meth `"wild"`, `"perm"` or `"power"`.
#' @param reps The draws: a whole number from 9; `NULL` for 9999 (200 for
#'   `"power"`).
#' @param wtype `"rad"` (Rademacher) or `"webb"`.
#' @param ties `"count"` or `"strict"`.
#' @param stat `"t"` or `"est"`: the statistic permuted.
#' @param ntrt The clusters treated in each power draw.
#' @param adopt The adoption periods drawn from, periods of `time`.
#' @param eff `NULL` (0 to two standard deviations of `y` in five steps) or
#'   the effects tried.
#' @param alpha The level of the power test.
#' @param power The power of the minimum detectable effect.
#' @inheritParams evtst
#' @return An object of class `wcbx`, a list: `res`, for `"wild"` one row
#'   (`est`, `stat`, `p`, `lo`, `hi`), for `"perm"` one row per outcome
#'   (`outcome`, `est`, `stat`, `p`, `pstep`), for `"power"` one row per
#'   effect (`eff`, `power`); `meth`; `ncl`; `draws`; `enum` (`TRUE` when
#'   every draw was enumerated); `mde` and `size` (`"power"`, else `NA`);
#'   and `level`.
#' @references
#'   Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2008).
#'   Bootstrap-based improvements for inference with clustered errors.
#'   *Review of Economics and Statistics*, 90(3), 414-427.
#'   \doi{10.1162/rest.90.3.414}
#'
#'   Webb, M. D. (2023). Reworking wild bootstrap-based inference for
#'   clustered errors. *Canadian Journal of Economics*, 56(3), 839-858.
#'   \doi{10.1111/caje.12661}
#'
#'   Burlig, F., Preonas, L. and Woerman, M. (2020). Panel data and
#'   experimental design. *Journal of Development Economics*, 144, 102458.
#'   \doi{10.1016/j.jdeveco.2020.102458}
#' @examples
#' d <- expand.grid(year = 2015:2024, town = LETTERS[1:8])
#' i <- match(d$town, LETTERS)
#' d$bylaw <- as.integer(i <= 3 & d$year >= 2019)
#' d$rate <- 30 + i + 0.5 * (d$year - 2015) - 2 * d$bylaw + sin(seq_len(80))
#' # eight towns: every one of the 2^8 sign patterns
#' wcbx(d$rate, d$town, d$year, d$bylaw)
#' # every assignment of the by-law to three of the eight towns
#' wcbx(d$rate, d$town, d$year, d$bylaw, meth = "perm")$res
#' @export
wcbx <- function(y, unit, time, trt = NULL, adj = NULL, clus = NULL,
                 meth = "wild", reps = NULL, wtype = "rad", ties = "count",
                 stat = "t", ntrt = 1, adopt = NULL, eff = NULL,
                 alpha = 0.05, power = 0.8, seed = NULL, level = 0.95) {
  .Call(C_clab_r_wcbx, y, unit, time, trt, adj, clus, meth, reps, wtype,
        ties, stat, ntrt, adopt, eff, alpha, power, seed, level)
}

#' @export
print.wcbx <- function(x, ...) {
  cat("wcbx: ", x[["meth"]], ", ", x[["ncl"]], " clusters, ", x[["draws"]],
      if (isTRUE(x[["enum"]])) " draws (all enumerated)\n" else " draws\n",
      sep = "")
  print(x[["res"]], row.names = FALSE)
  if (x[["meth"]] == "power") {
    cat("size ", format(x[["size"]]), ", minimum detectable effect ",
        format(x[["mde"]]), "\n", sep = "")
  }
  invisible(x)
}

#' Westfall-Young stepdown p-values
#'
#' Adjusted p-values for several hypotheses tested together, from the joint
#' null distribution of their statistics (permutations or bootstrap draws),
#' holding the family-wise error rate while using the correlation between
#' the tests.
#'
#' The statistics are ordered by size, `|t_(1)| >= ... >= |t_(K)|`. For each
#' null draw `b`, `q_bj = max_{r >= j} |t*_b(r)|`; the adjusted p-value of
#' the `j`-th is the share of draws with `q_bj >= |t_(j)|`, made monotone:
#' `p~_(j) = max(p~_(j-1), p_(j))` (Westfall and Young 1993, algorithm
#' 2.8). Draws with a missing or infinite statistic are left out.
#'
#' @param obs The observed statistics: finite numbers, one per hypothesis.
#' @param null The statistics under the null: a numeric matrix, one row per
#'   draw and one column per hypothesis.
#' @return A data frame, one row per hypothesis in the order of `obs`:
#'   `stat`, `p` (the share of draws of that statistic at least as large in
#'   absolute value) and `padj` (the stepdown p-value).
#' @references
#'   Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
#'   Testing*. Wiley, New York. ISBN 978-0-471-55761-6.
#' @examples
#' set.seed(1)
#' null <- matrix(rnorm(3000), ncol = 3)
#' wystep(c(2.5, 0.4, -2.1), null)
#' @export
wystep <- function(obs, null) {
  .Call(C_clab_r_wystep, obs, null)
}
