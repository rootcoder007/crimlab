#' Accuracy of a hot-spot map
#'
#' Asks how well a map of predicted hot spots (from [kdsurf()],
#' [prohs()], [rtmfit()] or any score per cell) caught the crimes that
#' followed. The hot cells are those of highest predicted score whose
#' areas together first reach the share `top` of the study area (or the
#' cells marked `TRUE` when `pred` is logical). Then:
#'
#' * the hit rate `n / N`, the share of the crimes of the evaluation
#'   period (`obs`) in hot cells, and the predictive accuracy index `PAI =
#'   (n / N) / (a / A)`, the hit rate over the hot share of the area
#'   (Chainey, Tompson and Uhlig 2008);
#' * the predictive efficiency index `PEI = n / nmax`, the crimes caught
#'   over the most an area of the same size could have caught, its cells
#'   taken by observed crimes per unit area, of equal density the larger
#'   first (Drawve 2016): 1 for a perfect map of that size;
#' * with `prev`, the counts of the period the map was made from, the
#'   recapture rate index of Levine (2008), `(n / n_prev) / (N / N_prev)`:
#'   whether the hot cells kept their crimes better than the whole area.
#'
#' @param pred The predicted score of each cell (higher is hotter), or
#'   `TRUE` for the cells of a map's hot spots.
#' @param obs The crimes of each cell in the evaluation period: numbers from
#'   0, not all 0.
#' @param area `NULL` for cells of equal area, or the area of each cell:
#'   positive.
#' @param top The shares of the area to mark hot: numbers above 0 and at
#'   most 1, one row of the result each (not read when `pred` is logical).
#' @param prev `NULL` or the crimes of each cell in the period the map was
#'   made from: numbers from 0.
#' @param ties Which of the cells whose score ties across the cut-off of
#'   `top` are hot: `"in"` (the default) all of them, so that the hot area
#'   may pass `top`; `"out"` none of them, so that it may fall short; or
#'   `"frac"` each by the share of their area that the cut-off still needs,
#'   so that the hot area is `top` and `cells`, `area`, `n` and `nmax` may
#'   be fractional. The sources cited do not say how ties are broken. With
#'   any of the three the result does not depend on the order of the
#'   cells; with no tie across the cut-off all three give the first cells
#'   by score whose areas reach `top`, as before `ties` was added (when ties
#'   were broken by position, so that reordering the cells could change
#'   the result).
#' @return A data frame with one row per `top`: `top`, `cells`, `area`
#'   (hot), `ashare`, `n`, `N`, `hit`, `pai`, `nmax`, `pei` and `rri` (`NA`
#'   without `prev`).
#' @references
#'   Chainey, S., Tompson, L. and Uhlig, S. (2008). The utility of hotspot
#'   mapping for predicting spatial patterns of crime. *Security Journal*,
#'   21(1-2), 4-28. \doi{10.1057/palgrave.sj.8350066}
#'
#'   Levine, N. (2008). The "hottest" part of a hotspot: comments on "The
#'   utility of hotspot mapping for predicting spatial patterns of crime".
#'   *Security Journal*, 21(4), 295-302. \doi{10.1057/sj.2008.5}
#'
#'   Drawve, G. (2016). A metric comparison of predictive hot spot
#'   techniques and RTM. *Justice Quarterly*, 33(3), 369-397.
#'   \doi{10.1080/07418825.2014.904393}
#' @examples
#' set.seed(12)
#' risk <- rgamma(400, 0.5)
#' before <- rpois(400, risk)
#' after <- rpois(400, risk)
#' hsacc(pred = before, obs = after, top = c(0.05, 0.1, 0.2), prev = before)
#' @export
hsacc <- function(pred, obs, area = NULL, top = 0.1, prev = NULL,
                  ties = "in") {
  .Call(C_clab_r_hsacc, pred, obs, area, top, prev, ties)
}
