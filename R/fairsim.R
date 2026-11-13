#' Simulated records with a known disparity
#'
#' Incidents in segregated areas, whose true occurrence is the same for
#' every group but whose recording, and the risk an algorithm predicts for
#' them, differ by group by a set amount: data on which an audit's answer
#' is known (after Barman and Barman 2026). The first group is the
#' reference.
#'
#' Each record draws its group by `props`, then an area among those its
#' group lives in (with `g` groups, areas 1, `g + 1`, `2 g + 1`, ... hold
#' the first group only, areas 2, `g + 2`, ... the second, and so on), then
#' `truth`, 1 with probability `base` whatever the
#' group; `det`, the recorded outcome, 1 with probability `base` in the
#' reference group and `base (1 - bias)` in the others (clamped to 0 and
#' 1), so that the other groups' incidents are recorded less when `bias` is
#' positive; and `risk`, normal with mean 250 in the reference group and
#' `250 + 100 bias` in the others and standard deviation 40, clamped to 0
#' and 500. So [fairdf()] on `det` should find a ratio near `1 - bias`,
#' [scrdsp()] on `risk` a gap near `100 bias`, and [calaud()] the other
#' groups' areas ranked as riskier than their recorded outcomes.
#'
#' @param n The number of records: a whole number from 1.
#' @param groups The group labels, at least two, the reference first.
#' @param props `NULL` for equal shares, or one weight from 0 per group (not
#'   all 0).
#' @param nareas The number of areas: a whole number, at least the number
#'   of groups.
#' @param base The true rate, and the reference group's recorded rate: one
#'   number from 0 to 1.
#' @param bias The disparity injected: one number from -1 to 1.
#' @param seed `NULL` to draw the seed from R's generator, or one whole
#'   number from 0; the records come from crimlab's own generator, the
#'   same on every platform.
#' @return A data frame with one row per record: `area` (a number from 1),
#'   `group`, `truth`, `det` and `risk`.
#' @references Barman, P. K. and Barman, P. K. (2026). Spatial disparities
#'   in simulated patrol detection: a reproducible generative-model audit
#'   across Baltimore and Chicago. \doi{10.48550/arXiv.2603.18987}
#' @examples
#' d <- fairsim(2000, seed = 1)
#' head(d)
#' # the injected ratio is 1 - 0.5
#' fairdf(d$det, d$group, base = "A")$tab[, c("group", "di", "dilo", "dihi")]
#' @export
fairsim <- function(n = 2000, groups = c("A", "B"), props = NULL,
                    nareas = 20, base = 0.3, bias = 0.5, seed = NULL) {
  .Call(C_clab_r_fairsim, n, groups, props, nareas, base, bias, seed)
}

#' Noisy-OR detection by patrols
#'
#' The chance that an incident is seen by patrols when each officer within
#' `radius` of it sees it independently with probability `p`: with `k`
#' officers in range, `1 - (1 - p)^k` (the detection model of Barman and
#' Barman 2026). With `draw = TRUE`, one detection outcome is drawn for
#' each incident. A model of who is seen, under stated assumptions; how
#' patrols are actually placed and how much they see are matters of fact
#' the model does not settle.
#'
#' @param x The incidents: a numeric matrix of two columns, or a data frame
#'   with columns `x` and `y` (or two columns), planar coordinates.
#' @param off The officers, in the same forms; `NULL` or no rows for none.
#' @param radius The range: one positive number, in the coordinates' units
#'   (distances are Euclidean; an officer exactly at `radius` counts).
#' @param p The chance one officer in range sees an incident: one number
#'   above 0 and at most 1.
#' @param draw `TRUE` to draw whether each incident is detected.
#' @param seed `NULL` to draw the seed of the detections from R's
#'   generator, or one whole number from 0; read only with `draw`.
#' @return A data frame with one row per incident: `k` (officers in range),
#'   `prob` and, with `draw`, `det` (1 detected, 0 not).
#' @examples
#' set.seed(1)
#' inc <- matrix(runif(20), ncol = 2)
#' off <- matrix(runif(10), ncol = 2)
#' nordet(inc, off, radius = 0.3, p = 0.7)
#' @export
nordet <- function(x, off, radius, p = 0.85, draw = FALSE, seed = NULL) {
  .Call(C_clab_r_nordet, x, off, radius, p, draw, seed)
}
