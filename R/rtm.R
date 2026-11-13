#' Circle hypothesis of a crime series
#'
#' Canter and Larkin's (1993) circle hypothesis: draw the circle whose
#' diameter joins the two crime sites of a series that lie farthest apart.
#' An offender whose home lies within it is a "marauder" (who travels out
#' from home), one whose home lies outside a "commuter". Gives the circle,
#' the share of the sites it holds (all of them, unless the series is far
#' from round) and, for each home given, whether it lies within.
#'
#' Of several pairs equally far apart, the first in the order of the sites
#' is used. A point lies within the circle when its distance from the
#' centre is at most \eqn{r + 10^{-12} \max(1, r)}.
#'
#' @param x,y The coordinates of the crime sites: finite numbers, at least
#'   two sites, in a planar projection (metres, feet).
#' @param home `NULL`, or the home of the offender: `c(x, y)`, or a matrix of
#'   two columns with one row per candidate home.
#' @return A list: `centre` (`x` and `y`), `r` the radius, `pair` the two
#'   sites (row numbers) on the diameter, `share` the share of the sites
#'   within the circle, and `mar`, `NULL` or for each home `TRUE` when it
#'   lies within (a marauder).
#' @references
#'   Canter, D. and Larkin, P. (1993). The environmental range of serial
#'   rapists. *Journal of Environmental Psychology*, 13(1), 63-69.
#' @examples
#' crchyp(c(0, 4, 2), c(0, 0, 1), home = c(2, -1))
#' @export
crchyp <- function(x, y, home = NULL) {
  .Call(C_clab_r_crchyp, x, y, home)
}

#' Hit score of a geographic profile
#'
#' A geographic profile ranks the cells of a search area; searched from the
#' highest score down, the hit score is the share of the area searched
#' before the offender's home is reached (Rossmo 2000), the smaller the
#' better. Cells scored the same as the home cell are counted both ways:
#' `hit` searches them after the home, `hitt` before it.
#'
#' @param s The scores of the cells: finite numbers, higher searched first
#'   (a profile from [jtsurf()], say).
#' @param home The cells holding the home: cell numbers from 1 to
#'   `length(s)`, one or several (several offenders, or candidate homes).
#' @return A data frame with one row per home: `home`, `above` (the cells
#'   scored higher), `ties` (the cells scored the same, the home included),
#'   `rank` (`above + 1`), and the hit scores in per cent, `hit` (`100 *
#'   above / length(s)`) and `hitt` (`100 * (above + ties) / length(s)`).
#' @references
#'   Rossmo, D. K. (2000). *Geographic Profiling*. CRC Press.
#' @examples
#' hitscr(c(0.9, 0.5, 0.7, 0.5, 0.1), home = 2)
#' @export
hitscr <- function(s, home) {
  .Call(C_clab_r_hitscr, s, home)
}

#' Risk terrain model
#'
#' Risk terrain modelling (Caplan, Kennedy and Miller 2011) relates the
#' crimes counted in the cells of a grid to risk layers, the presence or
#' density of features (bars, transit stops, vacant lots) around each cell,
#' by Poisson regression. The exponentiated coefficient of a layer is its
#' relative risk value; the cells' relative risk combines the layers.
#'
#' The layers are given either as columns (`cell` `NULL`: a matrix, a data
#' frame or a list with one value per cell each) or built here from the
#' coordinates of the features: `cell` holds the cells' centres and `lay`
#' a list of feature coordinates (one two-column matrix per layer), and a
#' cell's value is, by `kind`, 1 when a feature lies within `rad` and 0
#' otherwise (`"prox"`), the number of features within `rad` (`"dens"`), or
#' the Gaussian kernel density of the features with bandwidth `rad`,
#' standardised over the cells to mean 0 and standard deviation 1 (`"kde"`;
#' 0 everywhere when it does not vary). The fit is maximum likelihood by
#' iteratively reweighted least squares, as [stats::glm()] with
#' `family = poisson()` does; a layer that the others explain is aliased,
#' its coefficient `NA`.
#'
#' @param y The crime counts of the cells: numbers from 0, not all 0.
#' @param lay The layers: columns with one value per cell (`cell` `NULL`),
#'   or a list of feature coordinates, a two-column matrix for each layer.
#'   Names of the columns or list elements name the layers.
#' @param cell `NULL`, or the cells' centres: a matrix of two columns with
#'   one row per cell.
#' @param rad With `cell`: the distance within which features count, or the
#'   kernel bandwidth: one positive number.
#' @param kind With `cell`: `"prox"`, `"dens"` or `"kde"`.
#' @inheritParams cntreg
#' @return A list: `coef`, a data frame with one row per term (`term`, `est`
#'   the coefficient, `se`, `z`, `p`, `rrv` the relative risk value
#'   `exp(est)` and its interval `lo`, `hi`); `cells`, a data frame with one
#'   row per cell (`rr` the relative risk, the product of the layers'
#'   relative risk values to the power of the cell's layer values; `score`,
#'   `rr` over its smallest value; `fitted`; `y`); `lay`, the matrix of
#'   layers; `dev` the deviance; `aic`; `conv`, whether the fit converged;
#'   `level`.
#' @references
#'   Caplan, J. M., Kennedy, L. W. and Miller, J. (2011). Risk terrain
#'   modeling: brokering criminological theory and GIS methods for crime
#'   forecasting. *Justice Quarterly*, 28(2), 360-381.
#' @examples
#' rtmfit(c(1, 2, 4, 8, 3, 5), list(bars = c(0, 1, 0, 1, 0, 1),
#'                                  stops = c(0, 0, 1, 1, 1, 0)))
#' # layers from feature locations
#' cell <- as.matrix(expand.grid(x = 1:4, y = 1:4))
#' bars <- rbind(c(1, 1), c(3, 4))
#' y <- c(3, 1, 0, 0, 2, 1, 0, 1, 0, 0, 1, 2, 0, 1, 2, 3)
#' rtmfit(y, list(bars = bars), cell = cell, rad = 1.5, kind = "dens")$coef
#' @export
rtmfit <- function(y, lay, cell = NULL, rad = NULL, kind = "prox",
                   level = 0.95) {
  .Call(C_clab_r_rtmfit, y, lay, cell, rad, kind, level)
}
