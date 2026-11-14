#' Crosswalk between boundary versions
#'
#' Carries values kept on one set of areas (the 158 neighbourhoods of a
#' city's new boundaries, a census vintage's tracts) to another (the 140
#' old neighbourhoods) through a correspondence table: each link says what
#' share `w` of a source area's value falls in a target area. Counts of a
#' source wholly inside one target move whole (weight 1, exact); a source
#' split between targets is shared by area or population (uniform density
#' within the source, an assumption). The same table read the other way
#' disaggregates. `map` also gives each source's main target, the column to
#' add to a table keyed by the source areas.
#'
#' @param from,to The correspondence: the source and target area of each
#'   link, labels. A `clpre` list from [clprep()] is taken too (as `from`):
#'   its elements fill the arguments by name.
#' @param x The values to carry: finite numbers, one per record.
#' @param unit The source area of each record: labels (as in `from`).
#' @param w `NULL` (every link weight 1), or the share of each link: numbers
#'   from 0.
#' @return A list: `tab`, a data frame with one row per target (`to`,
#'   `value`); `map`, one row per source (`from`; `to`, the target of its
#'   largest share, the first of ties; `share`, that share; `wsum`, the sum
#'   of its shares, 1 when nothing is lost or counted twice); `lost`, the
#'   total of the records whose area is not in `from`; `unmatched`, those
#'   areas.
#' @examples
#' # old areas a, b merge into A; c splits between A and B by area
#' bndxw(from = c("a", "b", "c", "c"), to = c("A", "A", "A", "B"),
#'       x = c(3, 4, 5, 8), unit = c("a", "b", "b", "c"),
#'       w = c(1, 1, 0.25, 0.75))
#' @export
bndxw <- function(from, to, x, unit, w = NULL) {
  if (inherits(from, "clpre")) {
    return(.prep(bndxw, from, match.call(), parent.frame()))
  }
  .Call(C_clab_r_bndxw, from, to, x, unit, w)
}

#' Synthetic exposure of areas from a survey
#'
#' A small-area synthetic estimate (Gonzalez and Hoza 1978): the prevalence of a
#' trait (serious mental illness, say, for the benchmark of a disparity in
#' police use of force) is modelled on survey microdata by logistic
#' regression ([lgtreg()]); applied to each area's covariates it gives the
#' area's rate, and times its population the number of people with the
#' trait, the exposure a rate per person needs. The area covariates must be
#' on the scale of the survey's (shares of an area for indicators read as
#' the probability of the indicator, a plug-in that assumes the model
#' holds within areas).
#'
#' @param y The trait of each survey respondent: 0 or 1, both present.
#' @param x The respondents' covariates: a numeric matrix or data frame,
#'   one row per respondent, more rows than columns plus one.
#' @param ax The areas' covariates: a numeric matrix or data frame with the
#'   columns of `x`, one row per area.
#' @param pop The population of each area: numbers from 0.
#' @return A list: `coef`, the logistic fit (`term`, `est`, `se`); `areas`,
#'   a data frame with one row per area (`rate`, `pop`, `expo`); `total`,
#'   the summed exposure; `dev`, the fit's deviance.
#' @references
#'   Gonzalez, M. E. and Hoza, C. (1978). Small-area estimation with
#'   application to unemployment and housing estimates. *Journal of the
#'   American Statistical Association*, 73(361), 7-15.
#'   \doi{10.1080/01621459.1978.10479991}
#' @examples
#' set.seed(1)
#' s <- data.frame(age = runif(400, 18, 80), male = rbinom(400, 1, 0.5))
#' y <- rbinom(400, 1, plogis(-3 + 0.02 * s$age + 0.4 * s$male))
#' areas <- data.frame(age = c(35, 52), male = c(0.48, 0.51))
#' synexp(y, s, areas, pop = c(12000, 8000))$areas
#' @export
synexp <- function(y, x, ax, pop) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(ax)) ax <- as.matrix(ax)
  .Call(C_clab_r_synexp, y, x, ax, pop)
}
