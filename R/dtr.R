#' Offending chosen under certainty and severity of sanction
#'
#' An offender picks the level of offending `x` from a menu to maximise the
#' benefit less the expected sanction, `ben - p * sanc`, with the sanction
#' increasing in `x`. For any shape of the benefit (no concavity, no
#' smoothness) the chosen level never rises with the certainty `p`, and
#' never rises when the sanction schedule grows faster: convexity buys a
#' unique and smooth response, not its direction. Gives the chosen level at
#' each certainty; of tied optima (within `1e-12`), the largest.
#'
#' @param x The levels of offending on the menu: finite numbers, no level
#'   twice (any order).
#' @param ben The benefit of each level: finite numbers.
#' @param sanc The sanction at each level: finite numbers, increasing with
#'   `x`.
#' @param p The certainties of sanction to evaluate: numbers from 0.
#' @return A data frame with one row per certainty: `p`, `xopt` (the level
#'   chosen) and `value` (its net benefit).
#' @references
#'   Becker, G. S. (1968). Crime and punishment: an economic approach.
#'   *Journal of Political Economy*, 76(2), 169-217.
#' @examples
#' x <- 0:10
#' detrsp(x, ben = 3 * sqrt(x) - (x %% 3 == 0), sanc = x^1.5,
#'        p = c(0.1, 0.3, 0.5, 1))
#' @export
detrsp <- function(x, ben, sanc, p) {
  .Call(C_clab_r_detrsp, x, ben, sanc, p)
}

#' Which deterrence effects a design identifies
#'
#' For the response `y = b0 + b1 p + b2 s + b3 c` to the certainty `p`,
#' severity `s` and celerity (delay) `c` of the sanction regimes observed,
#' a dimension that never varies has no identified coefficient, and the
#' coefficients are identified exactly when the design `(1, p, s, c)` has
#' full rank. Gives the rank and, for each dimension, whether its effect is
#' identified given the others (full rank, or the rank falls without it).
#' With `log = TRUE` the design is `(1, log p, log s, c)`: the expected
#' cost of an offence `p s exp(-rho c)` is linear on that scale, and with a
#' response `y` the ratio of the celerity to the certainty coefficient,
#' `-b3 / b1`, is the discount rate `rho` the response implies.
#'
#' @param p,s,c The certainty, severity and delay of each regime: finite
#'   numbers, one per regime (`p` and `s` positive with `log = TRUE`).
#' @param y `NULL`, or the response of each regime: finite numbers (the
#'   design must then have full rank).
#' @param log `TRUE` for the log design, `FALSE` for the plain one.
#' @param rho `NULL`, or a discount rate at which to give the expected cost
#'   of each regime: one finite number.
#' @return A list: `rank`; `ident` and `const`, for `certainty`,
#'   `severity` and `celerity`, whether the effect is identified and
#'   whether the dimension is constant; `coef`, `NULL` or the least-squares
#'   coefficients; `disc`, `NULL`, or the implied discount rate (`NA` but
#'   for the log design); `cost`, `NULL` or `p * s * exp(-rho * c)`.
#' @references
#'   Nagin, D. S. (2013). Deterrence in the twenty-first century. *Crime
#'   and Justice*, 42(1), 199-263.
#' @examples
#' # certainty varied at fixed severity and celerity
#' dtrid(p = c(0.1, 0.3, 0.5), s = c(2, 2, 2), c = c(30, 30, 30))
#' set.seed(2)
#' p <- runif(40, 0.1, 0.9)
#' s <- runif(40, 1, 10)
#' c <- runif(40, 0, 12)
#' y <- 2 - 0.8 * log(p) - 0.8 * log(s) + 0.08 * c + rnorm(40, sd = 0.1)
#' dtrid(p, s, c, y, log = TRUE)$disc
#' @export
dtrid <- function(p, s, c, y = NULL, log = FALSE, rho = NULL) {
  .Call(C_clab_r_dtrid, p, s, c, y, log, rho)
}

#' The criminal, citizen and inspector game
#'
#' Replicator dynamics of a rock-paper-scissors game of three strategies
#' (criminals, citizens, inspectors), each beating the next with stake
#' `a`: `x_i' = a x_i (x_{i+1} - x_{i+2})`. The centroid is the only
#' interior rest point and `x1 x2 x3` is conserved, so the shares circle
#' the centroid for ever and more inspection does not steadily lower crime
#' (D'Orsogna and Perc 2015). Integrated by fourth-order Runge-Kutta; the
#' drift of the invariant measures the integration error.
#'
#' @param a The stake: one finite number, not 0.
#' @param x The starting shares: three numbers from 0 summing to 1.
#' @param steps The Runge-Kutta steps: one whole number from 0, at most
#'   1e7.
#' @param dt The step size: one positive number.
#' @return A list: `traj`, the shares after each step (`steps + 1` rows);
#'   `inv`, the invariant `x1 x2 x3` along the path; `drift`, its largest
#'   departure from the start; `dvdt`, its exact derivative at the start
#'   (0); `field`, the vector field at the start; `fixed`, the centroid.
#' @references
#'   D'Orsogna, M. R. and Perc, M. (2015). Statistical physics of crime: a
#'   review. *Physics of Life Reviews*, 12, 1-21.
#' @examples
#' g <- inspg(1, c(0.5, 0.3, 0.2), steps = 200)
#' g$drift
#' @export
inspg <- function(a, x, steps = 1000, dt = 0.01) {
  .Call(C_clab_r_inspg, a, x, steps, dt)
}

#' Phase diagram of the inspection game
#'
#' The inspection game of Helbing et al. (2010) with cooperators (C),
#' criminals (P) and inspectors (O), against C, P, O: C earns (1, 0, 1),
#' P earns (T, 0, T - 1) for the temptation `T`, O earns `1 - g` whatever
#' it meets for the inspection cost `g`. Discrete replicator steps
#' `x += rate * x * (f - mean f)` (negative shares set to 0, the shares
#' normalised) run from `x0` at each point of a grid of temptations and
#' costs; the final shares map the phases of the game, the criminals'
#' share standing for the crime rate.
#'
#' @param temp,cost The temptations and inspection costs of the grid:
#'   finite numbers.
#' @param steps The replicator steps at each point: one whole number from
#'   0.
#' @param rate The step size: one positive number.
#' @param x0 The starting shares of C, P and O: three numbers from 0 with a
#'   positive sum.
#' @return A list: `crime`, `coop` and `insp`, the final shares of P, C and
#'   O (matrices, one row per temptation and one column per cost); `temp`
#'   and `cost`.
#' @references
#'   Helbing, D., Szolnoki, A., Perc, M. and Szabo, G. (2010). Punish, but
#'   not too hard: how costly punishment spreads in the spatial public goods
#'   game. *New Journal of Physics*, 12, 083005.
#' @examples
#' ph <- insphs(seq(0.1, 1.8, by = 0.3), seq(0.1, 1.1, by = 0.25),
#'              steps = 200)
#' round(ph$crime, 2)
#' @export
insphs <- function(temp = seq(0.05, 1.8, length.out = 20),
                   cost = seq(0.05, 1.2, length.out = 20), steps = 600,
                   rate = 0.05, x0 = c(0.34, 0.33, 0.33)) {
  .Call(C_clab_r_insphs, temp, cost, steps, rate, x0)
}

#' Collider bias of arrest data
#'
#' Two causes of arrest with prevalences `a` and `e`, independent in the
#' population, and a background cause of prevalence `bg`; a person is
#' arrested when any is present. Among arrestees the two causes are no
#' longer independent: their odds ratio is exactly `bg`, below 1, and the
#' negative association grows as arrests for other reasons become rare.
#' Studying arrestees alone manufactures it (Hernan and Robins, Fine Point
#' 8.2).
#'
#' @param a,e The prevalences of the two causes: numbers strictly between 0
#'   and 1, one per scenario or one for all.
#' @param bg The prevalence of the background cause: numbers from 0 to 1,
#'   one per scenario or one for all.
#' @return A data frame with one row per scenario: `a`, `e`, `bg`, `or` the
#'   odds ratio among arrestees, `lor` its log, and the shares of arrestees
#'   with both causes (`c11`), `a` only (`c10`), `e` only (`c01`) and
#'   neither (`c00`).
#' @references
#'   Hernan, M. A. and Robins, J. M. (2020). *Causal Inference: What If*.
#'   Chapman and Hall/CRC.
#' @examples
#' colarr(a = 0.3, e = 0.2, bg = c(0.5, 0.1, 0.01))
#' @export
colarr <- function(a, e, bg) {
  .Call(C_clab_r_colarr, a, e, bg)
}
