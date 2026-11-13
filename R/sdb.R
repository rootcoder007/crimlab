#' Hot-spot dynamics: the Short-D'Orsogna-Brantingham model
#'
#' Simulates the model of residential burglary of Short et al. (2008) on
#' an `n` by `n` grid of houses that wraps at its edges. Each site has an
#' attractiveness `a0 + b`; `b` rises by `theta` with every burglary
#' there, spreads to the four neighbours at rate `eta` and decays at rate
#' `omega`, so that burglaries draw further burglaries (repeat and
#' near-repeat victimisation). In the `"agent"` model offenders
#' (`rho`, whole numbers per site) burgle with probability
#' \eqn{1 - e^{-(a_0 + b) dt}} and leave, or else move to a neighbouring
#' site chosen in proportion to its attractiveness; new offenders arrive
#' at `gam` per site per unit time. In the `"pde"` model the densities
#' follow the continuum limit, \eqn{b_t = (\eta D / z) \Delta b - \omega b
#' + \epsilon D \rho A} and \eqn{\rho_t = (D / z)(\Delta \rho - 2 \nabla
#' \cdot (\rho \nabla A / A)) - \rho A + \gamma} with \eqn{A = a_0 + b},
#' by explicit Euler steps on a grid of spacing `dx`; `sup` marks sites
#' where policing suppresses crime (no burglaries there). The step must
#' satisfy `dt <= z dx^2 / (4 dif max(1, eta))`.
#'
#' Started at the uniform steady state [sdbthr()] gives, the `"pde"`
#' model stays there; started near it, a mode grows or dies at the rate
#' [sdbthr()] predicts for its wavenumber. Without offenders the
#' attractiveness only decays: each step multiplies its largest value by
#' at most `1 - omega dt` (`proof("sdbsim")`) and its mean by exactly
#' that.
#'
#' @param n The side of the grid: one whole number from 3 (at most 2000).
#' @param steps The number of steps: one whole number from 0.
#' @param mode `"agent"` or `"pde"`.
#' @param a0 The static attractiveness: one positive number.
#' @param omega The decay rate of `b`: one number from 0 (`omega dt` at
#'   most 1).
#' @param eta The spread of `b` to neighbours: one number from 0 to 1.
#' @param gam The arrival rate of offenders per site: one number from 0.
#' @param theta The rise of `b` per burglary (`"agent"`): one number from
#'   0.
#' @param dt The time step: one positive number; `NULL` for 1 (`"agent"`)
#'   or 0.01 (`"pde"`).
#' @param b,rho `NULL` or `n` by `n` matrices of numbers from 0: the
#'   starting attractiveness (0 by default) and offenders (`"agent"`:
#'   whole numbers, none by default; `"pde"`: a density, `gam / a0` by
#'   default).
#' @param dx,dif,eps,z The grid spacing, the offenders' diffusion, the
#'   rise of `b` per unit of crime and the lattice's neighbour count of
#'   the `"pde"` model: positive numbers (`eps` from 0).
#' @param sup `NULL` or an `n` by `n` logical matrix of suppressed sites
#'   (`"pde"`).
#' @param seed `NULL` to draw the seed from R's generator, or one whole
#'   number from 0 (`"agent"`).
#' @return An object of class `sdbsim`, a list: `b` and `rho`, the fields
#'   at the end; `tot`, the burglaries at each site over the run (the
#'   integrated crime rate for `"pde"`); `burg` and `meanb`, the
#'   burglaries and the mean of `b` at each step; and `mode`.
#' @references
#'   Short, M. B., D'Orsogna, M. R., Pasour, V. B., Tita, G. E.,
#'   Brantingham, P. J., Bertozzi, A. L. and Chayes, L. B. (2008). A
#'   statistical model of criminal behavior. *Mathematical Models and
#'   Methods in Applied Sciences*, 18(supp01), 1249-1267.
#'   \doi{10.1142/S0218202508003029}
#' @examples
#' s <- sdbsim(20, 50, a0 = 1 / 30, omega = 1 / 15, eta = 0.03,
#'             gam = 0.02, theta = 0.5, seed = 1)
#' s
#' p <- sdbsim(16, 200, mode = "pde", eps = 2, dt = 0.05,
#'             rho = matrix(0.07, 16, 16))
#' range(p$b)
#' @export
sdbsim <- function(
  n,
  steps,
  mode = "agent",
  a0 = 1 / 30,
  omega = 1 / 15,
  eta = 0.03,
  gam = 0.002,
  theta = 1,
  dt = NULL,
  b = NULL,
  rho = NULL,
  dx = 1,
  dif = 1,
  eps = 1,
  z = 4,
  sup = NULL,
  seed = NULL
) {
  if (is.null(dt)) dt <- if (identical(mode, "pde")) 0.01 else 1
  .Call(C_clab_r_sdbsim, n, steps, mode, a0, omega, eta, gam, theta, dt, b,
        rho, dx, dif, eps, z, sup, seed)
}

#' @export
print.sdbsim <- function(x, ...) {
  k <- length(x[["burg"]])
  cat("sdbsim: ", x[["mode"]], " model, ", nrow(x[["b"]]), " x ",
      ncol(x[["b"]]), " grid, ", k, " steps; burglaries ",
      format(sum(x[["burg"]])), "; mean attractiveness at the end ",
      if (k) format(x[["meanb"]][k]) else format(mean(x[["b"]])), "\n",
      sep = "")
  invisible(x)
}

#' Hot-spot threshold of the Short-D'Orsogna-Brantingham model
#'
#' The continuum model of [sdbsim()] has one uniform steady state,
#' \eqn{\bar b = \epsilon D \gamma / \omega}, \eqn{\bar A = a_0 + \bar b},
#' \eqn{\bar\rho = \gamma / \bar A}. Linearised about it, a disturbance of
#' wavenumber-squared \eqn{k} evolves by a 2 by 2 matrix whose trace is
#' negative at every \eqn{k} and whose determinant is \eqn{P k^2 - Q k +
#' R}, \eqn{R = \omega \bar A > 0}: a mode grows exactly when the
#' determinant is negative. So hot spots form when \eqn{Q > 0}, that is
#' \eqn{3 \epsilon D \bar\rho > \omega + \eta \bar A} (repeat
#' victimisation outweighs decay and spreading), and \eqn{Q^2 > 4 P R};
#' the fastest-growing mode is near \eqn{k^* = Q / (2 P)}. Without
#' reinforcement (`eps = 0`) the uniform state is stable at every
#' wavelength (`proof("sdbthr")`). `z` places the lattice's neighbour
#' count on both diffusion terms, as [sdbsim()] does (`z = 1` is the
#' continuum form).
#'
#' @param a0,omega,eps,dif,gam,eta,z The model's constants, as for
#'   [sdbsim()]: one number each, `eps` and `eta` from 0, the others
#'   positive.
#' @param k `NULL` or wavenumbers-squared (numbers from 0) at which to
#'   tabulate the trace, the determinant and the growth rate.
#' @return An object of class `sdbthr`, a list: `steady` (`b`, `a`,
#'   `rho`); `p`, `q`, `r`, the determinant's coefficients; `thr`, the two
#'   sides of the threshold (`lhs` \eqn{3 \epsilon D \bar\rho}, `rhs`
#'   \eqn{\omega + \eta \bar A}); `unstable`; `kcrit`, `detcrit`,
#'   `grcrit`, the critical mode, its determinant and growth rate (`NA`
#'   unless `p` and `q` are positive); and `disp`, `NULL` or a data frame
#'   (`k`, `trace`, `det`, `growth`).
#' @references
#'   Short, M. B., Bertozzi, A. L. and Brantingham, P. J. (2010). Nonlinear
#'   patterns in urban crime: hotspots, bifurcations, and suppression.
#'   *SIAM Journal on Applied Dynamical Systems*, 9(2), 462-483.
#'   \doi{10.1137/090759069}
#' @examples
#' sdbthr(a0 = 1, omega = 1, eps = 2, dif = 1, gam = 1, eta = 0.01)
#' sdbthr(a0 = 1, omega = 1, eps = 0, dif = 1, gam = 1, eta = 0.01)$unstable
#' @export
sdbthr <- function(a0, omega, eps, dif, gam, eta, z = 1, k = NULL) {
  .Call(C_clab_r_sdbthr, a0, omega, eps, dif, gam, eta, z, k)
}

#' @export
print.sdbthr <- function(x, ...) {
  cat("sdbthr: steady b ", format(x[["steady"]][["b"]]), ", rho ",
      format(x[["steady"]][["rho"]]), "; 3 eps dif rho = ",
      format(x[["thr"]][["lhs"]]), " against omega + eta A = ",
      format(x[["thr"]][["rhs"]]), "\n",
      if (x[["unstable"]]) "hot spots form" else "the uniform state is stable",
      if (!is.na(x[["kcrit"]])) {
        paste0("; critical k ", format(x[["kcrit"]]))
      },
      "\n", sep = "")
  invisible(x)
}
