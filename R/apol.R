#' Atmospheric dispersion: sigmas, plume rise, plume, puff, grid and
#' particles
#'
#' @description
#' Six standard models of how a pollutant released into the air spreads
#' downwind, for exposure estimates around a stack, a fire or an industrial
#' site (for example when relating air quality to health or justice
#' outcomes by neighbourhood):
#' \itemize{
#'   \item \code{pgsig()}: how wide and how tall a plume has grown at a
#'     given distance, by atmospheric stability class;
#'   \item \code{briggs()}: how far a hot plume rises above the stack top;
#'   \item \code{gplume()}: the steady concentration at receptor points from
#'     a continuous release;
#'   \item \code{gspuff()}: the concentration from a single instantaneous
#'     release, some time after it;
#'   \item \code{advdif2d()}: a concentration field on a grid moved by wind
#'     and spread by diffusion, step by step;
#'   \item \code{lagrprt()}: the same physics as a random walk of many
#'     particles.
#' }
#' They are screening models with textbook assumptions (flat terrain,
#' steady uniform wind, no chemistry or deposition), not regulatory model
#' runs.
#'
#' In criminology these models serve environmental criminology in its two
#' established senses. In the first, pollution is itself the harm or the
#' offence: the study of environmental crime and harm (green criminology)
#' asks who bears the exposure from an emitting site, and the modelled
#' concentrations by area and group are the exposure behind burden and
#' environmental justice measures (\code{\link{polbrdn}},
#' \code{\link{polbrdarea}}, \code{\link{expconcidx}}). In the second, the
#' environment is a driver of offending and victimisation: studies relating
#' air pollution to crime join an exposure series to crime counts by place
#' and day (Bondy, Roth and Sager 2020, on London wards; Herrnstadt et al.
#' 2019, on Los Angeles), and these models supply the exposure where no
#' monitor stands.
#'
#' @details
#' \code{pgsig()}: Briggs (1973) fits of the Pasquill-Gifford dispersion
#' coefficients by stability class A-F, rural
#' \code{sigmay = a_y x (1 + 0.0001 x)^(-1/2)} (urban: 0.0004) and
#' \code{sigmaz = a_z x (1 + b_z x)^(e_z)}.
#' \code{briggs()}: buoyancy flux \code{F = g v_s d^2 (T_s - T_a) / (4 T_s)};
#' final rise \code{21.425 F^(3/4) / u} with \code{x_f = 49 F^(5/8)} when
#' \code{F < 55}, else \code{38.71 F^(3/5) / u} with \code{x_f = 119 F^(2/5)},
#' for classes A-D; \code{2.6 (F / (u s))^(1/3)} with
#' \code{s = g (dtheta/dz) / T_a} and \code{x_f = 2.0715 u / sqrt(s)} for E
#' and F; transitional rise \code{1.6 F^(1/3) x^(2/3) / u} before \code{x_f}
#' (EPA ISC3).
#' \code{gplume()}:
#' \code{Q / (2 pi u sigmay sigmaz) exp(-y^2 / (2 sigmay^2)) V},
#' with \code{V} the vertical term with ground reflection and, under a
#' mixing lid, \code{nimg} pairs of image terms.
#' \code{gspuff()}: an instantaneous release with \code{sigma_x = sigmay}
#' evaluated at the travel distance \code{u t}.
#' \code{advdif2d()}: explicit upwind advection, central diffusion, forward
#' Euler, zero boundaries. The explicit scheme is stable when
#' \code{|u dt/dx| + |v dt/dy| + 2 (kx dt/dx^2 + ky dt/dy^2) <= 1} (the von
#' Neumann bound for upwind advection with forward-time centred-space
#' diffusion); that sum is returned as \code{stabn} with \code{stable}. A
#' call outside the bound is not run as given: each step is split into
#' enough sub-steps of a smaller \code{dt} to satisfy it (reported as
#' \code{substeps} and \code{dt}), with a warning, so the field is a
#' solution over the same physical time rather than a numerical explosion.
#' The separate \code{cfl} and \code{difn} are reported for reference; each
#' below 1 is necessary, not sufficient.
#' \code{lagrprt()}: the random walk \code{x + u dt + sqrt(2 K dt) xi}
#' (Thomson 1987) with an optional counting grid. The standard normals
#' \code{xi} of step \code{k} (from 0) come from the splitmix64 uniform
#' stream (Steele, Lea and Flood 2014) seeded with
#' \code{(seed * 65536 + s) mod 2147483647}, \code{s = 2k} for x and
#' \code{2k + 1} for y, through the normal quantile; R's own random number
#' stream is not touched.
#'
#' The stability classes run from \code{"A"} (very unstable, sunny with
#' light wind) through \code{"D"} (neutral, overcast or windy) to \code{"F"}
#' (stable, clear night with light wind). Units are SI: metres, seconds,
#' kelvin; with \code{q} in g/s, \code{gplume()} returns g/m^3, and with
#' \code{mass} in g, \code{gspuff()} returns g/m^3. The plume is not defined
#' at or upwind of the source: a receptor with \eqn{x \le 0} gets 0. The
#' plume's \code{h} is the effective height, the stack height plus the rise
#' from \code{briggs()}.
#'
#' Every physical input is checked: emission rates, masses, heights and
#' distances are non-negative, wind speed, travel time and grid steps
#' positive, \code{nimg} a whole number up to 1,000, the stack as hot as
#' the air or hotter (temperatures below 150 are refused as probably
#' Celsius), and \code{x0}, \code{y0} single numbers; a release above
#' \code{mixh} warns, and so does a receptor below ground or above the lid
#' (the reflection terms give mirror values there). Unknown stability
#' classes or settings are errors. Both grid models refuse a run above
#' \code{maxwork} cell or particle updates.
#'
#' @param x Downwind distance(s) in metres.
#' @param stab Pasquill stability class, one of \code{"A"} to \code{"F"}.
#' @param set \code{"rural"} (the default) or \code{"urban"}: which set of
#'   Briggs coefficients to use.
#' @param u Wind speed in m/s (the x velocity for the grid models).
#' @param diam,exitvel,stktmp,ambtmp Stack diameter (m), exit velocity
#'   (m/s), stack gas and ambient temperatures (K). Only buoyant rise is
#'   modelled: a stack at the air's temperature (within half a kelvin) has
#'   zero rise, and momentum rise from the exit velocity is not added.
#' @param dthdz Potential temperature gradient (K/m) for the stable classes;
#'   the defaults are 0.020 (E) and 0.035 (F).
#' @param g Gravitational acceleration in m/s^2, default 9.80616.
#' @param q Emission rate (mass per second, for example g/s).
#' @param h Effective stack height (m): physical height plus plume rise.
#' @param rcpt A matrix with columns x, y, z, or a list of length-3 vectors,
#'   or a vector of x, y, z triples.
#' @param mixh Mixing-lid height in metres, or \code{NULL} (the default) for
#'   no lid. Under a lid the plume is reflected between the ground and the
#'   lid.
#' @param nimg Whole number, default 3: the number of lid image pairs (0 to
#'   1,000). The terms decay quickly; a few are enough.
#' @param mass Puff mass (for example g).
#' @param t Time since release (s).
#' @param c0 Initial concentration grid: a non-negative numeric matrix, one
#'   row per x cell and one column per y cell.
#' @param v Wind velocity along y (m/s), may be negative.
#' @param kx,ky Eddy diffusivities along x and y (m^2/s), non-negative.
#' @param dx,dy Grid spacings (m), positive.
#' @param dt Time step (s), positive.
#' @param nsteps Whole number of time steps (0 returns the start).
#' @param source Source grid (concentration per unit time), or \code{NULL}.
#' @param maxwork The most cell updates (cells x steps x sub-steps; for
#'   \code{lagrprt()}, particles x steps) a run may take, 100 million by
#'   default; \code{Inf} lifts it.
#' @param n Whole number of particles, at most 100 million.
#' @param x0,y0 Release point (m).
#' @param seed Number, default 0: the seed of the particle random walk, a
#'   whole number from -2147483647 to 2^53 - 1.
#' @param grid \code{NULL}, or \code{c(xmin, xmax, ymin, ymax, nx, ny)} to
#'   count particles into a concentration grid.
#' @param normals \code{NULL} (the default stream above), or a function
#'   \code{(n, seed, stream)} returning \code{n} standard normals; step
#'   \code{k} (from 0) draws streams \code{2k} for x and \code{2k + 1} for y.
#' @param ... Not used.
#' @return \code{pgsig()}: a list with \code{sigmay} and \code{sigmaz}.
#'   \code{briggs()}: a list with \code{flux}, \code{finalrise},
#'   \code{xfinal} and \code{rise} (at each \code{x}). \code{gplume()} and
#'   \code{gspuff()}: a numeric vector of concentrations, one per receptor.
#'   \code{advdif2d()}: a list of class \code{advdif2d} with \code{field},
#'   \code{mass}, \code{cfl}, \code{difn}, \code{stabn}, \code{stable},
#'   \code{dt}, \code{substeps}, the axes \code{x} and \code{y}, \code{dx},
#'   \code{dy}, \code{units} and \code{title}. \code{lagrprt()}: a list of
#'   class \code{lagrprt} with \code{x}, \code{y} (final particle
#'   positions), \code{mean_x}, \code{meany}, \code{grid} (as given) and,
#'   with a grid, \code{concentration} (an \code{nx} by \code{ny} matrix of
#'   particle density per unit area, summing to the share of particles
#'   inside the grid divided by the cell area).
#' @references Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the
#'   air: the contemporaneous relationship between air pollution and crime.
#'   Journal of the Association of Environmental and Resource Economists
#'   7(3), 555-585. \doi{10.1086/707127}. Herrnstadt, E., Heyes, A.,
#'   Muehlegger, E. and Saberian, S. (2019). Air pollution as a cause of
#'   violent crime. Working paper. Briggs, G. A. (1973). Diffusion
#'   estimation for small
#'   emissions. ATDL Contribution 79. Briggs, G. A. (1975). Plume rise
#'   predictions. American Meteorological Society. Turner, D. B. (1970).
#'   Workbook of Atmospheric Dispersion Estimates. EPA AP-26. Seinfeld, J.
#'   H. and Pandis, S. N. (2016). Atmospheric Chemistry and Physics, 3rd ed.
#'   Thomson, D. J. (1987). Criteria for the selection of stochastic models
#'   of particle trajectories in turbulent flows. Journal of Fluid Mechanics
#'   180, 529-556. Steele, G. L., Lea, D. and Flood, C. H. (2014). Fast
#'   splittable pseudorandom number generators. OOPSLA 2014, 453-472.
#'   \doi{10.1145/2660193.2660195}
#' @examples
#' # plume widths 1 km downwind: neutral, then very stable, urban
#' pgsig(1000, "D")
#' pgsig(c(500, 1000, 2000), "F", set = "urban")
#'
#' # rise of a hot plume, and the effective height it gives a 50 m stack
#' rise <- briggs(c(100, 1000), u = 4, diam = 1, exitvel = 12,
#'                stktmp = 420, ambtmp = 285)
#' rise$rise
#' heff <- 50 + rise$finalrise
#'
#' # ground-level concentration (g/m^3) from 100 g/s along the centreline
#' rec <- cbind(x = c(500, 1000, 2000, 8000), y = 0, z = 0)
#' gplume(100, u = 4, h = heff, rcpt = rec)
#'
#' # under a 120 m mixing lid the far receptor sees a trapped plume
#' gplume(100, u = 4, h = heff, rcpt = rec, mixh = 120)
#'
#' # a 1 kg puff, 10 minutes after release, at the puff centre
#' gspuff(1000, u = 3, h = 10, rcpt = rbind(c(1800, 0, 0)), t = 600)
#'
#' # a small grid: a blob moved by wind and spread by diffusion
#' c0 <- matrix(0, 15, 15)
#' c0[4, 8] <- 1
#' f <- advdif2d(c0, u = 1, v = 0, kx = 0.5, ky = 0.5,
#'               dx = 1, dy = 1, dt = 0.2, nsteps = 20)
#' f$stable
#' round(f$mass, 3)
#'
#' # 2,000 simulated particles, counted on a grid
#' p <- lagrprt(2000, x0 = 0, y0 = 0, u = 1, v = 0, kx = 1, ky = 1,
#'              dt = 1, nsteps = 10, seed = 1,
#'              grid = c(-10, 30, -20, 20, 8, 8))
#' c(p$mean_x, p$meany)
#'
#' # inputs that cannot be right are refused
#' try(pgsig(100, "G"))
#' try(briggs(100, u = 4, diam = 1, exitvel = 12, stktmp = 155, ambtmp = 20))
#' @export
pgsig <- function(x, stab = "D", set = "rural") {
  f <- "pgsig"
  cls <- .adcls(stab, f)
  urb <- .adset(set, f)
  x <- .adnum(x, "x", f, min = 0, scalar = FALSE)
  r <- .Call(C_clab_r_pgsig, x, cls, urb)
  list(sigmay = r[[1L]], sigmaz = r[[2L]])
}

# The checks of the dispersion inputs, in their words.
.adnum <- function(
  x,
  what,
  fn,
  min = -Inf,
  positive = FALSE,
  scalar = TRUE,
  whole = FALSE
) {
  if (is.logical(x) && !all(is.na(x))) {
    .ckstop(fn, "`%s` must be numeric, not logical", what)
  }
  if (is.null(x)) .ckstop(fn, "`%s` must be numeric, got NULL", what)
  if (is.factor(x) || is.list(x) || !(is.numeric(x) || is.logical(x))) {
    .ckstop(fn, "`%s` must be numeric, got %s", what, class(x)[1L])
  }
  v <- as.numeric(x)
  if (scalar && length(v) != 1L) {
    .ckstop(fn, "`%s` must be a single number", what)
  }
  if (!length(v) || anyNA(v) || !all(is.finite(v))) {
    .ckstop(fn, "`%s` must be finite, with no missing values", what)
  }
  if (positive && any(v <= 0)) .ckstop(fn, "`%s` must be positive", what)
  if (any(v < min)) .ckstop(fn, "`%s` must be at least %s", what, format(min))
  if (any(abs(v) > 1e300)) .ckstop(fn, "`%s` must be below 1e300 in size", what)
  if (whole && any(v != round(v))) {
    .ckstop(fn, "`%s` must be a whole number", what)
  }
  v
}

.adcls <- function(stab, fn) {
  cls <- toupper(as.character(stab))
  if (length(cls) != 1L || is.na(cls) || !cls %in% LETTERS[1:6]) {
    .ckstop(fn, "`stab` must be one of the Pasquill classes \"A\" to \"F\"")
  }
  match(cls, LETTERS[1:6]) - 1L
}

.adset <- function(set, fn) {
  if (!is.character(set) || length(set) != 1L ||
        !set %in% c("rural", "urban")) {
    .ckstop(fn, "`set` must be \"rural\" or \"urban\"")
  }
  as.integer(set == "urban")
}

.adimg <- function(nimg, fn) {
  nimg <- .adnum(nimg, "nimg", fn, min = 0, whole = TRUE)
  if (nimg > 1000) {
    .ckstop(
      fn,
      paste0(
        "`nimg` must be at most 1,000: the image terms decay ",
        "geometrically and more add only time"
      )
    )
  }
  as.integer(nimg)
}

.adlid <- function(mixh, h, fn) {
  if (is.null(mixh)) return(NA_real_)
  lid <- .adnum(mixh, "mixh", fn, positive = TRUE)
  if (h > lid) {
    warning(
      fn, ": release height ", format(h), " is above the mixing height ",
      format(lid), ": the plume is above the lid and the reflection terms ",
      "do not describe it.",
      call. = FALSE
    )
  }
  lid
}

.adrec <- function(r, fn) {
  bad <- function() {
    .ckstop(fn, "`rcpt` must give x, y, z for at least one receptor")
  }
  m <- if (is.list(r) && !is.data.frame(r)) {
    if (!length(r) || any(lengths(r) != 3L)) bad()
    do.call(rbind, lapply(r, as.numeric))
  } else if (is.matrix(r) || is.data.frame(r)) {
    if (ncol(r) != 3L) bad()
    matrix(as.numeric(as.matrix(r)), ncol = 3)
  } else {
    if (!length(r) || length(r) %% 3L != 0L) bad()
    matrix(as.numeric(r), ncol = 3, byrow = TRUE)
  }
  if (!nrow(m) || anyNA(m)) bad()
  m
}

.adhts <- function(R, lid, fn) {
  z <- R[, 3L]
  if (any(z < 0)) {
    warning(sprintf(
      paste(
        "%s: %d receptor(s) below ground level (z < 0): the plume is",
        "reflected at the ground, so these values are mirror images."
      ),
      fn, sum(z < 0)
    ), call. = FALSE)
  }
  if (is.finite(lid) && any(z > lid)) {
    warning(sprintf(
      paste(
        "%s: %d receptor(s) above the mixing height %s: the lid reflects",
        "the plume, so these values are not concentrations in the free",
        "atmosphere."
      ),
      fn, sum(z > lid), format(lid)
    ), call. = FALSE)
  }
}

#' @rdname pgsig
#' @export
briggs <- function(
  x,
  u,
  diam,
  exitvel,
  stktmp,
  ambtmp,
  stab = "D",
  dthdz = NULL,
  g = 9.80616
) {
  f <- "briggs"
  cls <- .adcls(stab, f)
  if (!is.numeric(u) || length(u) != 1L || !is.finite(u) || u <= 0) {
    .ckstop(f, "`u` must be a single positive wind speed")
  }
  x <- .adnum(x, "x", f, min = 0, scalar = FALSE)
  diam <- .adnum(diam, "diam", f, positive = TRUE)
  exitvel <- .adnum(exitvel, "exitvel", f, positive = TRUE)
  stktmp <- .adnum(stktmp, "stktmp", f, positive = TRUE)
  ambtmp <- .adnum(ambtmp, "ambtmp", f, positive = TRUE)
  if (stktmp < 150 || ambtmp < 150) {
    .ckstop(
      f,
      paste(
        "temperatures are in kelvin: stktmp %s and ambtmp %s look like",
        "degrees Celsius (ambient air is about 290 K)"
      ),
      format(stktmp),
      format(ambtmp)
    )
  }
  g <- .adnum(g, "g", f, positive = TRUE)
  if (stktmp < ambtmp - 0.5) {
    .ckstop(
      f,
      paste(
        "`stktmp` is below `ambtmp`: Briggs' formulas describe a buoyant",
        "plume, so the stack must be as hot as the air or hotter (equal",
        "gives no rise)"
      )
    )
  }
  # within half a kelvin of the air there is no buoyancy to speak of
  stktmp <- max(stktmp, ambtmp)
  dz <- if (is.null(dthdz)) {
    NA_real_
  } else {
    .adnum(dthdz, "dthdz", f, positive = TRUE)
  }
  r <- .Call(C_clab_r_briggs, x, u, diam, exitvel, stktmp, ambtmp, cls, dz, g)
  list(
    flux = r[[1L]][1L], finalrise = r[[1L]][2L], xfinal = r[[1L]][3L],
    rise = r[[2L]]
  )
}

#' @rdname pgsig
#' @export
gplume <- function(
  q,
  u,
  h,
  rcpt,
  stab = "D",
  set = "rural",
  mixh = NULL,
  nimg = 3
) {
  f <- "gplume"
  q <- .adnum(q, "q", f, min = 0)
  u <- .adnum(u, "u", f, positive = TRUE)
  h <- .adnum(h, "h", f, min = 0)
  nimg <- .adimg(nimg, f)
  lid <- .adlid(mixh, h, f)
  R <- .adrec(rcpt, f)
  .adhts(R, lid, f)
  cls <- .adcls(stab, f)
  urb <- .adset(set, f)
  .Call(C_clab_r_gplume, q, u, h, R, cls, urb, lid, nimg)
}

#' @rdname pgsig
#' @export
gspuff <- function(
  mass,
  u,
  h,
  rcpt,
  t,
  stab = "D",
  set = "rural",
  mixh = NULL,
  nimg = 3
) {
  f <- "gspuff"
  mass <- .adnum(mass, "mass", f, min = 0)
  u <- .adnum(u, "u", f, positive = TRUE)
  h <- .adnum(h, "h", f, min = 0)
  t <- .adnum(t, "t", f, positive = TRUE)
  nimg <- .adimg(nimg, f)
  lid <- .adlid(mixh, h, f)
  R <- .adrec(rcpt, f)
  .adhts(R, lid, f)
  cls <- .adcls(stab, f)
  urb <- .adset(set, f)
  .Call(C_clab_r_gspuff, mass, u, h, R, t, cls, urb, lid, nimg)
}

#' @rdname pgsig
#' @section Boundary: the outermost ring of cells is held at zero (an open
#'   boundary): mass that reaches it leaves the domain and \code{mass} falls
#'   accordingly. Pad \code{c0} with empty cells so the plume of interest
#'   stays inside the ring.
#' @export
advdif2d <- function(
  c0,
  u,
  v,
  kx,
  ky,
  dx,
  dy,
  dt,
  nsteps,
  source = NULL,
  maxwork = 1e8
) {
  f <- "advdif2d"
  if (!is.matrix(c0) || !is.numeric(c0)) {
    .ckstop(
      f,
      paste(
        "`c0` must be a numeric matrix (one row per x cell, one column per",
        "y cell)"
      )
    )
  }
  cm <- unname(c0) * 1
  if (anyNA(cm) || !all(is.finite(cm))) {
    .ckstop(f, "`c0` must be a finite numeric matrix")
  }
  if (any(cm < 0)) {
    .ckstop(f, "`c0` must be non-negative: it is a concentration field")
  }
  if (!is.null(source)) {
    if (!identical(dim(as.matrix(source)), dim(cm))) {
      .ckstop(f, "`source` must have the dimensions of `c0`")
    }
    if (!is.numeric(source) || anyNA(source) || !all(is.finite(source))) {
      .ckstop(f, "`source` must be a finite numeric matrix")
    }
    if (any(source < 0)) {
      .ckstop(
        f,
        "`source` must be non-negative: it is an emission rate, like `c0`"
      )
    }
    source <- unname(as.matrix(source)) * 1
  }
  nx <- nrow(cm)
  ny <- ncol(cm)
  u <- .adnum(u, "u", f)
  v <- .adnum(v, "v", f)
  kx <- .adnum(kx, "kx", f, min = 0)
  ky <- .adnum(ky, "ky", f, min = 0)
  dx <- .adnum(dx, "dx", f, positive = TRUE)
  dy <- .adnum(dy, "dy", f, positive = TRUE)
  dt <- .adnum(dt, "dt", f, positive = TRUE)
  nsteps <- .adnum(nsteps, "nsteps", f, min = 0, whole = TRUE)
  nb <- .adnb(u, v, kx, ky, dx, dy, dt)
  stab <- nb[[5L]]
  if (!is.finite(stab)) {
    .ckstop(
      f,
      "the stability number is not finite: check u, v, kx, ky, dx, dy and dt"
    )
  }
  sub <- if (stab > 1) ceiling(stab) else 1
  work <- as.numeric(nx) * ny * nsteps * sub
  if (!is.finite(work)) {
    .ckstop(f, "the amount of work is not finite: check nsteps and the grid")
  }
  if (!is.numeric(maxwork) || length(maxwork) != 1L || is.na(maxwork) ||
        maxwork <= 0) {
    .ckstop(f, "`maxwork` must be a single positive number (Inf lifts the cap)")
  }
  if (work > maxwork) {
    .ckstop(
      f,
      paste(
        "this run updates %s cells x steps x sub-steps (%d x %d cells, %s",
        "steps, %s sub-steps per step%s), above maxwork = %s. Pass",
        "maxwork = Inf to run it anyway, or use a coarser grid, fewer steps",
        "or a shorter total time"
      ),
      format(work, big.mark = ",", scientific = FALSE), nx, ny,
      format(nsteps, big.mark = ",", scientific = FALSE),
      format(sub, big.mark = ",", scientific = FALSE),
      if (sub > 1) {
        sprintf(" from a stability number of %s", format(stab, digits = 4))
      } else {
        ""
      },
      format(maxwork, big.mark = ",", scientific = FALSE)
    )
  }
  if (nsteps > .Machine$integer.max) {
    .ckstop(f, "`nsteps` must be at most 2^31 - 1")
  }
  sub <- as.integer(sub)
  if (stab > 1) {
    warning(sprintf(
      paste(
        "advdif2d: explicit scheme outside its stability bound",
        "(|u dt/dx| + |v dt/dy| + 2 (kx dt/dx^2 + ky dt/dy^2) = %.3f > 1):",
        "each step is taken as %d sub-steps of dt/%d; pass a smaller dt to",
        "silence this."
      ),
      stab, sub, sub
    ), call. = FALSE)
    dt <- dt / sub
    nb <- .adnb(u, v, kx, ky, dx, dy, dt)
    stab <- nb[[5L]]
    nsteps <- nsteps * sub
  }
  cm <- .Call(C_clab_r_advdif, cm, u, v, kx, ky, dx, dy, dt, nsteps, source)
  structure(
    list(
      field = cm,
      mass = sum(cm) * dx * dy,
      cfl = abs(nb[[1L]]) + abs(nb[[2L]]),
      difn = 2 * (nb[[3L]] + nb[[4L]]),
      stabn = stab,
      stable = stab <= 1,
      dt = dt,
      substeps = sub,
      x = seq_len(nx) * dx,
      y = seq_len(ny) * dy,
      dx = dx,
      dy = dy,
      units = "concentration",
      title = sprintf(
        "Advection-diffusion after %s steps of dt = %s",
        format(nsteps, big.mark = ",", scientific = FALSE),
        format(dt)
      )
    ),
    class = "advdif2d"
  )
}

# The advection and diffusion numbers and their stability number.
.adnb <- function(u, v, kx, ky, dx, dy, dt) {
  ax <- u * dt / dx
  ay <- v * dt / dy
  dxn <- kx * dt / (dx * dx)
  dyn <- ky * dt / (dy * dy)
  c(ax, ay, dxn, dyn, abs(ax) + abs(ay) + 2 * (dxn + dyn))
}

#' @rdname pgsig
#' @export
lagrprt <- function(
  n,
  x0,
  y0,
  u,
  v,
  kx,
  ky,
  dt,
  nsteps,
  seed = 0,
  grid = NULL,
  normals = NULL,
  maxwork = 1e8
) {
  f <- "lagrprt"
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != round(seed) || seed < -.Machine$integer.max || seed >= 2^53) {
    .ckstop(
      f,
      "`seed` must be a single whole number from -2147483647 to 2^53 - 1"
    )
  }
  n <- .adnum(n, "n", f, positive = TRUE, whole = TRUE)
  if (n > 1e8) .ckstop(f, "`n` must be at most 100 million particles")
  n <- as.integer(n)
  x0 <- .adnum(x0, "x0", f)
  y0 <- .adnum(y0, "y0", f)
  u <- .adnum(u, "u", f)
  v <- .adnum(v, "v", f)
  kx <- .adnum(kx, "kx", f, min = 0)
  ky <- .adnum(ky, "ky", f, min = 0)
  dt <- .adnum(dt, "dt", f, positive = TRUE)
  nsteps <- .adnum(nsteps, "nsteps", f, min = 0, whole = TRUE)
  if (!is.numeric(maxwork) || length(maxwork) != 1L || is.na(maxwork) ||
        maxwork <= 0) {
    .ckstop(f, "`maxwork` must be a single positive number (Inf lifts the cap)")
  }
  if (n * nsteps > maxwork) {
    .ckstop(
      f,
      paste(
        "this run moves %s particles x %s steps, above maxwork = %s; pass",
        "maxwork = Inf to run it anyway"
      ),
      format(n, big.mark = ","),
      format(nsteps, big.mark = ",", scientific = FALSE),
      format(maxwork, big.mark = ",", scientific = FALSE)
    )
  }
  if (nsteps > .Machine$integer.max) {
    .ckstop(f, "`nsteps` must be at most 2^31 - 1")
  }
  if (!is.null(normals) && !is.function(normals)) {
    .ckstop(f, "`normals` must be a function (n, seed, stream) or NULL")
  }
  g <- NULL
  if (!is.null(grid)) {
    if (length(grid) != 6L) {
      .ckstop(f, "`grid` must be c(xmin, xmax, ymin, ymax, nx, ny)")
    }
    g <- .adnum(grid, "grid", f, scalar = FALSE)
    if (g[2L] <= g[1L] || g[4L] <= g[3L] || g[5L] < 1 || g[6L] < 1 ||
          g[5L] != round(g[5L]) || g[6L] != round(g[6L])) {
      .ckstop(
        f,
        "`grid` needs xmin < xmax, ymin < ymax and whole positive nx, ny"
      )
    }
    if (g[5L] * g[6L] > maxwork) {
      .ckstop(
        f,
        paste(
          "a %s x %s grid has %s cells, above maxwork = %s; pass maxwork =",
          "Inf to allocate it anyway"
        ),
        format(g[5L], big.mark = ","), format(g[6L], big.mark = ","),
        format(g[5L] * g[6L], big.mark = ",", scientific = FALSE),
        format(maxwork, big.mark = ",", scientific = FALSE)
      )
    }
  }
  if (is.null(normals)) {
    p <- .Call(C_clab_r_lagr, n, x0, y0, u, v, kx, ky, dt, nsteps, seed, g)
  } else {
    xs <- rep(x0, n)
    ys <- rep(y0, n)
    fx <- sqrt(2 * kx * dt)
    fy <- sqrt(2 * ky * dt)
    draw <- function(s) {
      z <- tryCatch(normals(n, seed, s), error = function(e) {
        .ckstop(
          f,
          "`normals` failed when called as normals(n, seed, stream): %s",
          conditionMessage(e)
        )
      })
      if (!is.numeric(z) || length(z) != n || !all(is.finite(z))) {
        .ckstop(
          f,
          "`normals` must return %d finite numbers; it returned %s",
          n,
          if (is.numeric(z)) {
            sprintf("%d, %d of them finite", length(z), sum(is.finite(z)))
          } else {
            paste("a", class(z)[[1L]])
          }
        )
      }
      as.numeric(z)
    }
    for (k in seq_len(nsteps) - 1) {
      xy <- .Call(
        C_clab_r_lagstep, xs, ys, u, v, fx, fy, dt, draw(2 * k),
        draw(2 * k + 1)
      )
      xs <- xy[[1L]]
      ys <- xy[[2L]]
    }
    p <- c(list(xs, ys), .Call(C_clab_r_laggrid, xs, ys, g))
  }
  out <- list(
    x = p[[1L]], y = p[[2L]], mean_x = p[[3L]][1L], meany = p[[3L]][2L],
    grid = grid
  )
  if (!is.null(grid)) out$concentration <- p[[4L]]
  structure(out, class = c("lagrprt", "list"))
}

#' A Gaussian-plume concentration field on a grid
#'
#' Evaluates \code{\link{gplume}} at every point of an \code{x} by \code{y}
#' ground-level grid (or at height \code{z}) and returns the field as a
#' matrix with the axes. It maps where the plume from a single stack touches
#' the ground, how far downwind the highest concentration lies, and how
#' wide the affected area is.
#'
#' @details
#' Every grid point \eqn{(x_i, y_j, z)} is passed to \code{\link{gplume}}
#' in one call, with the source at the origin and the wind along \eqn{+x};
#' the field is therefore zero (or vanishingly small) upwind and at
#' \eqn{x \le 0}. The concentration is in the mass unit of \code{q} per
#' cubic metre (grams per cubic metre for \code{q} in grams per second).
#' \code{dx} and \code{dy} are the median grid spacings. The cost is one
#' plume evaluation per grid point.
#'
#' Errors: \code{x} or \code{y} not numeric or not strictly increasing;
#' \code{z} negative; and every check of \code{\link{gplume}}.
#'
#' @param q,u,h,stab,set,mixh,nimg As in \code{\link{gplume}}.
#' @param x,y Numeric vectors, strictly increasing: the downwind and
#'   crosswind coordinates of the grid, in metres from the source.
#' @param z Numeric scalar, the receptor height in metres; default \code{0}
#'   (ground level).
#' @param ... Not used.
#' @return A list of class \code{plmfld}: \code{field} (a numeric matrix of
#'   concentrations, rows along \code{x}, columns along \code{y}), \code{x},
#'   \code{y}, \code{dx} and \code{dy} (median spacings, \code{NA} for a
#'   single coordinate), \code{units} (\code{"concentration"}) and
#'   \code{title}.
#' @seealso \code{\link{gplume}}.
#' @examples
#' f <- plmfld(100, 5, 50, x = seq(100, 3000, by = 100),
#'             y = seq(-500, 500, by = 50))
#' dim(f$field)
#' # where the ground-level maximum lies, downwind
#' f$x[which.max(apply(f$field, 1, max))]
#' # the coordinates must increase
#' try(plmfld(100, 5, 50, x = c(300, 100), y = c(-50, 0, 50)))
#' @export
plmfld <- function(
  q,
  u,
  h,
  x,
  y,
  z = 0,
  stab = "D",
  set = "rural",
  mixh = NULL,
  nimg = 3
) {
  f <- "plmfld"
  x <- .adnum(x, "x", f, scalar = FALSE)
  y <- .adnum(y, "y", f, scalar = FALSE)
  if (is.unsorted(x, strictly = TRUE) || is.unsorted(y, strictly = TRUE)) {
    .ckstop(f, "`x` and `y` must be strictly increasing receptor coordinates")
  }
  z <- .adnum(z, "z", f, min = 0)
  rc <- cbind(rep(x, length(y)), rep(y, each = length(x)), z)
  v <- gplume(q, u, h, rc, stab = stab, set = set, mixh = mixh, nimg = nimg)
  structure(
    list(
      field = matrix(v, nrow = length(x), ncol = length(y)),
      x = x,
      y = y,
      dx = if (length(x) > 1L) stats::median(diff(x)) else NA_real_,
      dy = if (length(y) > 1L) stats::median(diff(y)) else NA_real_,
      units = "concentration",
      title = sprintf("Gaussian plume, class %s, h = %s m", stab, format(h))
    ),
    class = "plmfld"
  )
}

#' @rdname pgsig
#' @export
print.advdif2d <- function(x, ...) {
  fl <- x$field
  pk <- which(fl == max(fl), arr.ind = TRUE)[1L, ]
  cat(
    x$title, "\n",
    sprintf("%d x %d cells of %s by %s", nrow(fl), ncol(fl), format(x$dx),
            format(x$dy)),
    if (isTRUE(x$substeps > 1L)) {
      sprintf(", each step taken as %d sub-steps", x$substeps)
    },
    "\n",
    sprintf(
      "total mass %s; peak %s at x = %s, y = %s\n", format(x$mass),
      format(max(fl)), format(x$x[pk[[1L]]]), format(x$y[pk[[2L]]])
    ),
    sprintf(
      "Courant %.3g, diffusion %.3g, stability %.3f (%s)\n", x$cfl, x$difn,
      x$stabn, if (isTRUE(x$stable)) "within the bound of 1" else "above 1"
    ),
    sep = ""
  )
  invisible(x)
}

#' @rdname plmfld
#' @export
print.plmfld <- function(x, ...) {
  fl <- x$field
  pk <- which(fl == max(fl), arr.ind = TRUE)[1L, ]
  cat(
    x$title, "\n",
    sprintf(
      "%d x %d receptors; peak %s at x = %s, y = %s\n", nrow(fl), ncol(fl),
      format(max(fl)), format(x$x[pk[[1L]]]), format(x$y[pk[[2L]]])
    ),
    sep = ""
  )
  invisible(x)
}

#' @rdname pgsig
#' @export
print.lagrprt <- function(x, ...) {
  cat(sprintf(
    "Lagrangian particles: %d, mean position x = %s, y = %s\n",
    length(x$x), format(x$mean_x), format(x$meany)
  ))
  cat(sprintf(
    "x from %s to %s, y from %s to %s\n", format(min(x$x)), format(max(x$x)),
    format(min(x$y)), format(max(x$y))
  ))
  if (is.matrix(x$concentration)) {
    cat(sprintf(
      "grid %d x %d, peak density %s per unit area\n",
      nrow(x$concentration), ncol(x$concentration),
      format(max(x$concentration))
    ))
  }
  invisible(x)
}
