#' Search theory for search and rescue planning
#'
#' @description
#' The quantities a search and rescue planner works with (Koopman's search
#' theory, as used in the U.S. Coast Guard and IAMSAR planning methods):
#' \itemize{
#'   \item \code{ltrlrng()} (lateral range): the probability that a sensor
#'     passing at distance \code{x} detects the object, under the inverse
#'     cube, definite range or M-Beta model with sweep width \code{W};
#'   \item \code{swpwdth()} (sweep width): the area under a sampled lateral
#'     range curve;
#'   \item \code{srchpod()} (search POD): the probability of detection from
#'     coverage, for random search, a definite-range sweep and an inverse
#'     cube parallel sweep;
#'   \item \code{pswppod()} (parallel sweep POD): the POD of a parallel sweep
#'     at track spacing \code{S}, computed for any model, with optional
#'     normal navigation error;
#'   \item \code{srchscs()} (search success): the probability of success,
#'     the cumulative POD of repeated searches and the containment
#'     probabilities updated after an unsuccessful search;
#'   \item \code{datumpmap()} (datum probability map): the probability that
#'     each grid cell contains the object, for a point, line or area datum;
#'   \item \code{optcrcsrch()} (optimal circular search): Koopman's optimal
#'     effort for a circular normal datum;
#'   \item \code{opteffort()} (optimal effort): the allocation of a fixed
#'     effort over cells that maximises the probability of success.
#' }
#'
#' @details
#' Coverage is \code{C = W z / A}: sweep width times track length over the
#' area searched (for a parallel sweep, \code{W / S}). The lateral range
#' curves are: definite range, \code{p(x) = 1} for \code{|x| <= W / 2};
#' M-Beta, \code{p(x) = m} for \code{|x| <= W / (2 m)}; inverse cube,
#' \code{p(x) = 1 - exp(-W^2 / (4 pi x^2))}. Each has area \code{W}, the
#' sweep width. \code{swpwdth()} integrates by the trapezoid rule.
#'
#' With \code{navsd > 0}, \code{ltrlrng()} gives the expected curve seen
#' through a normal cross-track navigation error: the convolution of the
#' curve with that normal (Frost 1996, equation 5-14), exact for the
#' rectangular curves and by the trapezoid rule (641 points over 8 standard
#' deviations either side) for the inverse cube. Its area is still
#' \code{W}; the curve is lower and wider.
#'
#' \code{srchpod()}: random search \code{1 - exp(-C)}; definite range
#' \code{min(C, 1)}; inverse cube parallel sweep
#' \code{erf(sqrt(pi) C / 2)} (Koopman 1956); M-Beta parallel sweep
#' \code{C} for \code{C < M} and \code{1 - (1 - M)^(C / M)} above it
#' (Frost 1996, equation 5-11), which is exact wherever the tracks overlap
#' a whole number of times.
#'
#' \code{pswppod()} integrates over one track spacing (Simpson's rule on
#' \code{npts} points) the probability that at least one of the tracks
#' \code{-ntrk, ..., ntrk} detects the object; for the inverse cube model
#' the tracks beyond \code{ntrk} are added in closed form from
#' \code{sum_k 1 / (x - k S)^2 = (pi / S)^2 / sin^2(pi x / S)}. With
#' \code{navsd > 0} the lateral range curve is averaged over a normal
#' cross-track error (40-point Gauss-Hermite quadrature, nodes by the
#' Golub-Welsch method). Without navigation error the inverse cube result
#' equals \code{srchpod(W / S, "invcube")}.
#'
#' \code{srchscs()}: the probability of success is
#' \code{POS = sum(POC * POD)}; after an unsuccessful search the
#' containment probabilities become \code{POC (1 - POD) / (1 - POS)}
#' (Bayes). With several searches (\code{pod} a matrix, one row per search)
#' the update is applied in turn and the cumulative POD is
#' \code{1 - prod(1 - POD)} per area. When a search is certain to succeed
#' (\code{POS = 1}) the update is undefined: it is returned as \code{NA},
#' with a warning.
#'
#' \code{datumpmap()}: a point datum is a circular normal of standard
#' deviation \code{sigma} about \code{center}; a line datum the same spread
#' along the segment from \code{center} to \code{end} (Gauss-Legendre
#' quadrature on \code{nodes} points); an area datum a uniform distribution
#' over the rectangle with corners \code{center} and \code{end}; a corridor
#' datum normal across the segment from \code{center} to \code{end} and
#' uniform along it, with nothing beyond its ends (Frost 1996, figure 3-8;
#' the segment must run along x or along y); a general datum the cell
#' weights \code{wts} scaled to sum to 1. The share of the probability
#' outside the grid is returned as \code{outside}.
#'
#' \code{optcrcsrch()}: with effort (track length) \code{L}, the optimal
#' search of a circular normal datum covers a disc of radius
#' \code{R = (4 sigma^2 W L / pi)^(1/4)}, and succeeds with probability
#' \code{1 - (1 + u) exp(-u)}, \code{u = R^2 / (2 sigma^2)} (Koopman 1957).
#'
#' \code{opteffort()}: under random search detection the optimal effort in
#' cell \code{i} is \code{z_i = (A_i / W) log(r_i / lambda)} where
#' \code{r_i = POC_i W / A_i} exceeds \code{lambda}, and 0 elsewhere; the
#' Lagrange multiplier \code{lambda} is found by bisection on its logarithm
#' so that the efforts add up to \code{effort} (Stone 1975, chapter 2). It
#' assumes random search detection in every cell (the exponential POD
#' curve, \code{srchpod(model = "random")}), the conservative choice when
#' the tracks are not laid evenly; a plan searched by parallel sweeps does
#' at least as well.
#'
#' All distances (\code{x}, \code{W}, \code{S}, \code{navsd}, the grid
#' edges, \code{center}, \code{end}, \code{sigma}, \code{effort}) are in one
#' unit of the user's choice, and areas in its square.
#'
#' Every input is checked: numbers must be finite and of the right length,
#' sweep widths, spacings and standard deviations positive, probabilities
#' in \[0, 1\] and containment probabilities summing to at most 1; unknown
#' models and datum types are errors.
#'
#' @param x Numeric vector: lateral ranges (distances from the track), any
#'   sign. For \code{swpwdth()}, strictly increasing abscissae, or an
#'   \code{ltrlrng} object (then \code{p} is taken from it).
#' @param W Positive number: the sweep width, in the distance unit of the
#'   grid and the effort (nautical miles, kilometres: one unit throughout).
#'   For \code{srchpod()}, give it with \code{effort} and \code{area} in
#'   place of \code{coverage}; default \code{NULL}. For \code{opteffort()},
#'   \code{poc * W / area} and \code{area / W} must be finite, and positive
#'   where \code{poc} is (a width of 5e-324 makes every density 0).
#' @param model For \code{ltrlrng()} and \code{pswppod()}: \code{"invcube"}
#'   (the default), \code{"definite"} or \code{"mbeta"}. For
#'   \code{srchpod()}: \code{"random"} (the default), \code{"definite"},
#'   \code{"invcube"} or \code{"mbeta"}.
#' @param m Number in (0, 1]: the M-Beta detection probability inside its
#'   range. Default 1.
#' @param p Numeric vector in \[0, 1\]: the lateral range curve at \code{x}.
#' @param coverage Non-negative numeric vector of coverages (sweep width
#'   times track length over the area searched; 1 means the effort equals
#'   the area). Default \code{NULL}: then \code{W}, \code{effort} and
#'   \code{area} give it, \code{W * effort / area}. Give one or the other,
#'   not both.
#' @param S Positive number: the track spacing.
#' @param navsd Non-negative number: the standard deviation of the
#'   cross-track navigation error, in the units of \code{x}. Default 0.
#' @param npts,ntrk Whole numbers: the quadrature points across one
#'   spacing (at least 3, raised to the next odd number; default 2001) and
#'   the tracks on each side (at least 1; default 60).
#' @param poc Numeric vector or matrix in \[0, 1\], summing to at most 1:
#'   the probability that each area contains the object (the planners' own
#'   areas' shares, for example from a consensus of their scenarios). For
#'   \code{srchscs()} and \code{opteffort()} also a probability map from
#'   \code{datumpmap()}: its cells are the areas, in the order of
#'   \code{as.vector(map$poc)}.
#' @param pod Numeric vector in \[0, 1\], one per area, or a matrix with one
#'   row per search and one column per area. When \code{poc} is a map, also
#'   a matrix laid out like the map (one search, cell by cell), a
#'   \code{srchpod()} result with one value per cell, or an
#'   \code{opteffort()} result (its cells' PODs).
#' @param xedg,yedg Strictly increasing numeric vectors: the cell edges of
#'   the grid laid over the search area, east-west and north-south, in the
#'   distance unit (\code{seq(-5, 5, by = 0.5)}: half-unit cells over ten
#'   units). Cells may differ in size.
#' @param datum \code{"point"} (the default), \code{"line"}, \code{"area"},
#'   \code{"corridor"} or \code{"general"}.
#' @param center Two numbers, x then y in the grid's coordinates: the datum
#'   point (for example the last known position), the start of a line datum
#'   or a corner of an area datum. Default \code{c(0, 0)}.
#' @param sigma Positive number: the standard deviation of the circular
#'   normal position error, in the distance unit (a total probable error
#'   \code{E} is \code{1.1774 sigma}). Default 1.
#' @param end Two numbers: the end of a line datum or the opposite corner of
#'   an area datum; \code{NULL} for a point datum.
#' @param nodes Whole number from 2 to 1000: Gauss-Legendre nodes along a
#'   line datum. Default 64.
#' @param effort Non-negative number: the total effort, as track length
#'   (the distance the searchers travel while searching, in the unit of
#'   \code{W}). For \code{srchpod()}, the track length in each area, one
#'   number or one per area (a matrix laid out like the map when \code{area}
#'   is one); default \code{NULL}.
#' @param area Positive numeric vector, one per cell: the cell areas, in the
#'   square of the unit of \code{W}. May be \code{NULL} when \code{poc} is a
#'   \code{datumpmap} (its cells' areas are used). For \code{srchpod()}, the
#'   areas searched, or a map from \code{datumpmap()} whose cell areas are
#'   taken; default \code{NULL}.
#' @param tol Positive number: the bisection tolerance on log lambda.
#' @param wts For \code{datum = "general"}: a non-negative matrix of cell
#'   weights, one row per \code{y} cell and one column per \code{x} cell.
#' @inheritParams freqtbl
#' @return \code{ltrlrng()}: a numeric vector of class \code{ltrlrng} (the
#'   detection probabilities, with the inputs as attributes;
#'   \code{as.numeric()} drops them). \code{srchpod()}: a numeric vector of
#'   class \code{srchpod}. \code{swpwdth()} and \code{pswppod()}: one
#'   number. \code{srchscs()}: a list of class \code{srchscs} with
#'   \code{pos} (the probability of success), \code{posarea} (per area),
#'   \code{pocupd} (the updated containment probabilities), \code{cumpod}
#'   (the cumulative POD per area), \code{poc} and \code{nsrch}, and when
#'   \code{poc} is a map, \code{map}: the map after the unsuccessful search,
#'   ready for \code{opteffort()}. \code{datumpmap()}: a list of class
#'   \code{datumpmap} with \code{poc}, a matrix with one row per \code{y}
#'   cell and one column per \code{x} cell, \code{outside}, and the inputs.
#'   \code{optcrcsrch()}: a list of class \code{optcrcsrch} with
#'   \code{radius}, \code{pos} and \code{cntrdns} (the coverage at the
#'   centre). \code{opteffort()}: a list of class \code{opteffort} with
#'   \code{effort}, \code{coverage} and \code{pod} per cell (a matrix when
#'   \code{poc} is one), \code{pos} and \code{lambda}. Every class has a
#'   \code{print()} method.
#' @references Koopman, B. O. (1946). Search and Screening. OEG Report 56,
#'   Operations Evaluation Group, Office of the Chief of Naval Operations.
#'
#'   Koopman, B. O. (1956). The theory of search. I. Kinematic bases.
#'   Operations Research 4(3), 324-346. \doi{10.1287/opre.4.3.324}
#'
#'   Koopman, B. O. (1956). The theory of search. II. Target detection.
#'   Operations Research 4(5), 503-531. \doi{10.1287/opre.4.5.503}
#'
#'   Koopman, B. O. (1957). The theory of search. III. The optimum
#'   distribution of searching effort. Operations Research 5(5), 613-626.
#'   \doi{10.1287/opre.5.5.613}
#'
#'   Frost, J. R. (1996). The Theory of Search: A Simplified Explanation.
#'   Soza and Company and U.S. Coast Guard Office of Search and Rescue.
#'
#'   Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
#'   rules. Mathematics of Computation 23(106), 221-230.
#'   \doi{10.1090/S0025-5718-69-99647-1}
#'
#'   Stone, L. D. (1975). Theory of Optimal Search. Academic Press.
#' @srrstats {G1.0} Koopman (1946, 1956, 1957), Frost (1996) and Stone
#'   (1975) are cited in @references.
#' @srrstats {G1.3} @details defines coverage, the lateral range curves,
#'   sweep width, POD, POC and POS with their formulas.
#' @srrstats {G2.0} Scalars (W, S, m, navsd, sigma, effort, tol) must have
#'   length one, `center` and `end` two, `p` the length of `x`, `area` and
#'   the columns of `pod` the length of `poc`.
#' @srrstats {G2.13} Missing values are refused, with the argument named
#'   (test-sar.R).
#' @srrstats {G2.16} Inf, -Inf and NaN are refused, with the argument named
#'   (test-sar.R).
#' @srrstats {G5.4} test-sar.R compares every function with values computed
#'   by the reference implementation and with closed forms.
#' @examples
#' # inverse cube detection 0.5, 1 and 2 sweep widths off the track
#' ltrlrng(c(0.5, 1, 2), W = 1)
#'
#' # the area under the curve is the sweep width
#' xs <- seq(-200, 200, by = 0.01)
#' swpwdth(xs, ltrlrng(xs, W = 1, model = "definite"))
#'
#' # POD at coverage 1 under the three models
#' srchpod(1)
#' srchpod(c(0.5, 1, 2), model = "invcube")
#'
#' # a parallel sweep at spacing 1.3 W, without and with navigation error
#' pswppod(1, 1.3)
#' pswppod(1, 1.3, navsd = 0.3, npts = 201)
#'
#' # three areas searched once; the containment after finding nothing
#' srchscs(c(0.5, 0.3, 0.2), c(0.8, 0.5, 0))
#'
#' # a point datum on a 10 x 10 grid, and the best use of 30 units of track
#' m <- datumpmap(seq(-5, 5), seq(-5, 5), sigma = 2)
#' a <- opteffort(m, W = 1, effort = 30)
#' a$pos
#' optcrcsrch(sigma = 2, W = 1, effort = 30)$pos
#'
#' # a simulated search, step by step: a last known position 2 km east and
#' # 3 km north of the base, known to about 1.5 km, on a grid of 1 km cells
#' own <- datumpmap(seq(-3, 7), seq(-2, 8), center = c(2, 3), sigma = 1.5)
#' own
#' # 0.6 km sweep width; 12 km of track walked in each of the four cells
#' # round the position, none elsewhere
#' walked <- matrix(0, nrow(own$poc), ncol(own$poc))
#' walked[5:6, 5:6] <- 12
#' pd <- srchpod(W = 0.6, effort = walked, area = own, model = "invcube")
#' # the chance that the search found the person, and the map if it did not
#' day1 <- srchscs(own, pd)
#' day1
#' # the best use of the next day's 40 km on what is left
#' opteffort(day1$map, W = 0.6, effort = 40)
#'
#' # inputs that cannot be right are refused
#' try(ltrlrng(1, W = -1))
#' try(srchpod(1, W = 0.6))
#' try(srchscs(c(0.8, 0.8), c(0.5, 0.5)))
#' @export
ltrlrng <- function(
  x,
  W,
  model = c("invcube", "definite", "mbeta"),
  m = 1,
  navsd = 0
) {
  f <- "ltrlrng"
  model <- .ckarg(model, c("invcube", "definite", "mbeta"), "model", f)
  x <- as.vector(.cknum(x, "x", f, scalar = FALSE))
  W <- .cknum(W, "W", f, positive = TRUE)
  m <- .cknum(m, "m", f, positive = TRUE, min = 0, max = 1)
  navsd <- .cknum(navsd, "navsd", f, min = 0)
  structure(
    .Call(C_clab_r_sarlr, x, W, .sarlm(model), m, navsd),
    x = x,
    W = W,
    model = model,
    m = m,
    navsd = navsd,
    class = "ltrlrng"
  )
}

#' @rdname ltrlrng
#' @export
swpwdth <- function(x, p = NULL) {
  f <- "swpwdth"
  if (inherits(x, "ltrlrng") && is.null(p)) {
    p <- as.vector(unclass(x))
    x <- attr(x, "x")
  }
  x <- .sarincr(x, "x", f)
  p <- .cknum(p, "p", f, scalar = FALSE, min = 0, max = 1)
  if (length(p) != length(x)) {
    .ckstop(
      f,
      "`p` must have one value per `x` (%d), got %d",
      length(x),
      length(p)
    )
  }
  .Call(C_clab_r_sarsw, x, as.vector(p))
}

#' @rdname ltrlrng
#' @export
srchpod <- function(
  coverage = NULL,
  model = c("random", "definite", "invcube", "mbeta"),
  m = 1,
  W = NULL,
  effort = NULL,
  area = NULL
) {
  f <- "srchpod"
  ch <- c("random", "definite", "invcube", "mbeta")
  model <- .ckarg(model, ch, "model", f)
  parts <- !vapply(list(W, effort, area), is.null, NA)
  if (!is.null(coverage)) {
    if (any(parts)) {
      .ckstop(f, "give either `coverage` or `W`, `effort` and `area`, not both")
    }
    cv <- as.vector(.cknum(coverage, "coverage", f, scalar = FALSE, min = 0))
  } else {
    if (!all(parts)) {
      .ckstop(f, "give `coverage`, or all three of `W`, `effort` and `area`")
    }
    W <- .cknum(W, "W", f, positive = TRUE)
    if (inherits(area, "datumpmap")) {
      area <- outer(diff(area$yedg), diff(area$xedg))
    }
    effort <- as.vector(.cknum(effort, "effort", f, scalar = FALSE, min = 0))
    area <- as.vector(.cknum(area, "area", f, scalar = FALSE, positive = TRUE))
    if (
      length(effort) != length(area) && min(length(effort), length(area)) > 1
    ) {
      .ckstop(
        f,
        "`effort` and `area` must have the same length, or one of them length 1"
      )
    }
    cv <- W * effort / area
  }
  m <- .cknum(m, "m", f, positive = TRUE, min = 0, max = 1)
  structure(
    .Call(C_clab_r_sarpod, cv, match(model, ch) - 1L, m),
    coverage = cv,
    model = model,
    m = m,
    class = "srchpod"
  )
}

#' @rdname ltrlrng
#' @export
pswppod <- function(
  W,
  S,
  model = c("invcube", "definite", "mbeta"),
  m = 1,
  navsd = 0,
  npts = 2001,
  ntrk = 60
) {
  f <- "pswppod"
  model <- .ckarg(model, c("invcube", "definite", "mbeta"), "model", f)
  W <- .cknum(W, "W", f, positive = TRUE)
  S <- .cknum(S, "S", f, positive = TRUE, max = 1e12)
  m <- .cknum(m, "m", f, positive = TRUE, min = 0, max = 1)
  navsd <- .cknum(navsd, "navsd", f, min = 0)
  npts <- .cknum(npts, "npts", f, min = 3, max = 1e5, whole = TRUE)
  ntrk <- .cknum(ntrk, "ntrk", f, min = 1, max = 1e4, whole = TRUE)
  .Call(
    C_clab_r_sarpsw, W, S, .sarlm(model), m, navsd,
    as.integer(npts), as.integer(ntrk)
  )
}

#' @rdname ltrlrng
#' @export
srchscs <- function(poc, pod) {
  f <- "srchscs"
  map <- if (inherits(poc, "datumpmap")) poc
  if (!is.null(map)) {
    poc <- map$poc
  }
  if (inherits(pod, "opteffort")) {
    pod <- pod$pod
  }
  # a POD matrix laid out like the map is one search, cell by cell
  if (!is.null(map) && identical(dim(pod), dim(map$poc))) {
    pod <- as.vector(pod)
  }
  poc <- as.vector(.sarpoc(poc, f))
  pod <- .cknum(pod, "pod", f, scalar = FALSE, min = 0, max = 1)
  D <- if (is.matrix(pod)) pod else matrix(pod, 1)
  if (ncol(D) != length(poc)) {
    .ckstop(
      f,
      "`pod` must have one value (or column) per area of `poc` (%d), got %d",
      length(poc),
      ncol(D)
    )
  }
  r <- .Call(C_clab_r_sarscs, poc, D)
  if (r[[5L]]) {
    warning(
      "srchscs: the search is certain to find the object (POS = 1), ",
      "so the updated containment is undefined and returned as NA.",
      call. = FALSE
    )
  }
  out <- list(
    pos = r[[1L]],
    posarea = r[[2L]],
    pocupd = r[[3L]],
    cumpod = r[[4L]],
    poc = poc,
    nsrch = nrow(D)
  )
  if (!is.null(map)) {
    map$poc[] <- r[[3L]]
    map$outside <- 1 - sum(r[[3L]])
    out$map <- map
  }
  structure(out, class = "srchscs")
}

#' @rdname ltrlrng
#' @export
datumpmap <- function(
  xedg,
  yedg,
  datum = c("point", "line", "area", "corridor", "general"),
  center = c(0, 0),
  sigma = 1,
  end = NULL,
  nodes = 64,
  wts = NULL
) {
  f <- "datumpmap"
  ch <- c("point", "line", "area", "corridor", "general")
  datum <- .ckarg(datum, ch, "datum", f)
  xedg <- .sarincr(xedg, "xedg", f)
  yedg <- .sarincr(yedg, "yedg", f)
  center <- .cknum(center, "center", f, scalar = FALSE)
  if (length(center) != 2L) {
    .ckstop(f, "`center` must be two numbers, x and y")
  }
  sigma <- .cknum(sigma, "sigma", f, positive = TRUE)
  end <- .dmend(datum, end, f)
  nodes <- .cknum(nodes, "nodes", f, min = 2, max = 1000, whole = TRUE)
  .dmshape(datum, center, end, f)
  w <- if (datum == "general") .dmwts(wts, xedg, yedg, f)
  grid <- .Call(
    C_clab_r_sarmap, xedg, yedg, match(datum, ch) - 1L, center, sigma,
    if (is.null(end)) c(0, 0) else end, as.integer(nodes), w
  )
  structure(
    list(
      poc = grid,
      outside = 1 - sum(grid),
      xedg = xedg,
      yedg = yedg,
      datum = datum,
      center = center,
      end = end,
      sigma = sigma
    ),
    class = "datumpmap"
  )
}

#' @rdname ltrlrng
#' @export
optcrcsrch <- function(sigma, W, effort) {
  f <- "optcrcsrch"
  sigma <- .cknum(sigma, "sigma", f, positive = TRUE)
  W <- .cknum(W, "W", f, positive = TRUE)
  effort <- .cknum(effort, "effort", f, min = 0)
  r <- .Call(C_clab_r_sarcrc, sigma, W, effort)
  structure(
    list(
      radius = r[1L],
      pos = r[2L],
      cntrdns = r[3L],
      sigma = sigma,
      W = W,
      total = effort
    ),
    class = "optcrcsrch"
  )
}

#' @rdname ltrlrng
#' @export
opteffort <- function(poc, area = NULL, W, effort, tol = 1e-14) {
  f <- "opteffort"
  edges <- NULL
  if (inherits(poc, "datumpmap")) {
    edges <- list(xedg = poc$xedg, yedg = poc$yedg)
    if (is.null(area)) {
      area <- outer(diff(poc$yedg), diff(poc$xedg))
    }
    poc <- poc$poc
  }
  poc <- .sarpoc(poc, f)
  area <- .oparea(area, poc, f)
  W <- .cknum(W, "W", f, positive = TRUE)
  effort <- .cknum(effort, "effort", f, min = 0)
  tol <- .cknum(tol, "tol", f, positive = TRUE)
  .opdens(poc, area, W, f)
  r <- .Call(
    C_clab_r_saropt, as.vector(poc), as.vector(area), W, effort, tol
  )
  dm <- if (is.null(dim(poc))) dim(area) else dim(poc)
  shp <- function(v) {
    dim(v) <- dm
    v
  }
  structure(
    c(
      list(
        effort = shp(r[[1L]]),
        coverage = shp(r[[2L]]),
        pod = shp(r[[3L]]),
        pos = r[[4L]],
        lambda = r[[5L]],
        poc = poc,
        W = W,
        total = effort
      ),
      edges
    ),
    class = "opteffort"
  )
}

#' @rdname ltrlrng
#' @export
print.ltrlrng <- function(x, ...) {
  d <- data.frame(range = attr(x, "x"), pod = as.vector(unclass(x)))
  cat(
    "Lateral range curve: ", .sarname(attr(x, "model")), " model, sweep ",
    "width W = ", format(attr(x, "W")),
    if (attr(x, "model") == "mbeta") paste0(", m = ", format(attr(x, "m"))),
    if (attr(x, "navsd") > 0) {
      paste0(", navigation error sd ", format(attr(x, "navsd")))
    },
    "\n",
    sep = ""
  )
  print(utils::head(d, 20L), row.names = FALSE, digits = 6L)
  if (nrow(d) > 20L) cat("... and", nrow(d) - 20L, "more\n")
  cat("pod: the chance of detecting the object at that distance from the",
      "track;\nthe sweep width W is the area under this curve.\n")
  invisible(x)
}

#' @rdname ltrlrng
#' @export
print.srchpod <- function(x, ...) {
  d <- data.frame(coverage = attr(x, "coverage"), pod = as.vector(unclass(x)))
  cat(
    "Probability of detection: ", .sarname(attr(x, "model")), " model",
    if (attr(x, "model") == "mbeta") paste0(", M = ", format(attr(x, "m"))),
    "\n",
    sep = ""
  )
  print(utils::head(d, 20L), row.names = FALSE, digits = 6L)
  if (nrow(d) > 20L) cat("... and", nrow(d) - 20L, "more\n")
  invisible(x)
}

#' @rdname ltrlrng
#' @export
print.srchscs <- function(x, ...) {
  d <- data.frame(
    area = seq_along(x$poc),
    poc = x$poc,
    cumpod = x$cumpod,
    posarea = x$posarea,
    pocupd = x$pocupd
  )
  big <- !is.null(x$map) && nrow(d) > 20L
  if (big) {
    d <- d[utils::head(order(d$poc, decreasing = TRUE), 20L), ]
  }
  cat(sprintf(
    "Probability of success %.4g after %d search%s over %d areas\n",
    x$pos, x$nsrch, if (x$nsrch == 1L) "" else "es", length(x$poc)
  ))
  print(d, row.names = FALSE, digits = 4L)
  if (big) {
    cat(sprintf("The %d most likely of %d cells.\n", nrow(d), length(x$poc)))
  }
  if (!is.null(x$map)) {
    cat(sprintf(
      "Map: %d x %d cells (x by y), updated in $map\n",
      ncol(x$map$poc), nrow(x$map$poc)
    ))
  }
  invisible(x)
}

#' @rdname ltrlrng
#' @export
print.datumpmap <- function(x, ...) {
  p <- x$poc
  xc <- x$xedg[-1L] - diff(x$xedg) / 2
  yc <- x$yedg[-1L] - diff(x$yedg) / 2
  o <- utils::head(order(p, decreasing = TRUE), 5L)
  d <- data.frame(x = xc[col(p)[o]], y = yc[row(p)[o]], poc = p[o])
  cat(sprintf(
    paste0(
      "Datum probability map (%s datum), sigma %s: %d x %d cells (x by y);",
      " %.1f%% of the probability lies outside the grid\n"
    ),
    x$datum, format(x$sigma), ncol(p), nrow(p), 100 * x$outside
  ))
  cat("The most likely cells:\n")
  print(d, row.names = FALSE, digits = 4L)
  invisible(x)
}

#' @rdname ltrlrng
#' @export
print.optcrcsrch <- function(x, ...) {
  cat(sprintf(
    paste0(
      "Optimal search of a circular normal datum (sigma %s, W %s, effort %s)",
      ":\nspread the effort evenly over a circle of radius %.4g: the chance",
      " of finding\nthe object is then %.1f%%; the coverage at the centre ",
      "is %.4g.\n"
    ),
    format(x$sigma), format(x$W), format(x$total), x$radius, 100 * x$pos,
    x$cntrdns
  ))
  invisible(x)
}

#' @rdname ltrlrng
#' @export
print.opteffort <- function(x, ...) {
  used <- which(as.vector(x$effort) > 0)
  o <- utils::head(used[order(as.vector(x$poc)[used], decreasing = TRUE)], 20L)
  d <- data.frame(
    cell = o,
    poc = as.vector(x$poc)[o],
    effort = as.vector(x$effort)[o],
    coverage = as.vector(x$coverage)[o],
    pod = as.vector(x$pod)[o]
  )
  cat(sprintf(
    paste0(
      "Optimal allocation of effort: POS %.4g from effort %s over %d of %d",
      " cells (W = %s)\n"
    ),
    x$pos, format(x$total), length(used), length(x$poc), format(x$W)
  ))
  print(d, row.names = FALSE, digits = 4L)
  if (length(used) > 20L) cat("... and", length(used) - 20L, "more cells\n")
  invisible(x)
}

#' The optimal search square for a circular normal datum
#'
#' For a point datum with a circular normal distribution of position error,
#' the square centred on the datum that a given effort, spread evenly,
#' should cover to make the probability of success as large as possible
#' (Frost 1996, chapter 6). With several efforts, each search is optimized
#' in turn against the distribution as earlier unsuccessful searches left
#' it.
#'
#' @details The search factor \code{f} is the half side of the square in
#' units of the total probable error of position \code{E} (the radius of the
#' circle holding half the distribution, \code{E = 1.1774 sigma}). For
#' relative effort \code{z = Z / E^2} the coverage is \code{C = z / (2 f)^2},
#' the containment is the circular normal's mass in the square,
#' \code{(2 Phi(f E / sigma) - 1)^2} for a first search, and the probability
#' of success \code{POC x POD(C)} with the POD curve of \code{model}
#' (\code{\link{srchpod}}). The factor is found on a grid of 241 values from
#' 0.02 to 10 and refined by Brent's method. After a search the density
#' inside its square is multiplied by \code{1 - POD}; every later search is
#' optimized against what is left, and its \code{poc} and \code{pos} are
#' conditional on the earlier searches having failed.
#'
#' @param effort Positive numeric vector: the effort (sweep width times
#'   track length) of each search in turn, in the square of the units of
#'   \code{E}.
#' @param model \code{"invcube"} (the default; Frost's ideal search
#'   conditions), \code{"random"} (poor conditions), \code{"definite"} or
#'   \code{"mbeta"}.
#' @param m Number in (0, 1]: the M-Beta height. Default 1.
#' @param E Positive number: the total probable error of position. Default
#'   1.
#' @param x An \code{optsqr} result.
#' @inheritParams freqtbl
#' @return A data frame of class \code{optsqr} with one row per search:
#'   \code{search}, \code{effort}, \code{releff} (effort over \code{E^2}),
#'   \code{cumeff} (cumulative relative effort), \code{factor}, \code{side}
#'   (\code{2 f E}), \code{poc}, \code{coverage}, \code{pod}, \code{pos} and
#'   \code{cumpos} (the probability that one of the searches so far
#'   succeeded).
#' @references Frost, J. R. (1996). The Theory of Search: A Simplified
#'   Explanation. Soza and Company and U.S. Coast Guard Office of Search and
#'   Rescue. Chapter 6.
#'
#'   Koopman, B. O. (1957). The theory of search. III. The optimum
#'   distribution of searching effort. Operations Research 5(5), 613-626.
#'   \doi{10.1287/opre.5.5.613}
#'
#'   Brent, R. P. (1973). Algorithms for Minimization without Derivatives.
#'   Prentice-Hall. Chapter 5.
#' @srrstats {G5.4a} test-sar.R checks Frost's stated optima: 4 E^2 gives a
#'   factor of 1.02 and POS 0.4578, 2 E^2 gives 0.81 and 0.2875, 8 E^2 gives
#'   1.29 and 0.6588.
#' @examples
#' # Frost's example: an effort of 4 E^2 on a point datum
#' optsqr(4)
#' # poor search conditions: the random search curve
#' optsqr(4, model = "random")
#' # three searches in a row, each optimal for what the last ones left
#' optsqr(c(4, 4, 8), E = 2)
#' try(optsqr(-1))
#' @export
optsqr <- function(
  effort,
  model = c("invcube", "random", "definite", "mbeta"),
  m = 1,
  E = 1
) {
  f <- "optsqr"
  model <- .ckarg(
    model, c("invcube", "random", "definite", "mbeta"), "model", f
  )
  effort <- as.vector(
    .cknum(effort, "effort", f, scalar = FALSE, positive = TRUE)
  )
  m <- .cknum(m, "m", f, positive = TRUE, min = 0, max = 1)
  E <- .cknum(E, "E", f, positive = TRUE)
  z <- effort / E^2
  r <- .Call(
    C_clab_r_sarsq, z,
    match(model, c("random", "definite", "invcube", "mbeta")) - 1L, m
  )
  structure(
    data.frame(
      search = seq_along(z),
      effort = effort,
      releff = z,
      cumeff = r[, 1L],
      factor = r[, 2L],
      side = 2 * r[, 2L] * E,
      poc = r[, 3L],
      coverage = r[, 4L],
      pod = r[, 5L],
      pos = r[, 6L],
      cumpos = r[, 7L]
    ),
    class = c("optsqr", "data.frame")
  )
}

#' @rdname optsqr
#' @export
print.optsqr <- function(x, ...) {
  cat(sprintf(
    paste0(
      "Optimal square searches: %d, total effort %s, chance found by the",
      " end %.4g\n"
    ),
    nrow(x), format(sum(x$effort)), x$cumpos[nrow(x)]
  ))
  print(as.data.frame(unclass(x)), row.names = FALSE, digits = 4L)
  cat(
    "Each search covers the square (of side 'side') that gives it the best",
    "chance of\nsuccess after the searches before it. poc is the chance the",
    "object is in that\nsquare, pod the chance of detecting it there, pos",
    "their product and cumpos\nthe running chance it has been found.\n"
  )
  invisible(x)
}

#' The figures of The Theory of Search
#'
#' Lists the figures of Frost (1996), \dQuote{The Theory of Search: A
#' Simplified Explanation}, with how each is computed from the search
#' theory functions with the book's parameters, and one line of R that
#' computes a version of it from the reader's own numbers.
#'
#' @details The book's 36 figures (3-1 to 3-14, 4-1 to 4-7, 5-1 to 5-9 and
#' 6-1 to 6-6) are listed with how each is drawn: \code{"model"} figures
#' are computed from the functions with the book's parameters and reproduce
#' the values the text states; \code{"schematic"} figures are diagrams with
#' no data (3-6, 4-3, 4-6, 5-2); \code{"approximation"} figures are drawn
#' from the model where the book's plate cannot be recovered exactly: 3-7
#' (the density implied by boxing in two circles is drawn as the
#' \code{"line"} datum, a circular normal spread along the segment), 3-13
#' (the book's surface is illustrative; it is drawn from the cells of
#' figure 3-14), 5-9 (the plate is blank in the published file; it is drawn
#' from the text, an inverse cube curve with a smaller sweep width and
#' flattened by navigation error), 6-6 (the book's random efforts are not
#' given).
#' @param fig A figure number as text, \code{"5-7"}, or the chapter as a
#'   number with \code{num}; \code{NULL} (the default) lists every figure.
#' @param num The figure's number within chapter \code{fig}, or \code{NULL}.
#' @return A data frame of \code{fig}, \code{caption}, \code{how} and
#'   \code{yours} (the line of R behind the figure, with example inputs to
#'   replace; for a schematic, the computation it pictures): every figure,
#'   or the one asked for.
#' @references Frost, J. R. (1996). The Theory of Search: A Simplified
#'   Explanation. Soza and Company and U.S. Coast Guard Office of Search and
#'   Rescue. Chapters 3-6.
#' @examples
#' head(srchfig())
#' srchfig("3-4")
#' srchfig(5, 8)
#' try(srchfig("9-9"))
#' @export
srchfig <- function(fig = NULL, num = NULL) {
  f <- "srchfig"
  tab <- .Call(C_clab_r_sarfig)
  tab <- data.frame(
    fig = tab[[1L]], caption = tab[[2L]], how = tab[[3L]], yours = tab[[4L]]
  )
  if (is.null(fig)) {
    return(tab)
  }
  if (!is.null(num)) {
    fig <- paste0(
      .cknum(fig, "fig", f, whole = TRUE),
      "-",
      .cknum(num, "num", f, whole = TRUE)
    )
  }
  if (!is.character(fig) || length(fig) != 1L || is.na(fig) || !nzchar(fig)) {
    .ckstop(f, "`fig` must be a single non-empty string")
  }
  if (!fig %in% tab$fig) {
    .ckstop(
      f,
      "`fig` must be one of the book's figures, %s; got \"%s\"",
      toString(tab$fig),
      fig
    )
  }
  tab <- tab[tab$fig == fig, , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

.sarlm <- function(model) match(model, c("invcube", "definite", "mbeta")) - 1L

# Grid edges or abscissae: at least two, strictly increasing.
.sarincr <- function(x, what, fn) {
  if (!is.null(dim(x))) {
    .ckstop(fn, "`%s` must be a vector of edges, not a matrix", what)
  }
  v <- .cknum(x, what, fn, scalar = FALSE)
  if (length(v) < 2L || is.unsorted(v, strictly = TRUE)) {
    .ckstop(fn, "`%s` must hold at least two strictly increasing values", what)
  }
  v
}

# Containment probabilities: each in [0, 1], one positive, sum at most 1.
.sarpoc <- function(poc, fn) {
  poc <- .cknum(poc, "poc", fn, scalar = FALSE, min = 0, max = 1)
  if (sum(poc) > 1 + 1e-9) {
    .ckstop(fn, "`poc` must sum to at most 1, got %s", format(sum(poc)))
  }
  if (!any(poc > 0)) {
    .ckstop(
      fn,
      "`poc` must have at least one cell with a positive probability"
    )
  }
  poc
}

.sarname <- function(model) {
  c(
    invcube = "inverse cube",
    definite = "definite range",
    mbeta = "M-Beta",
    random = "random search",
    corridor = "corridor",
    general = "generalized",
    point = "point",
    line = "line",
    area = "area"
  )[[model]]
}

# datumpmap(): `end`, two numbers, for a line, area or corridor datum
.dmend <- function(datum, end, f) {
  if (datum %in% c("line", "area", "corridor") && is.null(end)) {
    .ckstop(f, "a %s datum needs `end`, two numbers", datum)
  }
  # checked whatever the datum: a point datum passes `end` on unread, but
  # the compiled code still reads it as numbers
  if (!is.null(end)) {
    end <- .cknum(end, "end", f, scalar = FALSE)
    if (length(end) != 2L) .ckstop(f, "`end` must be two numbers, x and y")
  }
  end
}

# datumpmap(): a corridor runs along an axis, an area spans both
.dmshape <- function(datum, center, end, f) {
  if (datum == "corridor") {
    if (center[2L] != end[2L] && center[1L] != end[1L]) {
      .ckstop(f, "a corridor datum must run along x or along y")
    }
    if (identical(center, end)) {
      .ckstop(f, "a corridor datum needs `center` and `end` to differ")
    }
  }
  if (datum == "area" && (center[1L] == end[1L] || center[2L] == end[2L])) {
    .ckstop(
      f,
      "an area datum needs `center` and `end` to differ in x and in y"
    )
  }
}

# datumpmap(): the cell weights of a general datum, a y-by-x matrix
.dmwts <- function(wts, xedg, yedg, f) {
  nx <- length(xedg) - 1L
  ny <- length(yedg) - 1L
  if (is.null(wts)) {
    .ckstop(f, "a general datum needs `wts`, a matrix of cell weights")
  }
  w <- .cknum(wts, "wts", f, scalar = FALSE, min = 0)
  if (!identical(dim(w), c(ny, nx))) {
    .ckstop(
      f,
      paste0(
        "`wts` must be a matrix with one row per y cell (%d) and one ",
        "column per x cell (%d)"
      ),
      ny,
      nx
    )
  }
  if (!any(w > 0)) {
    .ckstop(f, "`wts` must have at least one positive weight")
  }
  w
}

# opteffort(): the cells' areas, one positive number per cell of `poc`
# The densities poc * W / area and the ratios area / W the search works
# with: finite, and positive where poc is, or no root can be bracketed.
.opdens <- function(poc, area, W, f) {
  r <- poc * W / area
  bad <- !all(is.finite(r)) || !all(is.finite(area / W))
  if (bad || any(poc > 0 & r == 0)) {
    .ckstop(
      f,
      paste(
        "`W` = %s is out of range for these cells: poc * W / area and",
        "area / W must be finite, and positive where `poc` is"
      ),
      format(W)
    )
  }
}

.oparea <- function(area, poc, f) {
  if (is.null(area)) {
    .ckstop(f, "`area` must be given, one positive number per cell")
  }
  area <- .cknum(area, "area", f, scalar = FALSE, positive = TRUE)
  if (length(area) != length(poc)) {
    .ckstop(
      f,
      "`area` must have one value per cell of `poc` (%d), got %d",
      length(poc),
      length(area)
    )
  }
  area
}
