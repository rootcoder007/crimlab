#' Downwind exposure: places downwind of sources of pollution
#'
#' Whether each place lies downwind of sources of pollution in each period,
#' and how far off the wind's line it lies: the instruments of the papers
#' that let the wind carry smoke from fires that move from day to day (Bao
#' and Da 2026) or from a fixed set of sources (Ayesh 2023), for the
#' instrumental-variable estimators [ivfe()], [ppmlcf()] and [expgmm()].
#'
#' The sources are rows (`sx`, `sy`) with a period each (`sper`: a day, a
#' plume-day's date, a year) and the wind there (`sdir`, degrees clockwise
#' from north, where the wind blows from with `from = TRUE`: `sdir + 180`
#' is then the direction it blows to). The places are points (`rx`, `ry`)
#' or, for `rule = "wedge"`, polygons whose vertices are the rows of each
#' `rid` in order. Three rules:
#'
#' * `"sector"`: a place is downwind of a source when the bearing from the
#'   source to the place is within `half` degrees of the direction the wind
#'   blows to, and the place within `rad` of the source (a place at the
#'   source counts). This is the point ("centroid") instrument of Bao and
#'   Da (2026, sec. 3.3), each county collapsed to one point and "the
#'   great-circle bearing from the plume anchor to that point" within 30
#'   degrees either side of the mean wind to-direction (`half = 30`, `rad =
#'   Inf`).
#' * `"wedge"`: a place's polygon meets the wedge of half-angle `half` and
#'   radius `rad` about that direction from the source: their baseline
#'   instrument, "a wedge-shaped transport corridor" with a half-angle of
#'   30 degrees and a radius of 1,000 km, a county downwind "if its census
#'   polygon intersects this wedge" (sec. 3.3). Planar coordinates only
#'   (projected, `rad` in their unit). The polygon meets the wedge when a
#'   vertex lies in it, the apex lies in the polygon, or an edge crosses a
#'   side or the arc.
#' * `"angle"`: per period, the angle between the bearing from the mean
#'   point of the sources to the place and the direction of their mean wind
#'   vector, and the mean distance from the place to the sources (those at
#'   its own location left out): Ayesh's (2023, sec. 4) instrument, "the
#'   absolute value of `lat_j - DirT_t`" with `DirT_t` the "average annual
#'   wind direction ... from the twelve rice growing districts", weighted
#'   by "the mean distance from the rice growing districts" (for a source
#'   district "the mean distance from eleven other rice growing districts",
#'   n. 23).
#'
#' For `"sector"` and `"wedge"`, `n` counts the period's sources that reach
#' the place and `any` is 1 when one does; Bao and Da's instrument is the sum
#' of `any` over the days of a county-month (their eq 1), which
#' `stats::aggregate()` gives (see the examples). With several plumes a day,
#' give each its own row with that day's period.
#'
#' The choices the papers leave open, as arguments:
#'
#' * `from`: Bao and Da convert the "from-direction" to a "to-direction"
#'   by adding 180 degrees; Ayesh subtracts the wind direction from the
#'   azimuth without saying which way the direction is recorded. `from =
#'   TRUE` (the default) adds 180 degrees, so that the angle is 0 for a
#'   place straight downwind; `from = FALSE` uses `sdir` as given.
#' * `wrap`: the angle about the circle, at most 180 (the default), or the
#'   plain absolute difference of the two bearings in `[0, 360)`, Ayesh's
#'   formula as written.
#' * `dwt`: Ayesh does not say how the mean distance weights the angle:
#'   `"mult"` (the default) gives `z = angle * dist`, `"div"` `angle /
#'   dist`, `"none"` the angle alone.
#' * `lonlat`: with `TRUE`, `x` and `y` are longitude and latitude in
#'   degrees, bearings are initial great-circle bearings and distances
#'   kilometres (haversine, radius 6371.0088 km); else planar, distances in
#'   the coordinates' unit.
#'
#' The direction of the wind at a source is the user's: Bao and Da take the
#' speed-weighted circular mean of the to-direction over the wind grid
#' cells inside each plume ([windiv()] averages directions the same way),
#' and keep only plumes of medium and heavy density ([smkiv()] reads those
#' classes). As instruments these are valid only if the wind carries the
#' outcome nowhere but through the exposure.
#'
#' @param sx,sy The coordinates of the sources.
#' @param sper `NULL` (one period) or the period of each source (labels).
#' @param sdir The wind direction at each source, degrees clockwise from
#'   north.
#' @param rx,ry The coordinates of the places: points, or the vertices of
#'   their polygons for `rule = "wedge"`.
#' @param rid `NULL` (each row a place) or the place of each row of `rx`
#'   (labels); for `rule = "wedge"`, the vertices of a place in order.
#' @param rule `"sector"`, `"wedge"` or `"angle"`.
#' @param half The half-angle of the sector, degrees from 0 to 180.
#' @param rad The radius of the sector: a positive number or `Inf`.
#' @param from `TRUE` when `sdir` is where the wind blows from.
#' @param lonlat `TRUE` for longitude and latitude in degrees.
#' @param wrap `TRUE` for the angle about the circle (`rule = "angle"`).
#' @param dwt `"mult"`, `"div"` or `"none"`: the distance weight of the
#'   angle (`rule = "angle"`).
#' @return A data frame, one row per place and period (places fastest):
#'   `id`, `per`, and `n` and `any` for `"sector"` and `"wedge"`, or
#'   `angle`, `dist` and `z` for `"angle"` (`NA` where a period has no
#'   source at a distance).
#' @references
#'   Bao, J. and Da, Y. (2026). Wildfire smoke and crime: evidence from
#'   California. Working paper, Renmin University of China.
#'
#'   Ayesh, A. (2023). Burned agricultural biomass, air pollution and
#'   crime. *Journal of Environmental Economics and Management*, 122,
#'   102887. \doi{10.1016/j.jeem.2023.102887}
#' @seealso [windiv()], [smkiv()], [updwn()], [ivfe()]
#' @examples
#' # two plumes on day 1, one on day 2; four county points (planar, km)
#' s <- data.frame(x = c(0, 50, 0), y = c(0, 0, 0), day = c(1, 1, 2),
#'                 from = c(180, 270, 0))
#' p <- data.frame(x = c(0, 100, 0, -80), y = c(100, 10, -60, 0))
#' r <- dwnexp(s$x, s$y, s$day, s$from, p$x, p$y)
#' r
#' # the county-month count of downwind days (Bao and Da 2026, eq 1)
#' r$month <- 1
#' aggregate(any ~ id + month, r, sum)
#'
#' # Ayesh's instrument: the angle off the mean wind of the source
#' # districts, times the mean distance to them
#' dwnexp(c(0, 10), c(0, 0), c(2020, 2020), c(240, 250), c(30, -30, 5),
#'        c(20, -20, 0), rule = "angle")
#'
#' # the wedge against a square county
#' sq <- data.frame(id = 1, x = c(-10, 10, 10, -10), y = c(40, 40, 60, 60))
#' dwnexp(0, 0, NULL, 180, sq$x, sq$y, sq$id, rule = "wedge", rad = 100)
#' @export
dwnexp <- function(sx, sy, sper, sdir, rx, ry, rid = NULL, rule = "sector",
                   half = 30, rad = Inf, from = TRUE, lonlat = FALSE,
                   wrap = TRUE, dwt = "mult") {
  f <- "dwnexp"
  a <- .dwchk(sx, sy, sper, sdir, rx, ry, rid, rule, half, rad, f)
  from <- .ckflag(from, "from", f)
  lonlat <- .ckflag(lonlat, "lonlat", f)
  wrap <- .ckflag(wrap, "wrap", f)
  dwt <- .ckarg(dwt, c("mult", "div", "none"), "dwt", f)
  if (a$rule == 1L && lonlat) {
    .ckstop(f, "`rule = \"wedge\"` needs planar coordinates (`lonlat = FALSE`)")
  }
  to <- if (from) (as.numeric(sdir) + 180) %% 360 else as.numeric(sdir)
  k <- .Call(C_clab_r_dwnexp, a$rule, as.numeric(sx), as.numeric(sy),
             a$sp, to, a$rx, a$ry, a$ro, length(a$ids), length(a$per),
             as.numeric(half), as.numeric(rad), lonlat, wrap)
  r <- data.frame(id = rep(a$ids, length(a$per)),
                  per = rep(a$per, each = length(a$ids)))
  if (a$rule == 2L) {
    r$angle <- as.vector(k[[1L]])
    r$dist <- as.vector(k[[2L]])
    r$angle[is.na(r$angle)] <- NA_real_
    r$dist[is.na(r$dist)] <- NA_real_
    r$z <- switch(dwt, mult = r$angle * r$dist, div = r$angle / r$dist,
                  none = r$angle)
  } else {
    r$n <- as.vector(k[[1L]])
    r$any <- as.integer(r$n > 0)
  }
  r
}

#' The tailwind component of the wind along a flight course
#'
#' The part of the wind that blows along a course, positive for a
#' tailwind, from which Hener (2022) builds the instrument of the landing
#' direction at Frankfurt Airport: aircraft land westbound (course 250
#' degrees) unless the tailwind on that course exceeds 2.6 m/s, so the
#' days of the eastbound approach, and which municipalities lie under the
#' noise, follow from the wind.
#'
#' For wind from `dir` (degrees clockwise from north, where the wind
#' blows from) at speed `spd`, the component along `course` is `spd *
#' cos(dir + 180 - course)`: "positive for tailwinds and negative for
#' headwinds in the direction of the West25 landing approach", 0 for "no
#' wind at all or if it is blowing exactly perpendicular to the runway
#' (i.e., from 340 or 160 degrees)" (Hener 2022, sec. 3, for `course =
#' 250`). With `by` (the day of hourly readings), each group gives one
#' value: `agg = "max"`, the paper's "daily maximum partial wind speed
#' ... calculated from hourly wind data", or `"mean"`. The instrument is
#' then 1 for a municipality under the eastbound approach on a day whose
#' tailwind exceeds 2.6 m/s, 1 for one under the westbound approach on a
#' day whose tailwind is at most 2.6, and 0 otherwise (sec. 3; see
#' [ivfe()] for the estimate).
#'
#' @param dir Wind directions in degrees (`NA` allowed).
#' @param spd Wind speeds, one per direction (`NA` allowed).
#' @param course The course flown, degrees clockwise from north.
#' @param by `NULL` or the group of each reading.
#' @param agg `"max"` or `"mean"`: a group's value.
#' @return Without `by`, the tailwind component of each reading; with
#'   `by`, a data frame, one row per group (sorted): `by` and `tw` (`NA`
#'   for a group without a reading).
#' @references
#'   Hener, T. (2022). Noise pollution and violent crime. *Journal of
#'   Public Economics*, 215, 104748. \doi{10.1016/j.jpubeco.2022.104748}
#' @seealso [windiv()], [ivfe()], [dwnexp()]
#' @examples
#' # wind from the east (70 degrees) at 4 m/s is a 4 m/s tailwind westbound
#' tailwd(c(70, 340, 250), c(4, 4, 4))
#' # hourly readings to the daily maximum, then Hener's instrument for a
#' # municipality under the eastbound approach
#' tw <- tailwd(c(70, 90, 250, 260), c(3, 2, 5, 1), by = c(1, 1, 2, 2))
#' as.integer(tw$tw > 2.6)
#' @export
tailwd <- function(dir, spd, course = 250, by = NULL, agg = "max") {
  f <- "tailwd"
  if (!is.numeric(dir) || !is.numeric(spd) || length(dir) != length(spd) ||
        !length(dir)) {
    .ckstop(f, "`dir` and `spd` must be numeric vectors of one length")
  }
  course <- .cknum(course, "course", f)
  agg <- .ckarg(agg, c("max", "mean"), "agg", f)
  g <- .twgrp(by, length(dir), f)
  v <- .Call(C_clab_r_tailwd, as.numeric(dir), as.numeric(spd), course,
             if (is.null(g)) NULL else match(by, g) - 1L, length(g),
             match(agg, c("max", "mean")) - 1L)
  v[is.na(v)] <- NA_real_
  if (is.null(g)) v else data.frame(by = g, tw = v)
}

# The checks of dwnexp() in its order; the codes it passes on: rule
# (0-2), the period codes from 0 and their labels, the places' labels and,
# for the wedge, their vertices in order with the offsets.
.dwchk <- function(sx, sy, sper, sdir, rx, ry, rid, rule, half, rad, f) {
  rule <- match(.ckarg(rule, c("sector", "wedge", "angle"), "rule", f),
                c("sector", "wedge", "angle")) - 1L
  .dwlen(list(sx = sx, sy = sy, sdir = sdir), sper, "sper",
         "`sx`, `sy` and `sdir` must have one length",
         "the period of each source", f)
  .dwlen(list(rx = rx, ry = ry), rid, "rid",
         "`rx` and `ry` must have one length",
         "the place of each row of `rx`", f)
  .cknum(half, "half", f, min = 0, max = 180)
  if (!is.numeric(rad) || length(rad) != 1L || is.na(rad) || rad <= 0) {
    .ckstop(f, "`rad` must be a positive number or Inf")
  }
  per <- if (is.null(sper)) 1L else sort(unique(sper))
  sp <- if (is.null(sper)) integer(length(sx)) else match(sper, per) - 1L
  .dwpl(rx, ry, rid, rule, f, list(rule = rule, per = per, sp = sp))
}

# Numeric vectors of one length (named in v), and labels lab of each of
# their rows (NULL passes).
.dwlen <- function(v, lab, lnm, msg, what, f) {
  for (nm in names(v)) .cknum(v[[nm]], nm, f, scalar = FALSE)
  if (length(unique(lengths(v))) != 1L) .ckstop(f, msg)
  if (!is.null(lab) && (length(lab) != length(v[[1L]]) || anyNA(lab))) {
    .ckstop(f, "`%s` must give %s, with no NA", lnm, what)
  }
  invisible(NULL)
}

# The places of dwnexp(): labels, and for the wedge the vertices grouped
# by place in the order of first appearance, with offsets from 0.
.dwpl <- function(rx, ry, rid, rule, f, a) {
  if (rule != 1L) {
    a$ids <- if (is.null(rid)) seq_along(rx) else rid
    if (anyDuplicated(a$ids)) .ckstop(f, "`rid` must name each place once")
    a$rx <- as.numeric(rx)
    a$ry <- as.numeric(ry)
    a$ro <- NULL
    return(a)
  }
  if (is.null(rid)) .ckstop(f, "`rule = \"wedge\"` needs `rid`, the polygons")
  ids <- unique(rid)
  o <- order(match(rid, ids))
  cnt <- tabulate(match(rid, ids), length(ids))
  if (any(cnt < 3L)) {
    .ckstop(f, "each polygon in `rid` needs three vertices or more")
  }
  a$ids <- ids
  a$rx <- as.numeric(rx)[o]
  a$ry <- as.numeric(ry)[o]
  a$ro <- as.integer(c(0L, cumsum(cnt)))
  a
}

# The groups of tailwd(): the sorted values of `by`, or NULL without one.
.twgrp <- function(by, n, f) {
  if (is.null(by)) return(NULL)
  if (length(by) != n || anyNA(by)) {
    .ckstop(f, "`by` must give the group of each reading, with no NA")
  }
  sort(unique(by))
}
