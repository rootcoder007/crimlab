# Region maps: crosswalks from units (facilities, detachments, courts) to
# the statistical regions that contain them, checked before tables are
# built on top of them.

#' Internal soundness of a region map
#'
#' Checks that a region map (a crosswalk from units, such as correctional
#' facilities or police detachments, to the statistical regions that
#' contain them) lists each unit once, gives every unit a region, and names
#' only regions the reference geography knows. Run it on any crosswalk
#' before tables are built on it.
#'
#' These are the failures a comparison against published output cannot
#' see, because they would be present on both sides: a unit matched into
#' two regions, a unit matched into none, a region code that is a typo or
#' belongs to a neighbouring jurisdiction. None needs the geometry.
#'
#' The checks as computed: "units assigned more than one region" counts
#' rows whose (non-blank) unit repeats an earlier row's, whether or not the
#' region differs, so an exact duplicate row also fails; "units assigned no
#' region" counts rows whose region is `NA` or blank; "region codes not in
#' the reference geography" counts non-blank region codes not in `regs`,
#' compared as text after trimming spaces (so `"3506"` and `3506` match,
#' `"03506"` does not).
#'
#' @param map The region map, a data frame, one row per unit.
#' @param unit,reg The names of its unit and region columns.
#' @param regs `NULL`, or every valid region code (the identifier column of
#'   the population table); with `NULL` the third check is not run.
#' @return A data frame of class `rgnmapintg`, one row per check: `check`,
#'   `observed`, `expected` (always 0) and `pass`.
#' @seealso [rgnmapcmp()], [rgnmapalt()].
#' @examples
#' # a simulated crosswalk of facilities to census divisions
#' cw <- data.frame(inst = c("North Jail", "South Jail", "East Jail"),
#'                  cd = c("3557", "3520", "3506"))
#' rgnmapintg(cw, "inst", "cd", regs = c("3557", "3520", "3506"))
#' cw$cd[3] <- "2406"
#' rgnmapintg(cw, "inst", "cd", regs = c("3557", "3520", "3506"))
#' cw2 <- rbind(cw, data.frame(inst = c("North Jail", "West Jail"),
#'                             cd = c("3560", "")))
#' rgnmapintg(cw2, "inst", "cd")
#' try(rgnmapintg(cw, "institution", "cd"))
#' @export
rgnmapintg <- function(map, unit, reg, regs = NULL) {
  u <- as.character(.rmcol(map, unit, "map", "rgnmapintg"))
  r <- as.character(.rmcol(map, reg, "map", "rgnmapintg"))
  blank <- function(v) is.na(v) | !nzchar(trimws(ifelse(is.na(v), "", v)))
  obs <- c(sum(duplicated(u[!blank(u)])), sum(blank(r)))
  chk <- c("units assigned more than one region", "units assigned no region")
  if (!is.null(regs)) {
    .rmtxt(regs, "regs", "rgnmapintg")
    obs <- c(obs, sum(!blank(r) & !(trimws(r) %in% trimws(as.character(regs)))))
    chk <- c(chk, "region codes not in the reference geography")
  }
  out <- data.frame(check = chk, observed = as.numeric(obs), expected = 0)
  out$pass <- out$observed == out$expected
  structure(out, class = c("rgnmapintg", "data.frame"))
}

#' @export
print.rgnmapintg <- function(x, ...) {
  bad <- sum(!x$pass)
  cat("rgnmapintg:", if (bad) sprintf("%d of %d checks failed", bad, nrow(x))
      else "every check passed", "\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}

#' Compare a recomputed region map against a published one
#'
#' Matches two region maps on the unit and compares the named columns cell
#' by cell, reporting for each column how many cells disagree and the
#' first disagreement. Use it to confirm that a region map rebuilt from a
#' boundary file reproduces the one that was published.
#'
#' Units are matched as text with spaces trimmed; when a unit repeats, its
#' first row is used. Numeric columns are compared with [all.equal()] at
#' its default tolerance (relative \eqn{1.5 \times 10^{-8}}), so a
#' coordinate that survived a round trip through text is not reported as a
#' change; everything else is compared as text after trimming. Two missing
#' values agree; a missing value against a present one does not.
#'
#' What this establishes is that the recomputation reproduced the
#' published assignment, not that the assignment is right: the same method
#' on the same boundary file reproduces a definitional error perfectly.
#' That is what [rgnmapalt()] is for.
#'
#' @param pub,obs The published and the recomputed map (data frames).
#' @param unit The name of the unit column, in both.
#' @param cols The columns to compare; `NULL` compares every column the two
#'   share apart from `unit`.
#' @return A data frame of class `rgnmapcmp`: a first row `"rows"` for units
#'   on one side only, then one row per column, with `column`, `cells`,
#'   `mismatched` and `first` (the first disagreement written out, `""`
#'   when none). Zero `mismatched` throughout is agreement.
#' @seealso [rgnmapalt()], [rgnmapintg()].
#' @examples
#' pub <- data.frame(inst = c("North Jail", "South Jail"),
#'                   cd = c("3557", "3520"))
#' obs <- pub
#' obs$cd[2] <- "3521"
#' rgnmapcmp(pub, obs, "inst")
#' pub$lat <- c(44.39, 43.65)
#' obs2 <- data.frame(inst = c("North Jail", "Hill Jail"),
#'                    cd = c("3557", "3506"), lat = c(44.39 + 1e-12, 45.1))
#' rgnmapcmp(pub, obs2, "inst", cols = c("cd", "lat"))
#' try(rgnmapcmp(pub, obs, "inst", cols = "lat"))
#' @export
rgnmapcmp <- function(pub, obs, unit, cols = NULL) {
  pu <- trimws(as.character(.rmcol(pub, unit, "pub", "rgnmapcmp")))
  ou <- trimws(as.character(.rmcol(obs, unit, "obs", "rgnmapcmp")))
  both <- intersect(names(pub), names(obs))
  .rmtxt(cols, "cols", "rgnmapcmp")
  cols <- if (is.null(cols)) setdiff(both, unit) else as.character(cols)
  miss <- setdiff(cols, both)
  if (length(miss)) {
    stop("rgnmapcmp: both frames need the compared columns; missing: ",
         toString(miss), ".", call. = FALSE)
  }
  onlyp <- setdiff(pu, ou)
  onlyo <- setdiff(ou, pu)
  shared <- intersect(pu, ou)
  ip <- match(shared, pu)
  io <- match(shared, ou)
  first <- if (length(onlyp)) {
    paste0("published only: ", onlyp[1L])
  } else if (length(onlyo)) {
    paste0("recomputed only: ", onlyo[1L])
  } else {
    ""
  }
  out <- data.frame(column = "rows", cells = length(unique(c(pu, ou))),
                    mismatched = length(onlyp) + length(onlyo), first = first,
                    stringsAsFactors = FALSE)
  for (cl in cols) {
    a <- pub[[cl]][ip]
    b <- obs[[cl]][io]
    same <- if (is.numeric(a) && is.numeric(b)) {
      vapply(seq_along(a), function(i) isTRUE(all.equal(a[i], b[i])), NA)
    } else {
      trimws(as.character(a)) == trimws(as.character(b))
    }
    same[is.na(same)] <- is.na(a[is.na(same)]) & is.na(b[is.na(same)])
    bad <- which(!same)
    out <- rbind(out, data.frame(
      column = cl, cells = length(shared), mismatched = length(bad),
      first = if (length(bad)) {
        sprintf("%s: published %s, recomputed %s", shared[bad[1L]],
                a[bad[1L]], b[bad[1L]])
      } else {
        ""
      }
    ))
  }
  structure(out, class = c("rgnmapcmp", "data.frame"))
}

#' @export
print.rgnmapcmp <- function(x, ...) {
  bad <- sum(x$mismatched, na.rm = TRUE)
  cat("rgnmapcmp:", if (bad) "the maps disagree" else "the maps agree", "\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}

#' Check a region map against an independently derived assignment
#'
#' Compares the region each unit was assigned with the region a different
#' method assigns it (a name lookup, a postal geography, an administrative
#' list), and marks the disagreements that are already known and
#' explained. This is the only check here that can catch an error in the
#' method that built the map, because it does not use that method.
#'
#' A second route usually has understood weaknesses: a place name that is
#' a community rather than a municipality, an amalgamated city, a township
#' absorbed into a neighbour. Listing them in `known` by name keeps the
#' check sharp instead of loosening a tolerance until everything passes,
#' and an unexpected disagreement cannot hide inside an allowance. Units
#' and regions are compared as text with spaces trimmed; units the route
#' does not cover, and route entries that are `NA` or blank, are skipped.
#'
#' @param map The region map, a data frame.
#' @inheritParams rgnmapintg
#' @param route The second method's assignment: a vector named by unit, or
#'   a data frame with the columns `unit` and `reg`.
#' @param known The units whose disagreement is expected and documented.
#' @return A data frame of the disagreements, unexplained first: `unit`,
#'   `primary` (the map's region), `second` (the route's) and `known`.
#'   `sum(!x$known)` unexplained disagreements is the result, zero the
#'   pass; `nrow(x)` should equal the number of documented cases, as a
#'   documented case that stops disagreeing means the route changed.
#' @seealso [rgnmapcmp()], [rgnmapintg()].
#' @examples
#' cw <- data.frame(inst = c("North Jail", "South Jail", "Hill Jail"),
#'                  cd = c("3557", "3520", "3506"))
#' route <- c("North Jail" = "3557", "South Jail" = "3521",
#'            "Hill Jail" = "3519")
#' d <- rgnmapalt(cw, "inst", "cd", route, known = "Hill Jail")
#' d
#' sum(!d$known)
#' try(rgnmapalt(cw, "inst", "cd", c("3557", "3520")))
#' @export
rgnmapalt <- function(map, unit, reg, route, known = character()) {
  u <- trimws(as.character(.rmcol(map, unit, "map", "rgnmapalt")))
  r <- trimws(as.character(.rmcol(map, reg, "map", "rgnmapalt")))
  if (is.data.frame(route)) {
    route <- stats::setNames(
      trimws(as.character(.rmcol(route, reg, "route", "rgnmapalt"))),
      trimws(as.character(.rmcol(route, unit, "route", "rgnmapalt")))
    )
  }
  .rmtxt(route, "route", "rgnmapalt")
  .rmtxt(known, "known", "rgnmapalt")
  route <- route[!is.na(route) & nzchar(trimws(as.character(route)))]
  if (is.null(names(route))) {
    stop("rgnmapalt: `route` must be named by unit, or be a data frame.",
         call. = FALSE)
  }
  names(route) <- trimws(names(route))
  known <- trimws(as.character(known))
  i <- which(u %in% names(route))
  second <- trimws(as.character(route[u[i]]))
  bad <- i[second != r[i]]
  out <- data.frame(unit = u[bad], primary = r[bad],
                    second = trimws(as.character(route[u[bad]])),
                    known = u[bad] %in% known)
  out <- out[order(out$known, out$unit, method = "radix"), , drop = FALSE]
  row.names(out) <- NULL
  out
}

#' Checks that `nm` names one column of the data frame `d`.
#' @noRd
.rmcol <- function(d, nm, arg, fn) {
  if (!is.data.frame(d)) {
    stop(fn, ": `", arg, "` must be a data frame.", call. = FALSE)
  }
  if (!is.character(nm) || length(nm) != 1L) {
    .ckstop(fn, "`%s` must be one column name of `%s`",
            deparse(substitute(nm)), arg)
  }
  .cktxt(nm, deparse(substitute(nm)), fn)
  if (!nm %in% names(d)) {
    stop(fn, ": `", arg, "` has no column `", nm, "`.", call. = FALSE)
  }
  .cktxt(d[[nm]], arg, fn)
}

#' A vector of names or codes (NULL allowed): text, numbers or a factor,
#' valid in its encoding; what is not a vector (a function) is refused.
#' @noRd
.rmtxt <- function(v, what, fn) {
  if (!is.null(v) && !is.atomic(v) && !is.factor(v)) {
    .ckstop(fn, "`%s` must be a vector of names or codes", what)
  }
  .cktxt(v, what, fn)
}
