#' Banded categories: bounds, values and what the open band costs
#'
#' Justice statistics are often published in bands ("1-5 days", "16 to
#' 30", "Over 60", "Less than 1 year"). `bandx()` reads each label into a
#' lower and an upper end, places a value in each band (the midpoint, an
#' end, or the geometric mean), expands a banded frequency table into one
#' value per unit, and shows how much a statistic of those values depends
#' on the cap assumed for an open top band.
#'
#' A label is a range (`"5-10"`, `"5 to 10"`, `"5 through 10"`; any dash),
#' a single number, or a number with a word that opens it: `"60+"`, `"60
#' and over"`, `"60 or more"` (from 60), `"over 60"`, `"more than 60"`
#' (from 61 for whole numbers, `intsca`), `"10 or less"` (to 10), `"under
#' 10"`, `"less than 10"`, `"<10"` (to 9 for whole numbers). Thousands
#' separators are dropped. With `clsup = FALSE` a `"60+"` band starts at 61.
#' An open top band ends at `cap`, or at `fac` times its start; an open
#' bottom band starts at `floor`, or at 0.
#'
#' @param x The band labels: text or a factor.
#' @param cnts `NULL`, or the number of units in each band: numbers from 0.
#' @param rule Where in a band its value sits: `"midpoint"`, `"lower"`,
#'   `"upper"` or `"geometric"` (the midpoint where the lower end is 0).
#' @param cap `NULL`, or the end assumed for an open top band.
#' @param fac Without `cap`, the open top band ends at `fac` times its
#'   start: one positive number.
#' @param floor `NULL`, or the start assumed for an open bottom band.
#' @param clsup `TRUE` when `"60+"` includes 60.
#' @param intsca `TRUE` when the bands hold whole numbers, so that `"over
#'   60"` starts at 61 and `"under 10"` ends at 9.
#' @param caps `NULL`, or the caps at which to recompute `stat` (by
#'   default 12 caps from 1.05 to 20 times the open band's start).
#' @param stat `NULL`, or a function of the expanded values giving one
#'   number (with `cnts`, and an open top band).
#' @param drop `TRUE` to leave out the units of bands that cannot be read;
#'   otherwise they are an error.
#' @return An object of class `bandx`, a list: `bands`, a data frame
#'   (`label`, `lower`, `upper`, `openlo`, `openup`, `value`, `lo`, `hi` the
#'   ends used); `vals`, the expanded values (with `cnts`); `sens` (with
#'   `stat`): a data frame of `cap` and `value`, and `span` and `relspan`
#'   (the span over the caps, absolute and over the median absolute value).
#' @examples
#' b <- c("1-5", "6-15", "16-30", "31-60", "Over 60")
#' n <- c(40, 25, 15, 12, 8)
#' r <- bandx(b, n, stat = mean)
#' r$bands
#' r$sens
#' @export
bandx <- function(
  x,
  cnts = NULL,
  rule = "midpoint",
  cap = NULL,
  fac = 2,
  floor = NULL,
  clsup = TRUE,
  intsca = TRUE,
  caps = NULL,
  stat = NULL,
  drop = FALSE
) {
  b <- .bands(x, clsup, intsca)
  v <- .Call(
    C_clab_r_bandx,
    b$lower,
    b$upper,
    b$openlo,
    b$openup,
    rule,
    cap,
    fac,
    floor
  )
  out <- list(
    bands = cbind(b, v),
    vals = NULL,
    sens = NULL,
    span = NULL,
    relspan = NULL
  )
  if (!is.null(cnts)) {
    if (
      !is.numeric(cnts) ||
        length(cnts) != nrow(b) ||
        anyNA(cnts) ||
        any(cnts < 0)
    ) {
      stop(
        "bandx: `cnts` must be numbers from 0, one per band, ",
        nrow(b),
        ".",
        call. = FALSE
      )
    }
    lost <- is.na(v$value) & cnts > 0
    if (any(lost) && !isTRUE(drop)) {
      stop(
        "bandx: these bands cannot be read: ",
        paste(b$label[lost], collapse = ", "),
        " (drop = TRUE leaves out their units).",
        call. = FALSE
      )
    }
    keep <- !is.na(v$value)
    out$vals <- rep(v$value[keep], round(cnts[keep]))
  }
  if (!is.null(stat)) {
    if (!is.function(stat)) {
      stop("bandx: `stat` must be a function of the values.", call. = FALSE)
    }
    if (is.null(cnts) || !any(b$openup)) {
      stop("bandx: `stat` needs `cnts` and an open top band.", call. = FALSE)
    }
    base <- min(b$lower[b$openup])
    if (is.null(caps)) {
      caps <- unique(round(
        base * exp(seq(log(1.05), log(20), length.out = 12L)),
        6L
      ))
    }
    if (!is.numeric(caps) || !length(caps) || anyNA(caps)) {
      stop("bandx: `caps` must be numbers.", call. = FALSE)
    }
    caps <- sort(unique(caps))
    s <- vapply(
      caps,
      function(cp) {
        e <- bandx(x, cnts, rule, cp, fac, floor, clsup, intsca, drop = drop)
        r <- stat(e$vals)
        if (!is.numeric(r) || length(r) != 1L) {
          stop("bandx: `stat` must give one number.", call. = FALSE)
        }
        r
      },
      0
    )
    out$sens <- data.frame(cap = caps, value = s)
    f <- s[is.finite(s)]
    out$span <- if (length(f)) max(f) - min(f) else NA_real_
    out$relspan <- if (length(f) && all(f != 0)) {
      out$span / stats::median(abs(f))
    } else {
      NA_real_
    }
  }
  structure(out, class = "bandx")
}

# The band labels read into ends: a data frame of label, lower, upper,
# openlo and openup (NA ends where a label cannot be read).
.bands <- function(x, clsup, intsca) {
  if (!(is.character(x) || is.factor(x)) || !length(x)) {
    stop("bandx: `x` must be band labels: text or a factor.", call. = FALSE)
  }
  for (a in list(list(clsup, "clsup"), list(intsca, "intsca"))) {
    if (!is.logical(a[[1L]]) || length(a[[1L]]) != 1L || is.na(a[[1L]])) {
      stop("bandx: `", a[[2L]], "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  lab <- as.character(x)
  s <- tolower(trimws(lab))
  # the hyphens and dashes U+2010 to U+2015 and the minus sign U+2212
  dash <- intToUtf8(c(0x2010:0x2015, 0x2212), multiple = TRUE)
  s <- gsub(paste0("[", paste(dash, collapse = ""), "]"), "-", s)
  s <- gsub("\\s+", " ", s)
  s <- gsub("(?<=[0-9]),(?=[0-9]{3}(?![0-9]))", "", s, perl = TRUE)
  n <- length(s)
  lo <- up <- rep(NA_real_, n)
  ol <- ou <- rep(FALSE, n)
  num <- "(-?[0-9]*\\.?[0-9]+)"
  cue <- paste0(
    "\\b(over|above|more|plus|under|less|below|fewer|greater|",
    "than|gt|lt)\\b"
  )
  for (i in seq_len(n)) {
    z <- s[i]
    if (is.na(z) || !nzchar(z)) {
      next
    }
    m <- regmatches(
      z,
      regexec(
        paste0("^([a-z]+ )?", num, " ?(?:to|-|through) ?", num, "([^0-9]*)$"),
        z,
        perl = TRUE
      )
    )[[1L]]
    if (length(m) == 5L) {
      a <- as.numeric(m[3L])
      b <- as.numeric(m[4L])
      if (
        !grepl(cue, m[2L], perl = TRUE) &&
          !grepl(cue, m[5L], perl = TRUE) &&
          !grepl("+", m[5L], fixed = TRUE) &&
          a <= b
      ) {
        lo[i] <- a
        up[i] <- b
      }
      next
    }
    k <- regmatches(z, gregexpr("-?[0-9]*\\.?[0-9]+", z))[[1L]]
    if (length(k) != 1L) {
      next
    }
    v <- as.numeric(k)
    st <- if (intsca && !grepl("[0-9]\\.[0-9]", z)) 1 else 0
    if (grepl("(\\+\\s*$|and over|and above|or more|or over|plus\\s*$)", z)) {
      lo[i] <- if (clsup) v else v + st
      ou[i] <- TRUE
    } else if (grepl("(greater than|more than|over|above|\\bgt\\b)", z)) {
      lo[i] <- v + st
      ou[i] <- TRUE
    } else if (grepl("(or less|or fewer|and under|and below|or under)", z)) {
      up[i] <- v
      ol[i] <- TRUE
    } else if (grepl("(less than|fewer than|under|below|^<|\\blt\\b)", z)) {
      up[i] <- v - st
      ol[i] <- TRUE
    } else if (grepl("^-?[0-9]*\\.?[0-9]+$", z)) {
      lo[i] <- up[i] <- v
    }
  }
  data.frame(label = lab, lower = lo, upper = up, openlo = ol, openup = ou)
}

#' Bounds implied by rounding and suppression in a published table
#'
#' A published count stands for every true count that could have produced
#' it. Rounded to the nearest multiple of `rnd`, a true count lies within
#' `rnd / 2` of it (randomly rounded, within `rnd - 1`); suppressed below a
#' limit `k` (a cell written `"<k"`, or a marker such as `"x"` with
#' `suplim = k`), it lies from 0 to `k - 1`. Given the cells of an earlier
#' period, `chg` holds every change the two columns allow, and an
#' interval that spans 0 cannot tell a rise from a fall. See
#' `proof("pbbnd")`.
#'
#' @param x The published cells: numbers from 0, or text holding numbers
#'   (with thousands separators), `"<k"` and suppression markers.
#' @param rnd `NULL`, or the rounding base: one number above 0.
#' @param kind `"nearest"` or `"random"` rounding.
#' @param suplim `NULL`, or the suppression limit of the cells marked as
#'   suppressed without one.
#' @param supmar The text markers of a suppressed cell.
#' @param prev `NULL`, or the cells of an earlier period, in the form of
#'   `x`.
#' @param pop,prepop `NULL`, or the populations now and before: one number
#'   or one per cell, above 0.
#' @param per The rates are per `per` people: one positive number.
#' @return An object of class `pbbnd`, a list: `bnd`, a data frame (`value`
#'   the number read, `lower`, `upper`, `status` one of `"exact"`,
#'   `"rounded"`, `"suppressed"` and `"missing"`), and `chg`, with `prev`, a
#'   data frame (`chglo`, `chgup` the change's bounds, `pctlo`, `pctup` the
#'   percent change's, and with populations `ratelo`, `rateup`,
#'   `preratelo`, `prerateup`).
#' @examples
#' pbbnd(c("1,240", "35", "<5", "x"), rnd = 5, suplim = 10)
#' pbbnd(c(120, 15), rnd = 5, prev = c(100, 20), pop = 5e4)
#' @export
pbbnd <- function(
  x,
  rnd = NULL,
  kind = "nearest",
  suplim = NULL,
  supmar = c("x", "X", "s", "..", "F", "suppressed", "n/a"),
  prev = NULL,
  pop = NULL,
  prepop = NULL,
  per = 1000
) {
  a <- .cells(x, supmar, "x")
  p <- if (is.null(prev)) NULL else .cells(prev, supmar, "prev")
  r <- .Call(
    C_clab_r_pbbnd,
    a$value,
    a$st,
    a$lim,
    rnd,
    kind,
    suplim,
    p$value,
    p$st,
    p$lim,
    pop,
    prepop,
    per
  )
  structure(r, class = "pbbnd")
}

# Published cells read: the number (NA for none), the status (0 a number,
# 1 suppressed, 2 missing) and a suppressed cell's "<k" limit.
.cells <- function(x, supmar, nm) {
  if (!is.character(supmar)) {
    stop("pbbnd: `supmar` must be text markers.", call. = FALSE)
  }
  n <- length(x)
  if (!n || !(is.numeric(x) || is.character(x) || is.factor(x))) {
    stop("pbbnd: `", nm, "` must be numbers or text cells.", call. = FALSE)
  }
  lim <- rep(NA_real_, n)
  st <- rep(0L, n)
  if (is.numeric(x)) {
    v <- as.numeric(x)
    if (any(is.infinite(v))) {
      stop("pbbnd: `", nm, "` must be finite.", call. = FALSE)
    }
  } else {
    s <- trimws(as.character(x))
    lt <- !is.na(s) & grepl("^<\\s*[0-9]+$", s)
    lim[lt] <- as.numeric(sub("^<\\s*", "", s[lt]))
    st[lt | s %in% supmar] <- 1L
    v <- rep(NA_real_, n)
    ok <- st == 0L & !is.na(s) & nzchar(s)
    v[ok] <- suppressWarnings(as.numeric(gsub(",", "", s[ok], fixed = TRUE)))
    bad <- ok & is.na(v)
    if (any(bad)) {
      stop(
        "pbbnd: cannot read these cells as numbers or suppression marks: ",
        paste(unique(s[bad]), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  st[st == 0L & is.na(v)] <- 2L
  list(value = v, st = st, lim = lim)
}

#' Year-over-year change in a published series
#'
#' The change of each period against the period `lag` before it, within
#' groups (an agency, a region): the difference, the percent change (or,
#' for percentages, the change in points), and for counts the exact
#' interval of the percent change, from the conditional binomial of the
#' count now among both, the exact two-sided p-value of no change (as
#' [stats::binom.test()]), and the p-values adjusted for the number of
#' comparisons (as [stats::p.adjust()]). With `rnd`, the counts are taken
#' as rounded, and the envelope of the change they allow is added, with a
#' percent interval that covers both.
#'
#' Rows that share a group and a period are summed (counts) or averaged.
#' A flag marks a change with no earlier period, an earlier period of 0,
#' or one below `minbase`.
#'
#' @param x The values: counts from 0, or numbers (`units` other than
#'   `"count"`).
#' @param prd The period of each row: numbers (a year, a month index); the
#'   earlier period is `prd - lag`.
#' @param by `NULL`, or the group of each row: labels.
#' @param lag How many periods back: one positive number.
#' @param units `"count"`, `"continuous"` or `"percent"`.
#' @param minbase Changes from an earlier value below this are flagged: one
#'   number from 0.
#' @param level The level of the intervals: one number between 0 and 1.
#' @param rnd `NULL`, or the rounding base of the counts.
#' @param kind `"nearest"` or `"random"` rounding.
#' @param adj The adjustment of the p-values: `"none"`, `"bonferroni"`,
#'   `"holm"`, `"hochberg"`, `"BH"` or `"BY"`.
#' @param alpha The level at which `signif` marks an adjusted p-value.
#' @return A data frame of class `yoypx` with one row per group and period
#'   (in that order): `group` (with `by`), `prd`, `value`, `prev`,
#'   `change`, `pctchg` (`ppchg` for percentages), `flag`; for counts also
#'   `pctlo`, `pctup`, `p`, `padj`, `signif`, and with `rnd` `envlo`,
#'   `envup`, `envplo`, `envpup` (the envelope of the change and of the
#'   percent change) and `cpctlo`, `cpctup`.
#' @examples
#' d <- data.frame(year = rep(2019:2023, 2),
#'                 svc = rep(c("A", "B"), each = 5),
#'                 n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33))
#' yoypx(d$n, d$year, d$svc, rnd = 5)
#' @export
yoypx <- function(
  x,
  prd,
  by = NULL,
  lag = 1,
  units = "count",
  minbase = 0,
  level = 0.95,
  rnd = NULL,
  kind = "nearest",
  adj = "BH",
  alpha = 0.05
) {
  r <- .Call(
    C_clab_r_yoypx,
    x,
    prd,
    by,
    lag,
    units,
    minbase,
    level,
    rnd,
    kind,
    adj,
    alpha
  )
  structure(r, class = c("yoypx", "data.frame"))
}

#' @export
print.bandx <- function(x, ...) {
  cat("bandx:", nrow(x[["bands"]]), "bands\n")
  print(x[["bands"]], row.names = FALSE)
  if (!is.null(x[["sens"]])) {
    cat(sprintf(
      "over the caps tried the statistic spans %.4g (%.3g of it)\n",
      x[["span"]],
      x[["relspan"]]
    ))
    print(x[["sens"]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.pbbnd <- function(x, ...) {
  cat("pbbnd: what the published cells allow\n")
  print(x[["bnd"]], row.names = FALSE)
  if (!is.null(x[["chg"]])) {
    print(x[["chg"]], row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.yoypx <- function(x, ...) {
  cat("yoypx: change on the earlier period\n")
  print(as.data.frame(unclass(x)), row.names = FALSE)
  invisible(x)
}
