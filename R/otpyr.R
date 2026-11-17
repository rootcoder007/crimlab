#' Person-year table of a restrictive confinement release
#'
#' Reduces a placement-level release of Ontario's restrictive confinement
#' records (table a01 of the open data, one row per person, fiscal year,
#' pair of regions and alert configuration, with a count of placements) to
#' one row per person and fiscal year, carrying placement volatility
#' ([plcvol()]), alert complexity ([alrtcx()]) and the covariates of
#' Ruhela (2026): gender, age group, fiscal year and regional cluster (the
#' regions the person-year's records name, sorted and joined by `" + "`).
#'
#' `x` is the table, or a key of the open data service read by [odtab()].
#' Columns are found by `cols`, whose defaults are the names the release
#' has always used, matched without regard to case or underscores, so the
#' publisher's `end_fiscal_year` finds `EndFiscalYear`. The fiscal years
#' are those in the data; `years` keeps some of them. Identifiers are
#' reassigned every fiscal year in the release, so the person-year is the
#' unit, never the person.
#'
#' Rows missing any column used are dropped and counted (attribute
#' `drop`); an alert flag that is neither yes nor no (nor 1 or 0) is a
#' missing one. Alert complexity counts the six configurations the release
#' shows (`keep`), as the published analysis does.
#'
#' @param x A data frame, or the key of a table of the open data service.
#' @param cols The column names: a named character vector with any of
#'   `year`, `id`, `ra`, `rb`, `gender`, `age`, `mh`, `sr`, `sw`, `np`
#'   (the placements of the row), overriding the defaults.
#' @param years `NULL` for every fiscal year in `x`, or the years to keep.
#' @param ord The order of the records for [plcvol()].
#' @param keep The alert configurations counted by [alrtcx()].
#' @param base The data service, for a key (see [odtab()]).
#' @return A data frame of class `otpyr`, one row per person-year: `id`,
#'   `year`, `n` (rows), `np` (placements), `within`, `across`, `vm`, `ac`,
#'   `cx`, `anyal`, `gender`, `age`, `rc`; attributes `years`, `drop`,
#'   `ord`.
#' @inherit plcvol references
#' @examples
#' d <- data.frame(
#'   EndFiscalYear = c(2025, 2025, 2025, 2026, 2026),
#'   UniqueIndividual_ID = c("a", "a", "b", "a", "c"),
#'   Region_AtTimeOfPlacement = c("East", "West", "East", "North", "East"),
#'   Region_MostRecentPlacement = c("East", "West", "East", "East", "East"),
#'   Gender = "Male", Age_Category = "25 to 49",
#'   MentalHealth_Alert = c("Yes", "No", "No", "No", "Yes"),
#'   SuicideRisk_Alert = "No", SuicideWatch_Alert = "No",
#'   Number_Of_Placements = c(3, 1, 2, 4, 1)
#' )
#' otpyr(d)
#' @export
otpyr <- function(
  x,
  cols = NULL,
  years = NULL,
  ord = "region",
  keep = c("100", "010", "110", "011", "111", "000"),
  base = getOption("crimlab.data", "https://data.rmorie.com")
) {
  fn <- "otpyr"
  cn <- .otcn(cols, fn)
  if (is.character(x) && length(x) == 1L) {
    x <- .clas(fn, odtab(x, cols = unname(cn), base = base))
  }
  if (!is.data.frame(x)) {
    .ckstop(fn, "`x` must be a data frame or a table key")
  }
  d <- x[.odcol(unname(cn), names(x), fn)]
  names(d) <- names(cn)
  # a flag that is neither yes nor no is a missing one: the record is
  # dropped here, so that both counts below see the same records
  fl <- .acflags(d[c("mh", "sr", "sw")], fn)
  ok <- stats::complete.cases(d) & rowSums(is.na(fl)) == 0L
  if (!is.null(years)) {
    ok <- ok & d$year %in% years
  }
  if (!any(ok)) {
    .ckstop(fn, "no row has every column present in the years asked")
  }
  d <- d[ok, , drop = FALSE]
  key <- list(id = d$id, year = d$year)
  v <- .clas(fn, plcvol(d$ra, d$rb, key, ord = ord))
  a <- .clas(fn, alrtcx(d[c("mh", "sr", "sw")], key, keep = keep))
  pk <- c("id", "year", "n")
  if (!identical(unclass(v)[pk], unclass(a)[pk])) {
    .ckstop(fn, "the two counts do not cover the same person-years")
  }
  k <- paste(d$id, d$year, sep = "\r")
  first <- !duplicated(k)
  out <- data.frame(
    v[c("id", "year", "n")],
    np = as.vector(tapply(as.numeric(d$np), factor(k, unique(k)), sum)),
    v[c("within", "across", "vm")],
    a[c("ac", "cx", "anyal")],
    gender = d$gender[first],
    age = d$age[first],
    rc = .otrc(d$ra, d$rb, k),
    stringsAsFactors = FALSE
  )
  structure(
    out,
    class = c("otpyr", "data.frame"),
    years = sort(unique(d$year)),
    drop = sum(!ok),
    ord = ord
  )
}

#' @export
print.otpyr <- function(x, ...) {
  yrs <- attr(x, "years")
  cat(
    "otpyr: ",
    nrow(x),
    " person-years, fiscal years ",
    paste(range(yrs), collapse = "-"),
    " (",
    length(yrs),
    ")",
    if (attr(x, "drop") > 0) sprintf(", %d rows dropped", attr(x, "drop")),
    "\n",
    sep = ""
  )
  tab <- data.frame(
    year = yrs,
    py = as.vector(table(factor(x$year, yrs))),
    cx = as.vector(tapply(x$cx, factor(x$year, yrs), mean)),
    vm = as.vector(tapply(x$vm, factor(x$year, yrs), mean))
  )
  print(tab, row.names = FALSE, digits = 4)
  invisible(x)
}

# The column names, defaults overridden by cols.
.otcn <- function(cols, fn) {
  cn <- c(
    year = "EndFiscalYear",
    id = "UniqueIndividual_ID",
    ra = "Region_AtTimeOfPlacement",
    rb = "Region_MostRecentPlacement",
    gender = "Gender",
    age = "Age_Category",
    mh = "MentalHealth_Alert",
    sr = "SuicideRisk_Alert",
    sw = "SuicideWatch_Alert",
    np = "Number_Of_Placements"
  )
  if (is.null(cols)) {
    return(cn)
  }
  if (
    !is.character(cols) ||
      is.null(names(cols)) ||
      !all(names(cols) %in% names(cn))
  ) {
    .ckstop(fn, "`cols` must be named by %s", toString(names(cn)))
  }
  cn[names(cols)] <- cols
  cn
}

# The regional cluster of each person-year: its regions, sorted, joined.
.otrc <- function(ra, rb, k) {
  f <- factor(k, unique(k))
  r <- split(c(as.character(ra), as.character(rb)), c(f, f))
  vapply(
    r,
    function(z) paste(sort(unique(z), method = "radix"), collapse = " + "),
    "",
    USE.NAMES = FALSE
  )
}
