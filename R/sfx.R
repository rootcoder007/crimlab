#' Synthetic example data in the shape of published tables
#'
#' Draws a data frame with the columns of a published administrative table
#' and values of the kinds it holds, for examples, tests and teaching
#' where the real table cannot be shipped or fetched. Each column is
#' drawn on its own, uniformly from its values: nothing in the result is a
#' record of anyone, and no relation between columns is kept.
#'
#' A preset gives the columns of one table: the Ontario corrections
#' tables on segregation, restrictive confinement and deaths in custody
#' (`"otis-a01"`, `"otis-b01"` to `"otis-b09"`, `"otis-c01"` to
#' `"otis-c12"`, `"otis-d01"` to `"otis-d07"`) and on use of force in
#' correctional institutions (`"cuof-incidents"`, `"cuof-incident_type"`,
#' `"cuof-indigenous"`, `"cuof-inmate_incident"`,
#' `"cuof-inmate_participant"`, `"cuof-institution_summary"`,
#' `"cuof-location_summary"`, `"cuof-race"`, `"cuof-ethnic_origin"`,
#' `"cuof-religion"`, `"cuof-select_incident_summary"`,
#' `"cuof-staff_incident"`). A column of a preset takes one of a set of
#' labels, a whole number from a range, or a synthetic identifier
#' (`"syn-00001"` to `"syn-"` and `max(2, n %/% 4)` in five digits, so
#' that identifiers repeat as people's do). `"otis"` or `"cuof"` alone
#' gives the list of every table of the family, the `i`-th drawn with the
#' seed `(seed + i) %% 2^31`. `kind = NULL` lists the presets.
#'
#' Any other table is a named list, one element per column, each a vector
#' of the values the column takes (labels or numbers), drawn with
#' replacement; a value given once is a constant column.
#'
#' The draws are those of R's `sample.int(m, n, replace = TRUE)` after
#' `set.seed(seed)`, made by the package's own copy of R's generator, so
#' the result is the same on every system and in every language, and the
#' session's random numbers are not touched.
#'
#' @param kind The name of a preset or of a family of presets, a named
#'   list of the values of each column, or `NULL` for the list of presets.
#' @param n The number of rows: one whole number from 1 to 1e7.
#' @param seed The seed: one whole number from 0 to 2147483647.
#' @return A data frame of `n` rows; for a family, a named list of them
#'   (named `"a01"`, ..., or `"incidents"`, ...); for `kind = NULL`, a data
#'   frame with one row per column of every preset: `kind`, `column`,
#'   `type` (`"labels"`, `"range"` or `"id"`) and `values` (the labels
#'   separated by `"|"`, or the ends of the range).
#' @examples
#' head(synfix("otis-b02", n = 20, seed = 3))
#' x <- synfix("cuof", n = 5)
#' names(x)
#' synfix(list(year = 2019:2023, sex = c("F", "M")), n = 4, seed = 9)
#' head(synfix())
#' @export
synfix <- function(kind = NULL, n = 100, seed = 1) {
  f <- "synfix"
  n <- .cknum(n, "n", f, min = 1, max = 1e7, whole = TRUE)
  seed <- .cknum(seed, "seed", f, min = 0, max = 2147483647, whole = TRUE)
  tab <- as.data.frame(.Call(C_clab_r_sfxtab), stringsAsFactors = FALSE)
  if (is.null(kind)) {
    ty <- c("labels", "range", "id")[tab$type + 1L]
    va <- ifelse(tab$type == 1L, paste(tab$lo, tab$hi, sep = "|"), tab$lev)
    return(data.frame(kind = tab$kind, column = tab$column, type = ty,
                      values = va, stringsAsFactors = FALSE))
  }
  if (is.list(kind) && !is.data.frame(kind)) {
    return(.sfxdraw(.sfxlist(kind, f), n, seed))
  }
  if (!is.character(kind) || length(kind) != 1L || is.na(kind)) {
    .ckstop(f, paste("`kind` must be the name of a preset, a named list of",
                     "the values of each column, or NULL"))
  }
  fam <- sub("-.*", "", tab$kind)
  if (kind %in% fam) {
    ks <- unique(tab$kind[fam == kind])
    out <- lapply(seq_along(ks), function(i) {
      .sfxdraw(.sfxpre(tab, ks[i], n), n, (seed + i) %% 2^31)
    })
    return(stats::setNames(out, sub("^[^-]*-", "", ks)))
  }
  if (!kind %in% tab$kind) {
    .ckstop(f, "`kind` names no preset: \"%s\"; synfix() lists them", kind)
  }
  .sfxdraw(.sfxpre(tab, kind, n), n, seed)
}

# The columns of a preset as values to draw from.
.sfxpre <- function(tab, kind, n) {
  t <- tab[tab$kind == kind, ]
  v <- lapply(seq_len(nrow(t)), function(j) {
    switch(t$type[j] + 1L,
           strsplit(t$lev[j], "|", fixed = TRUE)[[1L]],
           seq.int(t$lo[j], t$hi[j]),
           sprintf("syn-%05d", seq_len(max(2, n %/% 4))))
  })
  stats::setNames(v, t$column)
}

# A named list of columns, each the values it takes.
.sfxlist <- function(x, f) {
  nm <- names(x)
  if (!length(x) || is.null(nm) || anyNA(nm) || !all(nzchar(nm)) ||
        anyDuplicated(nm)) {
    .ckstop(f, "a list `kind` must name each column once")
  }
  for (j in seq_along(x)) {
    v <- x[[j]]
    if (is.factor(v)) v <- as.character(v)
    if (!(is.character(v) || is.numeric(v) || is.logical(v)) ||
          !is.null(dim(v)) || !length(v) || anyNA(v)) {
      .ckstop(f, paste("column `%s` of `kind` must be a vector of the values",
                       "it takes, with no missing values"), nm[j])
    }
    x[[j]] <- v
  }
  x
}

# n rows drawn from the values of each column.
.sfxdraw <- function(cols, n, seed) {
  dn <- as.numeric(lengths(cols))
  k <- .Call(C_clab_r_synfix, dn, n, seed)
  out <- lapply(seq_along(cols), function(j) cols[[j]][k[, j] + 1L])
  as.data.frame(stats::setNames(out, names(cols)), stringsAsFactors = FALSE,
                check.names = FALSE)
}
