#' Check a table against a preset before use
#'
#' `prechk()` reads a table the way [clprep()] would under a preset and
#' reports, in one error, every problem that would stop it or quietly change
#' its result: the columns the preset reads (directly, in a filter or in a
#' recipe) that the table lacks, the labels of a column outside the codes
#' the publisher lists for it, and the labels a mapping recode with no `*`
#' rule does not cover. A table that passes is returned unchanged, so the
#' check can sit in front of the call: `clprep(prechk(d, p), p)`.
#'
#' The codes of the shipped presets are rows of `codes.tsv` (in the
#' installed package's `presets` folder): `agency`, `column`, `value` and
#' `source`, one row per label. A preset named `otis-...` is checked against
#' the rows of agency `otis`: the labels of the data dictionary of the
#' Ontario Ministry of the Solicitor General for its segregation,
#' restrictive confinement and deaths in custody tables, one list per column
#' over all the tables that carry it (so `Gender` allows `All`, which only
#' the tables of medians and modes use). Labels are compared without case or
#' surrounding space, as `clprep()` compares them, and missing values are
#' not codes (`clprep()` drops their rows). Columns with no codes listed,
#' counts and days among them, are not checked for labels.
#'
#' The presets describe the layouts of published tables, aggregate counts
#' and de-identified records; crimlab ships none of the tables.
#'
#' @param data The table: a data frame.
#' @param preset The preset's name, as `clprep()` lists them, or a preset
#'   of one's own laid out as the rows of `presets.tsv` (see [clprep()]).
#' @param codes The labels each column may hold: a data frame with columns
#'   `column` and `value`, one row per label. `NULL` (the default) takes
#'   the rows of `codes.tsv` whose `agency` is the preset's name up to its
#'   first `-`; none for a preset of one's own unless given here.
#' @return `data`, invisibly, when it passes; otherwise an error naming each
#'   missing column and each unknown label with its column.
#' @examples
#' # a table laid out as the published OTIS b03 table, simulated
#' d <- data.frame(
#'   EndFiscalYear = rep(2021:2022, each = 4),
#'   Region_AtTimeOfPlacement = rep(c("Central", "Toronto"), 4),
#'   Institution_AtTimeOfPlacement = rep(c("Quinte Detention Centre",
#'                                         "Toronto South Detention Centre"),
#'                                       4),
#'   Alert_Type = "Mental Health Alert",
#'   Alert_Presence = rep(c("Yes", "Yes", "No", "No"), 2),
#'   Number_SegregationPlacements = c(31, 52, 140, 260, 28, 61, 133, 248)
#' )
#' prechk(d, "otis-b03-cnt")
#' # a region the publisher does not list, and a column missing
#' bad <- d[-6]
#' bad$Region_AtTimeOfPlacement[1] <- "Southern"
#' try(prechk(bad, "otis-b03-cnt"))
#' # a preset of one's own, and a label its recode does not cover
#' mine <- data.frame(preset = "my-cnt", generic = "cntreg",
#'                    arg = c("y", "grp"), column = c("count", "sex"),
#'                    recode = c("", "f=1;m=0"))
#' e <- data.frame(count = c(4, 30, 6, 41), sex = c("f", "m", "f", "x"))
#' try(prechk(e, mine))
#' @export
prechk <- function(data, preset, codes = NULL) {
  r <- .clas("prechk", .pkrow(preset))
  if (!is.data.frame(data)) {
    stop("prechk: `data` must be a data frame.", call. = FALSE)
  }
  .cktxtdf(data, "data", "prechk")
  name <- r$preset[1L]
  codes <- .pkcodes(codes, name)
  rec <- r$arg == "recipe"
  bad <- setdiff(r$column[rec], c(names(.rout), "keep"))
  if (length(bad)) {
    stop("prechk: unknown recipe ", toString(bad), " (preset ", name, ").",
         call. = FALSE)
  }
  need <- .pkneed(r, rec)
  msg <- .pkmsg(data, r, rec, need, codes)
  if (length(msg)) {
    stop("prechk: the table does not fit preset ", name, ":\n",
         paste0("* ", msg, collapse = "\n"), call. = FALSE)
  }
  invisible(data)
}

#' prechk(): the rows of the preset, named or given as a table
#' @noRd
.pkrow <- function(preset) {
  if (is.data.frame(preset)) {
    return(.ptab(preset))
  }
  if (!is.character(preset) || length(preset) != 1L || is.na(preset)) {
    stop("prechk: `preset` must be one name or a preset table.",
         call. = FALSE)
  }
  tab <- .presets()
  r <- tab[tab$preset == preset, , drop = FALSE]
  if (!nrow(r)) {
    stop("prechk: no preset named ", preset, "; clprep() lists them.",
         call. = FALSE)
  }
  r
}

#' prechk(): the publisher's codes, by default those of the preset's agency
#' @noRd
.pkcodes <- function(codes, name) {
  if (is.null(codes)) {
    k <- .presets("codes.tsv")
    codes <- k[k$agency == sub("-.*$", "", name), , drop = FALSE]
  } else if (!is.data.frame(codes) ||
               !all(c("column", "value") %in% names(codes))) {
    stop("prechk: `codes` must be a data frame with columns column and ",
         "value.", call. = FALSE)
  }
  codes
}

#' prechk(): the columns the preset reads from the table
#' @noRd
.pkneed <- function(r, rec) {
  ins <- unlist(Map(.rin, r$column[rec], r$recode[rec]), use.names = FALSE)
  ins <- ins[!is.na(ins) & nzchar(ins)]
  outs <- unlist(.rout[r$column[rec]])
  read <- setdiff(c(r$column[!rec & nzchar(r$column)], .fcols(r$filter)),
                  outs)
  unique(c(ins, read))
}

#' prechk(): what does not fit: missing columns, labels outside the codes,
#' labels a mapping leaves without a value
#' @noRd
.pkmsg <- function(data, r, rec, need, codes) {
  msg <- character()
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    msg <- paste0("`data` has no column ", toString(miss), ".")
  }
  # labels outside the publisher's codes
  for (col in intersect(intersect(need, names(data)), codes$column)) {
    ok <- tolower(trimws(as.character(codes$value[codes$column == col])))
    msg <- c(msg, .unk(data[[col]], col, "not among its codes",
                       function(l) l %in% ok))
  }
  # labels a mapping without a `*` rule leaves without a value
  map <- !rec & grepl("=", r$recode, fixed = TRUE) &
    !startsWith(r$recode, "const:") &
    !grepl("(^|;)\\*=", r$recode) & r$column %in% names(data)
  for (i in which(map)) {
    msg <- c(msg, .unk(data[[r$column[i]]], r$column[i], "with no recode",
                       .known(r$recode[i])))
  }
  msg
}

#' The labels of x (column col) that known(), given labels in lower case
#' without surrounding space, rejects: one line of prechk()'s message, or
#' none when all are known.
#' @noRd
.unk <- function(x, col, why, known) {
  u <- unique(trimws(as.character(x[!is.na(x)])))
  u <- u[!known(tolower(u))]
  if (!length(u)) {
    return(character())
  }
  sprintf("column %s has labels %s: %s.", col, why, toString(u))
}

#' The labels a mapping rule s covers, as clprep() matches them: a function
#' of labels in lower case.
#' @noRd
.known <- function(s) {
  keys <- vapply(strsplit(strsplit(s, ";", fixed = TRUE)[[1L]], "=",
                          fixed = TRUE), `[`, "", 1L)
  has <- startsWith(keys, "has:")
  eq <- unlist(strsplit(keys[!has], "|", fixed = TRUE))
  part <- tolower(substring(keys[has], 5L))
  function(l) {
    vapply(l, function(v) {
      v %in% eq || any(vapply(part, grepl, NA, x = v, fixed = TRUE))
    }, NA, USE.NAMES = FALSE)
  }
}

#' The columns of the table a recipe reads, from its name and its recode.
#' @noRd
.rin <- function(name, arg) {
  sp <- function(s) strsplit(arg, s, fixed = TRUE)[[1L]]
  switch(
    name,
    long = ,
    keep = sp(","),
    gent = sp(",")[1:4],
    pop = sp(";")[1:2],
    pyr = .otcols,
    acvm = ,
    acnv = c(.otcols, "Gender", "Age_Category"),
    pairb = c(.otcols[-(3:4)], "Number_Of_Placements"),
    pairc = c("Region_AtTimeOfPlacement", "Region_MostRecentPlacement",
              "NumberConsecutiveDays_Segregation")
  )
}

# The columns each recipe adds (see clprep()).
.rout <- list(
  long = c("count", "with"),
  acvm = c("ac", "vm", "hiac"),
  acnv = c("acn", "vmany", "hiacn"),
  pyr = c("np", "nreg", "ac"),
  pairb = c("tb", "yb"),
  pairc = c("tc", "yc"),
  gent = "gent",
  pop = "pop"
)
