#' Presets: a published dataset's columns onto a generic function's arguments
#'
#' A preset says how the columns of one published dataset (an Ontario
#' segregation table, the Special Investigations Unit's cases, a use-of-force
#' file) feed the arguments of one of crimlab's generic functions, so that
#' `generic(clprep(data, preset))` gives the published analysis while the
#' same generic runs on anyone's data given as plain vectors. A preset
#' describes a layout and ships no data: the user supplies the table.
#'
#' The presets are rows of `presets.tsv` (in the installed package's `presets`
#' folder), one per preset and argument: `preset`, `generic`, `arg`, `column`
#' (the dataset's column), `recode`, `filter`, `units`, `note`, `source`. Python
#' reads the same file. A `recode` is one of: empty (the column as it is); `fac`
#' (a factor, one effect per value); `num` (numbers); `pct` (numbers over 100);
#' `one` (a 1 for every row); `flag` (1 for yes, y, true, t, 1 or a positive
#' number, 0 otherwise); `dum` (one 0 or 1 column per value but the first, named
#' by the column and the value); `sub:a,b` (characters a to b, as a year from an
#' ISO date); `const:v` (the argument is v, a number when it reads as one, no
#' column); or a mapping `x|y=1;has:abc=1;*=0` from labels, compared without
#' case or surrounding space (`has:` a label containing the text), to values,
#' `*` covering every other label and `*=NA` dropping those rows. A `filter`
#' `col=value` keeps the rows whose column has that value. An argument given by
#' several rows (or by `dum`) is a data frame of their columns, as `adj` always
#' is. An `arg` of `recipe` derives columns first, in the order of the rows:
#' `keep` drops the rows missing any of the columns named in `recode`; `long`
#' stacks the two count columns named in `recode` as `count` with `with` 1 for
#' the first; `acvm` and `pyr` reduce one row per placement to one per person
#' and year (see Details of the returned attributes), and `acnv` does so the
#' naive way (`acn`, the most alert flags on one row; `vmany`, 1 when any row
#' changes region; `hiacn`, `acn >= 2`); `pairb` adds `tb` (two or more of the
#' three alert flags) and `yb` (`Number_Of_Placements` of 2 or more); `pairc`
#' adds `tc` (the region at placement differs from the most recent) and `yc`
#' (`NumberConsecutiveDays_Segregation` capped at its 99th percentile); `gent`
#' takes the columns `inc,rent,colg,rentg` named in `recode` (and optionally the
#' two quantiles, `0.5,0.667` by default) and adds `gent`, a factor:
#' `ineligible` unless baseline income and rent are both at most their first
#' quantile, then `gentrified` when the growth in college graduates and in rent
#' both reach their second quantile, else `eligible`; `pop` takes
#' `count;rate;per` (separated by semicolons, as column names may hold commas)
#' and adds `pop`, the population `count / rate * per` that a published rate
#' implies.
#'
#' Rows with a missing value in any column the preset reads are dropped,
#' as the source analyses drop them; their number is the attribute
#' `drop`.
#'
#' @param data The dataset: a data frame with the columns the preset names.
#'   `NULL` (with `preset` `NULL`) lists the presets.
#' @param preset The preset's name, as `clprep()` lists them, or a preset
#'   of one's own: a data frame of rows laid out as those of `presets.tsv`
#'   (`preset`, `generic` and `arg` needed, the other columns optional);
#'   `NULL` for the list.
#' @return With a preset, an object of class `clpre`: a list of the
#'   generic's arguments, which the generic takes as its first argument,
#'   with attributes `preset`, `generic`, `source`, `units`, `note` and
#'   `drop`. Without one, a data frame of the presets: `preset`, `generic`,
#'   `args` and `source`.
#' @examples
#' head(clprep())
#' # a table laid out as the published OTIS c01 table (numbers made up)
#' d <- data.frame(EndFiscalYear = rep(2020:2022, each = 2),
#'                 Gender = rep(c("Female", "Male"), 3),
#'                 NumberIndividuals_RestrictiveConfinement =
#'                   c(40, 310, 52, 335, 47, 360))
#' p <- clprep(d, "otis-c01-cnt")
#' cntreg(p)
#' # the same model on plain vectors
#' cntreg(d$NumberIndividuals_RestrictiveConfinement, d$Gender == "Female",
#'        adj = factor(d$EndFiscalYear))
#' # a preset of one's own, for another agency's columns
#' mine <- data.frame(preset = "my-cnt", generic = "cntreg",
#'                    arg = c("y", "grp", "adj"),
#'                    column = c("count", "sex", "year"),
#'                    recode = c("", "f=1;*=0", "fac"))
#' e <- data.frame(count = c(4, 30, 6, 41), sex = c("f", "m", "f", "m"),
#'                 year = c(1, 1, 2, 2))
#' cntreg(clprep(e, mine))
#' @export
clprep <- function(data = NULL, preset = NULL) {
  tab <- .presets()
  if (is.data.frame(preset)) {
    tab <- .ptab(preset)
    preset <- tab$preset[1L]
  }
  if (is.null(preset)) {
    if (!is.null(data)) {
      stop("clprep: `preset` must name a preset.", call. = FALSE)
    }
    return(.plist(tab))
  }
  r <- .prows(tab, preset, data)
  rec <- r$arg == "recipe"
  for (i in which(rec)) {
    data <- .recipe(r$column[i], r$recode[i], data)
  }
  r <- r[!rec, , drop = FALSE]
  used <- unique(c(r$column[nzchar(r$column)], .fcols(r$filter)))
  miss <- setdiff(used, names(data))
  if (length(miss)) {
    stop("clprep: `data` has no column ", toString(miss), " (preset ",
         preset, ").", call. = FALSE)
  }
  ok <- rowSums(is.na(data[used])) == 0
  for (f in unlist(strsplit(r$filter[nzchar(r$filter)], ";", fixed = TRUE))) {
    kv <- strsplit(f, "=", fixed = TRUE)[[1L]]
    ok <- ok & trimws(as.character(data[[kv[1L]]])) == kv[2L]
  }
  # each argument's values; rows a mapping sends to NA are dropped after
  n <- nrow(data)
  vals <- lapply(seq_len(nrow(r)), function(i) {
    .recode(data[[r$column[i]]], r$recode[i], n)
  })
  for (v in vals) {
    if (length(v) == n) ok <- ok & !is.na(v)
  }
  structure(.pargs(r, vals, ok, n), class = "clpre", preset = preset,
            generic = r$generic[1L], source = r$source[1L],
            units = r$units[nzchar(r$units)][1L],
            note = r$note[nzchar(r$note)], drop = sum(!ok))
}

#' @export
print.clpre <- function(x, ...) {
  cat(
    "clpre: preset ",
    attr(x, "preset"),
    " for ",
    attr(x, "generic"),
    "()\n",
    "source: ",
    attr(x, "source"),
    "\n",
    sep = ""
  )
  cat(
    "arguments: ",
    toString(names(x)),
    "; rows dropped: ",
    attr(x, "drop"),
    "\n",
    sep = ""
  )
  invisible(x)
}

# The list of presets: one row each.
.plist <- function(tab) {
  p <- unique(tab$preset)
  first <- match(p, tab$preset)
  args <- vapply(p, function(k) {
    paste(unique(tab$arg[tab$preset == k]), collapse = " ")
  }, "")
  data.frame(preset = p, generic = tab$generic[first], args = unname(args),
             source = tab$source[first])
}

# The rows of one preset, with the checks of its name and of data.
.prows <- function(tab, preset, data) {
  if (!is.character(preset) || length(preset) != 1L || is.na(preset)) {
    stop("clprep: `preset` must be one name.", call. = FALSE)
  }
  r <- tab[tab$preset == preset, , drop = FALSE]
  if (!nrow(r)) {
    stop("clprep: no preset named ", preset, "; clprep() lists them.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("clprep: `data` must be a data frame.", call. = FALSE)
  }
  r
}

# The generic's arguments: each one's values on the rows kept (constants as
# they are), the adj columns as a data frame named by column.
.pargs <- function(r, vals, ok, n) {
  out <- list()
  for (a in unique(r$arg)) {
    j <- which(r$arg == a)
    v <- lapply(vals[j], function(x) {
      if (is.matrix(x)) {
        x[ok, , drop = FALSE]
      } else if (length(x) == n) {
        x[ok]
      } else {
        x
      }
    })
    if (a == "adj" || length(j) > 1L || is.matrix(v[[1L]])) {
      out[[a]] <- .pframe(v, r$column[j])
    } else {
      out[[a]] <- v[[1L]]
    }
  }
  out
}

# One data frame of the argument's columns v, named nm (a dum matrix gives
# one column per value, named column and value).
.pframe <- function(v, nm) {
  parts <- lapply(seq_along(v), function(i) {
    x <- v[[i]]
    d <- as.data.frame(x)
    names(d) <- if (is.matrix(x)) paste0(nm[i], colnames(x)) else nm[i]
    d
  })
  do.call(cbind, parts)
}

# A preset of one's own: the rows of one preset, laid out as presets.tsv.
.ptab <- function(p) {
  cols <- c("preset", "generic", "arg", "column", "recode", "filter",
            "units", "note", "source")
  if (!all(cols[1:3] %in% names(p)) || !nrow(p)) {
    stop("clprep: a preset table needs the columns preset, generic and ",
         "arg, and a row.", call. = FALSE)
  }
  out <- lapply(cols, function(k) {
    v <- if (k %in% names(p)) as.character(p[[k]]) else rep("", nrow(p))
    v[is.na(v)] <- ""
    v
  })
  out <- as.data.frame(out)
  names(out) <- cols
  if (length(unique(out$preset)) != 1L) {
    stop("clprep: a preset table must hold one preset.", call. = FALSE)
  }
  out
}

# The preset table, read with base R from the installed presets.tsv.
.presets <- function() {
  l <- readLines(
    system.file("presets", "presets.tsv", package = "crimlab"),
    encoding = "UTF-8"
  )
  m <- do.call(
    rbind,
    lapply(strsplit(l, "\t", fixed = TRUE), function(x) {
      length(x) <- 9L
      x[is.na(x)] <- ""
      x
    })
  )
  d <- as.data.frame(m[-1L, , drop = FALSE])
  names(d) <- m[1L, ]
  d
}

# The columns a filter string reads
.fcols <- function(f) {
  f <- as.character(unlist(strsplit(f[nzchar(f)], ";", fixed = TRUE)))
  vapply(strsplit(f, "=", fixed = TRUE), `[`, "", 1L)
}

# One argument's values from column x under the recode rule s (n rows).
.recode <- function(x, s, n) {
  if (startsWith(s, "const:")) {
    v <- substring(s, 7L)
    if (v %in% c("TRUE", "FALSE")) {
      return(as.logical(v))
    }
    num <- suppressWarnings(as.numeric(v))
    return(if (is.na(num)) v else num)
  }
  if (s == "one") {
    return(rep(1, n))
  }
  if (!nzchar(s)) {
    return(x)
  }
  if (s == "fac") {
    return(factor(x))
  }
  if (s == "num") {
    return(as.numeric(x))
  }
  if (s == "flag") {
    return(.flag(x))
  }
  if (s == "pct") {
    return(as.numeric(x) / 100)
  }
  if (s == "dum") {
    f <- factor(x)
    lv <- levels(f)[-1L]
    m <- matrix(0, length(f), length(lv), dimnames = list(NULL, lv))
    for (k in seq_along(lv)) m[, k] <- as.numeric(f == lv[k])
    return(m)
  }
  if (startsWith(s, "sub:")) {
    ab <- as.integer(strsplit(substring(s, 5L), ",", fixed = TRUE)[[1L]])
    return(substr(as.character(x), ab[1L], ab[2L]))
  }
  .remap(x, s)
}

# A mapping rule s ("a|b=1;has:c=2;*=0") applied to the labels x.
.remap <- function(x, s) {
  lab <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(lab))
  dflt <- NULL
  for (rule in strsplit(s, ";", fixed = TRUE)[[1L]]) {
    kv <- strsplit(rule, "=", fixed = TRUE)[[1L]]
    if (kv[1L] == "*") {
      dflt <- kv[2L]
      next
    }
    hit <- if (startsWith(kv[1L], "has:")) {
      grepl(tolower(substring(kv[1L], 5L)), lab, fixed = TRUE)
    } else {
      lab %in% strsplit(kv[1L], "|", fixed = TRUE)[[1L]]
    }
    out[is.na(out) & hit] <- kv[2L]
  }
  if (!is.null(dflt) && dflt != "NA") {
    out[is.na(out) & !is.na(lab)] <- dflt
  }
  if (is.null(dflt) && anyNA(out[!is.na(lab)])) {
    stop(
      "clprep: a value of the column has no recode: ",
      lab[!is.na(lab) & is.na(out)][1L],
      ".",
      call. = FALSE
    )
  }
  num <- suppressWarnings(as.numeric(out))
  if (anyNA(num[!is.na(out)])) out else num
}

# Recipes: derived columns a preset needs.
.recipe <- function(name, arg, d) {
  if (!name %in% c("long", "acvm", "pyr", "keep", "acnv", "pairb", "pairc",
                    "gent", "pop")) {
    stop("clprep: unknown recipe ", name, ".", call. = FALSE)
  }
  switch(
    name,
    long = {
      ab <- strsplit(arg, ",", fixed = TRUE)[[1L]]
      a <- d[setdiff(names(d), ab)]
      cbind(
        rbind(a, a),
        count = c(d[[ab[1L]]], d[[ab[2L]]]),
        with = rep(1:0, each = nrow(d))
      )
    },
    acvm = .acvm(d),
    pyr = .pyr(d),
    keep = .keep(d, strsplit(arg, ",", fixed = TRUE)[[1L]]),
    acnv = .acnv(d),
    pairb = .pairb(d),
    pairc = .pairc(d),
    gent = .gent(d, strsplit(arg, ",", fixed = TRUE)[[1L]]),
    pop = .pop(d, strsplit(arg, ";", fixed = TRUE)[[1L]])
  )
}

# pop: the population implied by a count and its published rate per `per`
.pop <- function(d, a) {
  .keep(d[0L, , drop = FALSE], a[1:2])
  d$pop <- as.numeric(d[[a[1L]]]) / as.numeric(d[[a[2L]]]) * as.numeric(a[3L])
  d
}

# the rows of d with none of the columns cols missing
.keep <- function(d, cols) {
  miss <- setdiff(cols, names(d))
  if (length(miss)) {
    stop("clprep: `data` has no column ", toString(miss), ".", call. = FALSE)
  }
  d[rowSums(is.na(d[cols])) == 0, , drop = FALSE]
}

# acnv: the naive reduction per person-year, the most alert flags on one
# row and whether any row changes region
.acnv <- function(d) {
  p <- .pyrows(d, c(.otcols, "Gender", "Age_Category"))
  x <- p$d
  nal <- .flag(x$MentalHealth_Alert) + .flag(x$SuicideRisk_Alert) +
    .flag(x$SuicideWatch_Alert)
  ch <- as.integer(as.character(x$Region_AtTimeOfPlacement) !=
                     as.character(x$Region_MostRecentPlacement))
  first <- vapply(p$g, `[`, 1L, 1L)
  out <- x[
    first,
    c("UniqueIndividual_ID", "EndFiscalYear", "Gender", "Age_Category")
  ]
  out$acn <- vapply(p$g, function(i) max(nal[i]), 1L)
  out$vmany <- vapply(p$g, function(i) max(ch[i]), 1L)
  out$hiacn <- as.integer(out$acn >= 2L)
  rownames(out) <- NULL
  out
}

# pairb: two or more alert flags, and two or more placements
.pairb <- function(d) {
  d <- .keep(d, c(.otcols[-(3:4)], "Number_Of_Placements"))
  d$tb <- as.integer(.flag(d$MentalHealth_Alert) +
                       .flag(d$SuicideRisk_Alert) +
                       .flag(d$SuicideWatch_Alert) >= 2L)
  np <- suppressWarnings(as.numeric(d$Number_Of_Placements))
  d$yb <- as.integer(!is.na(np) & np >= 2)
  d
}

# pairc: a change of region, and the segregation days capped at their 99th
# percentile
.pairc <- function(d) {
  d <- .keep(d, c("Region_AtTimeOfPlacement", "Region_MostRecentPlacement",
                  "NumberConsecutiveDays_Segregation"))
  d$tc <- as.integer(as.character(d$Region_AtTimeOfPlacement) !=
                       as.character(d$Region_MostRecentPlacement))
  y <- suppressWarnings(as.numeric(d$NumberConsecutiveDays_Segregation))
  y[is.na(y)] <- 0
  d$yc <- pmin(y, stats::quantile(y, 0.99, names = FALSE))
  d
}

# gent: the three-level gentrification factor of baseline-conditional
# designs (eligible areas poor at baseline; gentrified when both growth
# measures are high)
.gent <- function(d, a) {
  .keep(d[0L, , drop = FALSE], a[1:4])
  q <- if (length(a) > 4L) as.numeric(a[5:6]) else c(0.5, 0.667)
  v <- lapply(a[1:4], function(k) as.numeric(d[[k]]))
  cut <- vapply(1:4, function(k) {
    stats::quantile(v[[k]], q[(k + 1L) %/% 2L], na.rm = TRUE, names = FALSE)
  }, 0)
  el <- v[[1L]] <= cut[1L] & v[[2L]] <= cut[2L]
  gr <- v[[3L]] >= cut[3L] & v[[4L]] >= cut[4L]
  # 1 ineligible, 2 eligible, 3 gentrified; NA where ifelse() would give it
  lv <- c("ineligible", "eligible", "gentrified")
  d$gent <- factor(lv[1L + el + (el & gr)], lv)
  d
}

.otcols <- c(
  "UniqueIndividual_ID",
  "EndFiscalYear",
  "Region_AtTimeOfPlacement",
  "Region_MostRecentPlacement",
  "MentalHealth_Alert",
  "SuicideRisk_Alert",
  "SuicideWatch_Alert"
)

# alert flags as 0/1: yes, y, true, t or 1 (any case), a positive number, TRUE
.flag <- function(s) {
  if (is.logical(s)) {
    return(as.integer(s))
  }
  if (is.numeric(s)) {
    return(as.integer(!is.na(s) & s > 0))
  }
  as.integer(
    tolower(trimws(as.character(s))) %in% c("yes", "y", "true", "t", "1")
  )
}

# the rows of a placement table grouped by person and fiscal year, in the
# order of the first row of each, with each row's alert combination and
# its region changes (within the row, and from the previous row of the
# same person and year)
.pyrows <- function(d, need) {
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop(
      "clprep: `data` has no column ",
      toString(miss),
      ".",
      call. = FALSE
    )
  }
  d <- d[rowSums(is.na(d[need])) == 0, , drop = FALSE]
  d <- d[order(d$UniqueIndividual_ID, d$EndFiscalYear), , drop = FALSE]
  key <- paste(d$UniqueIndividual_ID, d$EndFiscalYear, sep = "\r")
  d$combo <- 4L *
    .flag(d$MentalHealth_Alert) +
    2L * .flag(d$SuicideRisk_Alert) +
    .flag(d$SuicideWatch_Alert)
  a <- as.character(d$Region_AtTimeOfPlacement)
  b <- as.character(d$Region_MostRecentPlacement)
  prev <- c(NA, a[-length(a)])
  prev[c(TRUE, key[-1L] != key[-length(key)])] <- NA
  d$vmrow <- as.integer(a != b) + as.integer(!is.na(prev) & a != prev)
  list(d = d, g = split(seq_len(nrow(d)), factor(key, unique(key))))
}

# acvm: alert complexity and regional volatility per person-year
.acvm <- function(d) {
  p <- .pyrows(d, c(.otcols, "Gender", "Age_Category"))
  first <- vapply(p$g, `[`, 1L, 1L)
  out <- p$d[
    first,
    c("UniqueIndividual_ID", "EndFiscalYear", "Gender", "Age_Category")
  ]
  out$ac <- vapply(p$g, function(i) length(unique(p$d$combo[i])), 1L)
  out$vm <- vapply(p$g, function(i) sum(p$d$vmrow[i]), 1L)
  out$hiac <- as.integer(out$ac >= 2L)
  rownames(out) <- NULL
  out
}

# pyr: placements, distinct regions and distinct alert combinations per
# person-year
.pyr <- function(d) {
  p <- .pyrows(d, .otcols)
  first <- vapply(p$g, `[`, 1L, 1L)
  out <- p$d[first, c("UniqueIndividual_ID", "EndFiscalYear")]
  out$np <- lengths(p$g)
  out$nreg <- vapply(
    p$g,
    function(i) {
      length(unique(c(
        as.character(p$d$Region_AtTimeOfPlacement[i]),
        as.character(p$d$Region_MostRecentPlacement[i])
      )))
    },
    1L
  )
  out$ac <- vapply(p$g, function(i) length(unique(p$d$combo[i])), 1L)
  rownames(out) <- NULL
  out
}

# The generic f called on a clpre list p from clprep(): p's elements fill
# f's arguments by name, and the arguments of the original call cl (but the
# first, p itself), evaluated in env, override them.
.prep <- function(f, p, cl, env) {
  a <- lapply(as.list(cl)[-(1:2)], eval, envir = env)
  b <- unclass(p)
  attributes(b) <- list(names = names(b))
  b[names(a)] <- a
  do.call(f, b)
}
