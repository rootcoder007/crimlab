# Categorical integrity: category mappings made explicit and checked.
# Numeric codes imported without their value labels, positional recodes
# that reassign groups wholesale and alphabetical releveling that changes
# the reference category can each relabel whole demographic groups
# without a warning. These functions make every recode name-based, refuse
# anything unmapped, prove a recode with a before/after cross-tabulation
# and check recoded counts against the counts a release published.

#' Recode categories by name
#'
#' `grdrcd()` recodes a categorical column by an explicit mapping from old
#' to new labels and refuses any observed value the mapping does not name:
#' silent pass-through is how group labels get corrupted. `deccodes()`
#' decodes numeric category codes with a code-to-label dictionary (codes
#' as names), the step that turns the codes of an imported file into its
#' categories. `relbl()` relabels a variable by name, never by position,
#' and returns a factor whose levels follow the mapping.
#'
#' Assigning labels by position is how category codes get rotated: when
#' the labels are listed in alphabetical order and the codes are not,
#' every group receives its neighbour's label. Mapping by name makes the
#' assignment independent of any order.
#'
#' In `grdrcd()` a value listed in `keep` passes through unchanged and
#' missing values stay missing. In `relbl()` every level (the levels of a
#' factor, or the sorted distinct values) must be mapped or kept; the
#' levels of the result are the new labels in the order of the mapping,
#' then the kept values that occur. The result carries an attribute
#' `recode`: the mapping and the kept values.
#'
#' @param x The labels or codes: a vector or a factor.
#' @param map,vlabs A named character vector: old label (or code) =
#'   new label.
#' @param keep Values that pass through unchanged.
#' @param keepna `FALSE` refuses missing codes.
#' @return `grdrcd()` and `deccodes()`: a character vector the length of
#'   `x`. `relbl()`: a factor.
#' @seealso [grdlvls()], [declbld()], [vrfrcd()], [vrfmrg()].
#' @examples
#' # simulated race codes of an imported file
#' codes <- c(1, 2, 2, 3, 1, NA)
#' deccodes(codes, c("1" = "White", "2" = "Black", "3" = "Indigenous"))
#' grdrcd(c("M", "F", "F", "X"), c(M = "Male", F = "Female"), keep = "X")
#' relbl(factor(c("a", "b", "a")), c(a = "Adult", b = "Youth"))
#' try(grdrcd(c("M", "F", "U"), c(M = "Male", F = "Female")))
#' @name grdrcd
#' @export
grdrcd <- function(x, map, keep = character()) {
  x <- .catx(x, "grdrcd")
  .catmap(map, "grdrcd")
  .crec(as.character(x), map, as.character(keep), "grdrcd")
}

#' @rdname grdrcd
#' @export
deccodes <- function(x, vlabs, keepna = TRUE) {
  x <- .catx(x, "deccodes")
  .catmap(vlabs, "deccodes", "vlabs")
  if (!isTRUE(keepna) && anyNA(x)) {
    stop("deccodes: `x` has missing codes and `keepna` is FALSE.",
         call. = FALSE)
  }
  .crec(as.character(x), vlabs, character(), "deccodes")
}

#' @rdname grdrcd
#' @export
relbl <- function(x, map, keep = character()) {
  .rlargs(x, map)
  old <- if (is.factor(x)) {
    levels(x)
  } else {
    as.character(sort(unique(x[!is.na(x)]), method = "radix"))
  }
  keep <- as.character(keep)
  un <- setdiff(old, c(names(map), keep))
  if (length(un)) {
    stop("relbl: levels with no mapping: ", toString(.cq(un)),
         ". Name every level or list it in `keep`.", call. = FALSE)
  }
  map <- stats::setNames(as.character(map), names(map))
  out <- .crec(as.character(x), map, keep, "relbl")
  lev <- unique(c(unname(map), keep))
  f <- factor(as.character(out), levels = lev[lev %in% out | lev %in% map])
  attr(f, "recode") <- attr(out, "recode")
  f
}

#' A factor with declared, verified levels
#'
#' Builds a factor whose levels are exactly the ones declared, in the
#' declared order, and refuses any value outside them. A declared
#' reference must be the first level, so the group every model contrast is
#' measured against is stated rather than left to alphabetical order.
#'
#' @param x The labels: a vector or a factor.
#' @param lvls The levels, a non-empty character vector without `NA`, the
#'   reference first.
#' @param ref `NULL`, or the reference level, which must be `lvls[1]`.
#' @return A factor with levels `lvls`.
#' @seealso [grdrcd()].
#' @examples
#' grdlvls(c("Male", "Female", "Male"), c("Male", "Female"), ref = "Male")
#' try(grdlvls(c("Male", "Female", "X"), c("Male", "Female")))
#' try(grdlvls("Male", c("Male", "Female"), ref = "Female"))
#' @export
grdlvls <- function(x, lvls, ref = NULL) {
  x <- .catx(x, "grdlvls")
  if (!is.character(lvls) || !length(lvls) || anyNA(lvls)) {
    stop("grdlvls: `lvls` must be a non-empty character vector.",
         call. = FALSE)
  }
  x <- as.character(x)
  stray <- setdiff(unique(x[!is.na(x)]), lvls)
  if (length(stray)) {
    stop("grdlvls: values outside the declared levels: ",
         toString(.cq(stray)), ".", call. = FALSE)
  }
  if (!is.null(ref) && !identical(as.character(ref), lvls[1L])) {
    stop("grdlvls: the declared reference ", .cq(ref), " is not the first ",
         "level (", .cq(lvls[1L]), "); put the reference first in `lvls`.",
         call. = FALSE)
  }
  factor(x, levels = lvls)
}

#' Decode a labelled import by its value labels
#'
#' Decodes a column of codes by the value labels it carries (the `labels`
#' attribute that readers of SPSS, Stata and SAS files attach, a named
#' vector label = code), or by `vlabs` (code = label) when it carries
#' none. The codes alone are not the categories: decoding by code keeps
#' every observation with its own label whatever order the labels are
#' listed in.
#'
#' The levels of the result are the labels in the order of their codes
#' (numerically when the codes are numbers). A column with neither a
#' `labels` attribute nor `vlabs` is refused.
#'
#' @param x The codes, with or without a `labels` attribute.
#' @param vlabs `NULL`, or a named character vector code = label, used when
#'   `x` has no `labels` attribute.
#' @return A factor.
#' @seealso [deccodes()], [xfrvrf()].
#' @examples
#' x <- structure(c(2, 1, 3, 2), labels = c(Female = 1, Male = 2, Other = 3))
#' declbld(x)
#' declbld(c(10, 20, 10), c("20" = "Remand", "10" = "Sentenced"))
#' try(declbld(c(1, 2)))
#' @export
declbld <- function(x, vlabs = NULL) {
  if (!is.atomic(x) || is.null(x)) {
    stop("declbld: `x` must be a vector of codes.", call. = FALSE)
  }
  lab <- attr(x, "labels", exact = TRUE)
  if (!is.null(lab)) {
    if (is.null(names(lab))) {
      stop("declbld: the `labels` attribute has no names.", call. = FALSE)
    }
    vlabs <- stats::setNames(names(lab), as.character(unname(lab)))
  }
  if (is.null(vlabs)) {
    stop("declbld: `x` carries no `labels` attribute and `vlabs` was not ",
         "given; the codes alone are not the categories.", call. = FALSE)
  }
  .catmap(vlabs, "declbld", "vlabs")
  out <- .crec(as.character(unclass(x)), vlabs, character(), "declbld")
  o <- order(suppressWarnings(as.numeric(names(vlabs))), names(vlabs),
             method = "radix")
  f <- factor(as.character(out), levels = unique(unname(vlabs[o])))
  attr(f, "recode") <- attr(out, "recode")
  f
}

#' Refuse a categorical column where a 0/1 indicator is required
#'
#' Checks that a column meant as a numeric 0/1 indicator (a treatment, an
#' outcome) is one. A factor or text column is refused rather than
#' converted: `as.numeric()` on a factor returns the level indices (1, 2,
#' ...), not the data, and a mis-ordered level relabels every
#' observation. Encode it explicitly first, for example
#' `as.integer(x == "treated")`.
#'
#' @param x The column.
#' @param col Its name, for the message.
#' @return `TRUE`, invisibly; an error otherwise.
#' @seealso [grdrcd()].
#' @examples
#' grdbin(c(0, 1, 1, NA), "custody")
#' try(grdbin(factor(c("no", "yes")), "custody"))
#' try(grdbin(c(0, 1, 2), "custody"))
#' @export
grdbin <- function(x, col) {
  if (is.null(x) || !is.atomic(x)) {
    stop("grdbin: column `", col, "` must be an atomic vector.", call. = FALSE)
  }
  if (is.factor(x) || is.character(x)) {
    lv <- if (is.factor(x)) levels(x) else unique(as.character(x))
    stop("grdbin: column ", .cq(col), " is categorical (",
         toString(.cq(utils::head(lv, 4))), "...). It is not coerced: ",
         "as.numeric() on a factor returns level indices, not the data. ",
         "Encode it first, for example as.integer(x == \"treated\").",
         call. = FALSE)
  }
  ux <- unique(x[!is.na(x)])
  if (!all(ux %in% c(0, 1))) {
    stop("grdbin: column ", .cq(col), " must be 0 or 1 (saw: ",
         toString(utils::head(ux, 5)), ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Prove a recode with a before/after cross-tabulation
#'
#' Cross-tabulates the original and the recoded column and checks that the
#' recode is a function of the original label (no category was split
#' between two new ones), that it follows the declared mapping (a value
#' the mapping does not name must be unchanged), that no row was lost or
#' added and that missingness did not change. This is how swapped groups
#' are caught before any model runs.
#'
#' @param orig,recd The column before and after the recode.
#' @param decl The declared mapping, a named vector old = new.
#' @return The cross-tabulation as a data frame (`original`, `recoded`,
#'   `Freq`), invisibly; an error at the first failed check.
#' @seealso [grdrcd()], [vrfmrg()].
#' @examples
#' o <- c("M", "F", "F", NA)
#' m <- c(M = "Male", F = "Female")
#' vrfrcd(o, grdrcd(o, m), m)
#' try(vrfrcd(o, c("Female", "Male", "Male", NA), m))
#' @export
vrfrcd <- function(orig, recd, decl) {
  orig <- .catx(orig, "vrfrcd", "orig")
  recd <- .catx(recd, "vrfrcd", "recd")
  if (length(orig) != length(recd)) {
    stop("vrfrcd: length mismatch (", length(orig), " vs ", length(recd),
         "): rows were lost or duplicated during the recode.", call. = FALSE)
  }
  if (!identical(is.na(orig), is.na(recd))) {
    stop("vrfrcd: missingness changed during the recode (values silently ",
         "became NA, or NAs were filled).", call. = FALSE)
  }
  ok <- !is.na(orig)
  tab <- table(original = as.character(orig[ok]),
               recoded = as.character(recd[ok]))
  fan <- rowSums(tab > 0)
  if (any(fan > 1L)) {
    stop("vrfrcd: original category mapped to MULTIPLE new categories: ",
         toString(.cq(names(fan)[fan > 1L])), ". The recode is not a ",
         "function of the category label.", call. = FALSE)
  }
  for (old in rownames(tab)) {
    got <- colnames(tab)[tab[old, ] > 0]
    want <- if (old %in% names(decl)) unname(decl[old]) else old
    if (!identical(got, as.character(want))) {
      stop("vrfrcd: ", .cq(old), " was mapped to ", .cq(got), " but the ",
           "declared mapping says ", .cq(want), ". THIS is how groups get ",
           "swapped; fix the recode before any model runs.", call. = FALSE)
    }
  }
  invisible(as.data.frame(tab))
}

#' Verify recoded counts against the counts a release published
#'
#' Counts each label of `x` and compares the counts with those the source
#' published (a table's margins, an annual report). When they differ and
#' some permutation of the labels would make them agree (up to seven
#' labels are searched), the labels are attached to the wrong groups, and
#' the message says which.
#'
#' A label of `x` that the published counts do not have is a failure.
#' With `strict` (the default) a failure is an error; otherwise the check
#' is returned with `ok` `FALSE` and its `message`.
#'
#' @param x The recoded labels.
#' @param pub The published counts, a named numeric vector (label = count).
#' @param tol The largest difference counted as agreement, one number from
#'   0.
#' @param strict `TRUE` stops at a failure.
#' @return A list of class `vrfmrg` (invisibly when the counts agree):
#'   `counts`, `published`, `ok`, `permutation` (the labels that would fit,
#'   or `NULL`) and `message`.
#' @seealso [vrfrcd()], [xfrvrf()].
#' @examples
#' x <- rep(c("White", "Black", "Indigenous"), c(50, 30, 20))
#' vrfmrg(x, c(White = 50, Black = 30, Indigenous = 20))
#' vrfmrg(x, c(White = 30, Black = 50, Indigenous = 20), strict = FALSE)
#' @export
vrfmrg <- function(x, pub, tol = 0, strict = TRUE) {
  x <- .catx(x, "vrfmrg")
  if (inherits(x, c("Date", "POSIXt"))) {
    x <- as.character(x)
  }
  strict <- isTRUE(strict)
  .mrgargs(pub, tol)
  labs <- names(pub)
  obs <- vapply(labs, function(l) sum(!is.na(x) & x == l), 0)
  pubv <- pub[labs]
  out <- structure(list(counts = obs, published = pubv, ok = FALSE,
                        permutation = NULL, message = NULL),
                   class = "vrfmrg")
  extra <- setdiff(unique(x[!is.na(x)]), labs)
  if (length(extra)) {
    out$message <- paste0("vrfmrg: labels present in the data but not in the ",
                          "published counts: ", toString(.cq(extra)), ".")
    if (strict) {
      stop(out$message, call. = FALSE)
    }
    return(out)
  }
  out$ok <- all(abs(obs - pubv) <= tol)
  if (out$ok) {
    return(invisible(out))
  }
  out <- .mrgfail(out, labs, tol)
  if (strict) {
    stop(out$message, call. = FALSE)
  }
  out
}

#' @export
print.vrfmrg <- function(x, ...) {
  cat("vrfmrg:", if (isTRUE(x$ok)) "the counts match" else "the counts differ",
      "\n")
  print(data.frame(label = names(x$counts), observed = unname(x$counts),
                   published = unname(x$published),
                   difference = unname(x$counts - x$published)),
        row.names = FALSE)
  if (!is.null(x$message)) cat(strwrap(x$message), sep = "\n")
  invisible(x)
}

#' Verify a categorical variable that crossed from one program to another
#'
#' Checks a categorical column after it was exported by one program and
#' imported by another: the value labels that arrived must agree with the
#' source's code book, and the decoded values must reproduce the counts
#' the source reported. When the column carries a `labels` attribute it
#' is decoded by it ([declbld()]); otherwise by `codebook`; otherwise the
#' values are taken as labels.
#'
#' @param imported The column as imported.
#' @param soucnts The counts the source program reported, a named numeric
#'   vector (label = count).
#' @param codebook `NULL`, or the source's code book, a named vector code =
#'   label.
#' @inheritParams vrfmrg
#' @return A list of class `xfrvrf`: `ok`, `decoded` (a factor), `marginals`
#'   (the [vrfmrg()] check), `codebook` (`TRUE`, `FALSE`, or `NA` when there
#'   was nothing to compare) and `reasons`.
#' @seealso [declbld()], [vrfmrg()].
#' @examples
#' imp <- structure(c(1, 2, 2, 1, 1), labels = c(No = 1, Yes = 2))
#' xfrvrf(imp, c(No = 3, Yes = 2), codebook = c("1" = "No", "2" = "Yes"))
#' xfrvrf(c("No", "Yes", "Yes"), c(No = 1, Yes = 2))
#' @export
xfrvrf <- function(imported, soucnts, codebook = NULL, tol = 0, strict = TRUE) {
  strict <- isTRUE(strict)
  if (!is.null(codebook) && (!is.atomic(codebook) || is.null(names(codebook)) ||
                               inherits(codebook, c("Date", "POSIXt")))) {
    stop("xfrvrf: `codebook` must be a named vector (code = label), or NULL.",
         call. = FALSE)
  }
  lab <- attr(imported, "labels", exact = TRUE)
  cbok <- if (!is.null(lab) && !is.null(codebook)) {
    .xfcb(lab, codebook, strict)
  } else {
    NA
  }
  dec <- .clas("xfrvrf", .xfdec(imported, lab, codebook, soucnts))
  m <- .clas("xfrvrf", vrfmrg(as.character(dec), soucnts, tol = tol,
                              strict = strict))
  ok <- isTRUE(m$ok) && !isFALSE(cbok)
  reasons <- c(if (isFALSE(cbok)) {
    "the value labels disagree with the code book"
  },
               if (!isTRUE(m$ok)) m$message)
  structure(list(ok = ok, decoded = dec, marginals = m, codebook = cbok,
                 reasons = reasons), class = "xfrvrf")
}

#' @export
print.xfrvrf <- function(x, ...) {
  cat("xfrvrf:", if (isTRUE(x$ok)) "as the source reported" else
    "does not match the source", "\n")
  cat("code book:", if (is.na(x$codebook)) "nothing to compare" else
    if (isTRUE(x$codebook)) "agrees" else "disagrees", "\n")
  print(x$marginals)
  if (length(x$reasons)) cat(paste("failed:", x$reasons), sep = "\n")
  invisible(x)
}

.cq <- function(v) paste0("'", v, "'")

.catx <- function(x, fn, nm = "x") {
  if (is.null(x) || !is.atomic(x)) {
    stop(fn, ": `", nm, "` must be a vector of labels.", call. = FALSE)
  }
  if (is.factor(x)) as.character(x) else x
}

.catmap <- function(map, fn, nm = "map") {
  if (!is.character(map) || is.null(names(map)) || anyNA(names(map)) ||
        !all(nzchar(names(map)))) {
    stop(fn, ": `", nm, "` must be a named character vector (old = new).",
         call. = FALSE)
  }
}

.crec <- function(x, map, keep, fn) {
  seen <- unique(x[!is.na(x)])
  un <- setdiff(seen, c(names(map), keep))
  if (length(un)) {
    stop(fn, ": values with no mapping: ", toString(.cq(un)), ". Every ",
         "observed category must be mapped by name (or listed in `keep`).",
         call. = FALSE)
  }
  out <- ifelse(x %in% names(map), unname(map[x]), x)
  out[is.na(x)] <- NA_character_
  attr(out, "recode") <- list(mapping = map, kept = keep)
  out
}

.mrgargs <- function(pub, tol) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("vrfmrg: `tol` must be one number from 0.", call. = FALSE)
  }
  if (!is.numeric(pub) || !length(pub) || anyNA(pub) || is.null(names(pub)) ||
        anyNA(names(pub)) || !all(nzchar(names(pub)))) {
    stop("vrfmrg: `pub` must be a named numeric vector.", call. = FALSE)
  }
}

# vrfmrg(): counts that differ, with the permutation of the labels that
# would make them agree (searched for up to seven labels) and the message
.mrgfail <- function(out, labs, tol) {
  obs <- out$counts
  pubv <- out$published
  if (length(labs) <= 7L) {
    out["permutation"] <- list(.mrgperm(obs, pubv, labs, tol))
  }
  hint <- ""
  p <- out$permutation
  if (!is.null(p)) {
    mv <- names(p)[p != names(p)]
    hint <- paste0(" The observed counts match the published ones if the ",
                   "labels are permuted (",
                   toString(paste(mv, "->", p[mv])),
                   "): the labels are attached to the wrong groups.")
  }
  out$message <- paste0(
    "vrfmrg: recoded counts do not match the published counts. ",
    paste(sprintf("%s: observed %s, published %s", labs, obs, pubv),
          collapse = "; "), ".", hint)
  out
}

# the first permutation p of labs (in the order of the recursion) under
# which the observed counts, relabelled by p, match the published ones
.mrgperm <- function(obs, pub, labs, tol) {
  perms <- function(v) {
    if (length(v) <= 1L) {
      return(list(v))
    }
    unlist(lapply(seq_along(v), function(i) {
      lapply(perms(v[-i]), function(r) c(v[i], r))
    }), recursive = FALSE)
  }
  for (p in perms(labs)) {
    r <- stats::setNames(obs, p)[labs]
    if (!identical(p, labs) && all(abs(r - pub) <= tol)) {
      return(stats::setNames(p, labs))
    }
  }
  NULL
}

# relbl(): `x` a vector of labels, `map` named
.rlargs <- function(x, map) {
  if (is.null(x) || !is.atomic(x) || is.complex(x) || is.raw(x)) {
    stop("relbl: `x` must be a vector of labels (text, a factor or ",
         "numbers).", call. = FALSE)
  }
  .cktxt(if (is.factor(x)) levels(x) else x, "x", "relbl")
  if (is.null(names(map)) || anyNA(names(map)) || !all(nzchar(names(map)))) {
    stop("relbl: `map` must be named (old label = new label); labels are ",
         "never assigned by position.", call. = FALSE)
  }
}

# xfrvrf(): do the value labels that arrived agree with the code book?
.xfcb <- function(lab, codebook, strict) {
  got <- stats::setNames(names(lab), as.character(unname(lab)))
  want <- stats::setNames(as.character(codebook), names(codebook))
  bad <- setdiff(union(names(got), names(want)),
                 names(got)[names(got) %in% names(want) &
                              got == want[names(got)]])
  cbok <- !length(bad)
  if (!cbok && strict) {
    stop("xfrvrf: the value labels that arrived disagree with the source ",
         "code book at codes ", toString(bad), ": arrived ",
         toString(paste0(bad, "=", got[bad])), "; source ",
         toString(paste0(bad, "=", want[bad])), ".", call. = FALSE)
  }
  cbok
}

# xfrvrf(): the column decoded by its labels, the code book, or as labels
.xfdec <- function(imported, lab, codebook, soucnts) {
  if (!is.null(lab)) {
    declbld(imported)
  } else if (!is.null(codebook)) {
    declbld(imported,
            stats::setNames(as.character(codebook), names(codebook)))
  } else {
    v <- as.character(imported)
    factor(v, levels = unique(c(names(soucnts), v[!is.na(v)])))
  }
}
