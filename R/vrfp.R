#' Run the pollution-to-health pipeline and report it
#'
#' Resolves the exposure (simulated demo data, a CSV with an \code{exposure}
#' column and an optional \code{income} column, a NAPS-style pull with a
#' \code{value} column and a \code{unit}, or scalar arguments), logs the
#' assumptions, and when they all hold runs the concentration-response,
#' attributable-fraction, displaced-mortality, burden and equity stages.
#' \code{polrpttxt()} renders the result as plain text. It is the one-call
#' version of the air pollution functions, made for scripts: for bad data or
#' a value it cannot use it does not stop with an R error, but returns a
#' report whose status says what went wrong. Only \code{poll}, which names
#' the report, is checked first and refused with an R error (as the command
#' line refuses it) when it is not one pollutant name: missing, not text,
#' \code{NULL}, several names, or text that is not valid UTF-8.
#'
#' In environmental criminology the report is the burden side of green
#' criminology: the harm a pollutant does to a population and how it falls
#' by income, from monitoring data an agency or a court can be given.
#'
#' @details
#' The assumption log checks that the exposure mean is finite and
#' non-negative and at most 10,000 micrograms per cubic metre, the prevalence
#' is in \eqn{[0, 1]}, the baseline rate (converted to per 100,000) is zero or
#' between 1 and 10,000 per 100,000 a year (from 1e-6 with \code{rare}, up to
#' 100,000 with \code{hirate}), the population is positive and at most 10
#' billion, the reference is finite and non-negative, and the pollutant has a
#' concentration-response function. Any failure gives status
#' \code{"assumption_failure"} and skips the stages.
#'
#' With data (a CSV or the demo), the \code{"exposed"} are the rows above the
#' reference: their mean is \code{expomean} and their share \code{expoprev}.
#' Empty cells and missing-value tokens (\code{NA}, \code{n/a}, \code{null},
#' \code{-} ...) are dropped with a message; text that is not a number (a
#' detection-limit flag such as \code{"<2"}) and negative values stop the run
#' with an error report. A gzip, bzip2 or xz file is decompressed (to at most
#' 1 GiB, and at most 200 times its size or 64 MiB, whichever is more); UTF-16
#' and Latin-1 files are decoded. An \code{income} column, when present and
#' numeric, gives the equity stage (\code{\link{expconcidx}}).
#'
#' The stages: \code{\link{crfpm25}} or \code{\link{crfno2}};
#' \code{\link{attrfrac}} with the curve's interval; \code{\link{mortdisp}}
#' for removing the excess exposure over the exposed population, with
#' \eqn{\beta = \log(RR) / (\bar x - x_{ref})}; and \code{\link{polbrdn}}.
#' The demo data are 1,000 simulated units with lognormal exposures at the
#' scale of Toronto NO2 or PM2.5 levels and five income groups, drawn from R's
#' generator with a fixed seed (R's own random stream is restored).
#'
#' @param poll \code{"no2"} or \code{"pm25"} (also \code{"PM2.5"}), in either
#'   case: one name (another name gives a report whose assumption on the
#'   pollutant fails; anything but one name is an R error).
#' @param y Outcome key for the concentration-response function; the
#'   mortality outcomes count deaths, \code{"chiast"} incident cases.
#' @param reg,years Labels for the report.
#' @param demo Use simulated demo data (\code{TRUE} or \code{FALSE}).
#' @param expcsv Path of a CSV with an \code{exposure} column in micrograms
#'   per cubic metre, or a NAPS pull with a \code{value} column and a
#'   \code{unit} column (required with \code{value}): each row is converted by
#'   its own unit, NO2 in ppb at 1.88 micrograms per cubic metre per ppb (the
#'   WHO 2021 conversion at 25 C and 1 atm), and any other unit is an error.
#'   With a CSV or the demo data the exposed are the rows above \code{ref}:
#'   \code{expoprev} is their share and \code{expomean} their mean, so Levin's
#'   formula is not diluted twice.
#' @param expomean,expoprev Scalar inputs used when neither \code{demo} nor
#'   \code{expcsv} is given: the mean exposure among the exposed and the share
#'   of the population exposed. They are required together: a prevalence of
#'   zero silently zeroes the burden.
#' @param ref Counterfactual reference concentration; \code{NULL} takes the
#'   pollutant's default (PM2.5 5.8, NO2 10; see \code{\link{crfpm25}}).
#' @param brate Baseline outcome rate per 100,000 per year (unlike
#'   \code{\link{polbrdn}}, which takes a rate per person-year).
#' @param basrateper The denominator \code{brate} is quoted per: 100,000 (the
#'   default), or 1 for a rate per person-year such as a crude death rate of
#'   0.012.
#' @param pop Population at risk (persons).
#' @param rare \code{TRUE} when the outcome really is rarer than one case per
#'   100,000 per year. Without it a rate that converts to less than that is an
#'   assumption failure, because it is also what a rate per person-year passed
#'   as per 100,000 looks like (0.008 per person-year is 800 per 100,000; as
#'   0.008 per 100,000 the burden is 100,000 times too small).
#' @param hirate \code{TRUE} when the outcome really is more common than
#'   10,000 per 100,000 a year (no population's all-cause death rate is);
#'   without it such a rate is an assumption failure.
#' @param report A result of \code{vrfpol()} (or a list with the same parts);
#'   anything else is refused with an R error.
#' @param x A result of \code{vrfpol()}.
#' @inheritParams freqtbl
#' @return \code{vrfpol()}: the report as a list of class \code{vrfpol}; its
#'   \code{status} is \code{"ok"}, \code{"assumption_failure"} or
#'   \code{"error"}, and the attribute \code{exitstat} carries a command-line
#'   exit code (0, 1 or 2). An argument the pipeline cannot use (no exposure
#'   data at all, a vector where a single number is expected, an unknown
#'   outcome, a missing reference, a directory as \code{expcsv}) is an
#'   \code{"error"} report with exit status 2, not an R error (\code{poll}
#'   aside, see the description); an exposure
#'   at or below the reference is a valid run with zero attributable cases.
#'   On \code{"ok"} the report also holds \code{inputs} (the resolved exposure
#'   mean and prevalence, the rate per 100,000, the population, the reference
#'   and the rows dropped), \code{assumptions} (a list of \code{assumption},
#'   \code{ok}, \code{note}) and \code{pipeline} (\code{crf}, \code{paf},
#'   \code{pafci}, \code{displaced}, \code{burden} and \code{equity}, which is
#'   \code{NULL} without income data). \code{polrpttxt()}: a single character
#'   string with embedded newlines, ending in a \code{STATUS:} line (or
#'   \code{"ERROR:"} and the message for an error report).
#' @seealso \code{\link{polbrdn}} and the other stage functions to run one
#'   step on its own.
#' @examples
#' r <- vrfpol("no2", demo = TRUE)
#' r$status
#' r$pipeline$paf
#' cat(polrpttxt(r))
#'
#' # Scalar inputs: a city of 2.8 million at a PM2.5 mean of 9.
#' s <- vrfpol("pm25", expomean = 9, expoprev = 1, brate = 750, pop = 2.8e6,
#'             reg = "Toronto")
#' s$pipeline$burden$attrcases
#' attr(s, "exitstat")
#'
#' # A CSV of simulated measurements with an income column.
#' f <- tempfile(fileext = ".csv")
#' utils::write.csv(data.frame(exposure = c(12, 18, 25, 31, 9), income = 1:5),
#'                  f, row.names = FALSE)
#' vrfpol("no2", expcsv = f)$pipeline$equity$concidx
#'
#' # A rate per person-year passed as per 100,000 fails an assumption ...
#' bad <- vrfpol("no2", expomean = 20, expoprev = 1, brate = 0.008)
#' bad$status
#' # ... and an argument it cannot use is an error report, not an R error.
#' vrfpol("no2")$error
#' @export
vrfpol <- function(
  poll,
  y = "allcause",
  reg = NULL,
  years = NULL,
  demo = FALSE,
  expcsv = NULL,
  expomean = NULL,
  expoprev = NULL,
  ref = NULL,
  brate = 500,
  basrateper = 1e5,
  pop = 1e6,
  rare = FALSE,
  hirate = FALSE
) {
  f <- "vrfpol"
  poll <- .vpkey(.cktxt(poll, "poll", f), f)
  if (!is.logical(demo) || length(demo) != 1L || is.na(demo)) {
    return(.vpfail("demo must be TRUE or FALSE"))
  }
  inp <- .vpinputs(
    demo, expcsv, poll, ref, expomean, expoprev,
    list(brate = brate, basrateper = basrateper, pop = pop), rare, hirate, y
  )
  if (is.character(inp)) {
    return(.vpfail(inp))
  }
  brate <- as.numeric(brate) * 1e5 / basrateper # per 100,000 from here on
  .vprun(poll, y, reg, years, inp$ex, brate, as.numeric(pop), inp$ref, rare,
         hirate)
}

#' @rdname vrfpol
#' @export
polrpttxt <- function(report) {
  if (
    !inherits(report, "vrfpol") &&
      !(is.list(report) && is.character(report$status))
  ) {
    .ckstop("polrpttxt", "`report` must be a report from vrfpol()")
  }
  if (identical(report$status, "error")) {
    return(paste0("ERROR: ", report$error, "\n"))
  }
  if (!.vpok(report)) {
    .ckstop(
      "polrpttxt",
      paste0(
        "`report` is not a complete vrfpol() report (its inputs, ",
        "assumptions or pipeline are missing)"
      )
    )
  }
  lines <- c(.vphead(report), .vpasm(report$assumptions))
  if (identical(report$status, "assumption_failure")) {
    return(paste0(
      paste(lines, collapse = "\n"),
      "\n\nSTATUS: assumption_failure (pipeline skipped)\n"
    ))
  }
  lines <- c(lines, .vppipe(report$pipeline))
  paste0(paste(lines, collapse = "\n"), "\n\nSTATUS: ok\n")
}

#' @rdname vrfpol
#' @export
print.vrfpol <- function(x, ...) {
  cat(polrpttxt(x))
  invisible(x)
}

#' polrpttxt(): does the report have its inputs, assumptions and pipeline?
#' @noRd
.vpok <- function(report) {
  is.list(report$inputs) &&
    is.list(report$assumptions) &&
    all(vapply(
      report$assumptions,
      function(a) is.list(a) && is.logical(a$ok),
      logical(1)
    )) &&
    (!identical(report$status, "ok") || is.list(report$pipeline))
}

#' `b` when `a` is NULL
#' @noRd
.vpor <- function(a, b) if (is.null(a)) b else a

#' The pollutant `poll` as the key the pollution check uses.
#' @noRd
.vpkey <- function(poll, fn) {
  if (!is.character(poll) || length(poll) != 1L || is.na(poll)) {
    .ckstop(fn, "`poll` must be one pollutant name (\"PM2.5\" or \"NO2\")")
  }
  pol <- .vppol(poll)
  if (is.na(pol)) {
    tolower(as.character(poll))
  } else if (pol == "NO2") {
    "no2"
  } else {
    "pm25"
  }
}

#' The pollutant name as its standard name, NA when it is not recognised.
#' @noRd
.vppol <- function(p) {
  if (!is.character(p) || length(p) != 1L || is.na(p)) {
    return(NA_character_)
  }
  key <- gsub("[^a-z0-9]", "", tolower(p))
  if (key %in% c("pm25", "pm2.5")) {
    return("PM2.5")
  }
  if (identical(key, "no2")) {
    return("NO2")
  }
  NA_character_
}

#' The result of a pollution check that could not run: the error and exit
#'   status 2.
#' @noRd
.vpfail <- function(msg) {
  structure(
    list(status = "error", error = msg),
    exitstat = 2L,
    class = c("vrfpol", "list")
  )
}

#' One row of the table of assumptions: its name, whether it holds, a note.
#' @noRd
.vprow <- function(name, ok, note) {
  list(assumption = name, ok = isTRUE(ok), note = note)
}

#' Whether `x` is one finite number.
#' @noRd
.vpfin <- function(x) length(x) == 1L && is.finite(x)

#' The assumptions of the pollution check, each with whether it holds.
#' @noRd
.vpasmlist <- function(poll, em, ep, br, pop, ref, rare, hirate) {
  fin <- .vpfin
  list(
    .vprow(
      "exposure finite and non-negative",
      fin(em) && em >= 0,
      .vpexpnote(em, ref)
    ),
    .vprow(
      "exposure within the concentration-response range (at most 10,000 ug/m3)",
      fin(em) && em <= 1e4,
      sprintf("exposure_mean=%s", format(em))
    ),
    .vprow(
      "prevalence in [0,1]",
      fin(ep) && ep >= 0 && ep <= 1,
      sprintf("exposure_prevalence=%s", format(ep))
    ),
    .vprate(br, rare, hirate),
    .vprow(
      "population finite, positive and at most 10 billion",
      fin(pop) && pop > 0 && pop <= 1e10,
      if (fin(pop) && pop > 1e10) {
        sprintf(
          "population=%s is more people than live on Earth",
          format(pop, big.mark = ",", scientific = FALSE)
        )
      } else {
        sprintf(
          "population=%s",
          format(pop, big.mark = ",", scientific = FALSE)
        )
      }
    ),
    .vprow(
      "reference finite and non-negative",
      fin(ref) && ref >= 0,
      sprintf("reference=%s", format(ref))
    ),
    .vprow(
      "pollutant supported by the CRFs",
      !is.na(.vppol(poll)),
      paste(
        "Current CRFs: NO2 (log-linear), PM2.5 (log-linear all-cause;",
        "Burnett IER for IHD and stroke).",
        "Other pollutants reject."
      )
    )
  )
}

#' The note on the exposure mean.
#' @noRd
.vpexpnote <- function(em, ref) {
  if (!.vpfin(em) || em < 0) {
    sprintf(
      "expomean=%s is not a finite non-negative concentration",
      format(em)
    )
  } else if (isTRUE(em <= ref)) {
    sprintf(
      paste0(
        "mean %s is at or below the counterfactual floor %s: no ",
        "excess risk, PAF and cases are 0."
      ),
      format(em, digits = 10),
      format(ref)
    )
  } else {
    sprintf("mean %s vs ref %s", format(em, digits = 10), format(ref))
  }
}

#' The row of assumptions for the baseline rate.
#' @noRd
.vprate <- function(br, rare, hirate) {
  top <- if (isTRUE(hirate)) 1e5 else 1e4
  fin <- .vpfin(br)
  .vprow(
    sprintf(
      "brate finite: zero, or between %s and %s per 100k%s",
      if (isTRUE(rare)) "1e-6" else "1",
      format(top, big.mark = ",", scientific = FALSE),
      if (isTRUE(rare)) " (a rare outcome)" else ""
    ),
    fin && br >= 0 && br <= top &&
      (br == 0 || br >= 1 || (isTRUE(rare) && br >= 1e-6)),
    .vpratenote(br, fin, rare, top)
  )
}

#' The note on the baseline rate.
#' @noRd
.vpratenote <- function(br, fin, rare, top) {
  if (fin && br < 1e-6 && br > 0) {
    sprintf(
      paste(
        "brate=%s per 100k is below one case per 100 billion",
        "person-years: not a",
        "usable rate (check basrateper)"
      ),
      format(br)
    )
  } else if (fin && br < 1 && br > 0 && !isTRUE(rare)) {
    sprintf(
      paste(
        "brate=%s per 100k is below one case per 100,000. If it is",
        "a rate per person-year",
        "(0.012 is 1,200 per 100k), pass basrateper = 1; if the",
        "outcome really is this",
        "rare, say so with rare = TRUE"
      ),
      format(br)
    )
  } else if (fin && br > top && br <= 1e5) {
    sprintf(
      paste(
        "brate=%s per 100k a year is above a whole population's",
        "all-cause death rate",
        "(an age-specific rate for the oldest ages can be this",
        "high): if the outcome really is",
        "this common, say so with hirate = TRUE"
      ),
      format(br, big.mark = ",", scientific = FALSE)
    )
  } else if (fin && br == 0) {
    "brate=0: a measured zero rate, so no cases"
  } else {
    sprintf("brate=%s per 100k per year", format(br))
  }
}

#' The simulated demo data: lognormal exposures and five income groups,
#' drawn from R's generator with a fixed seed.
#' @noRd
.vpdemo <- function(poll) {
  .clseed(20260417L, {
    n <- 1000L
    ml <- switch(poll, no2 = log(22), pm25 = log(9), log(15))
    sl <- if (poll == "pm25") 0.30 else 0.35
    exposure <- stats::rlnorm(n, ml, sl)
    income <- sample(1:5, n, replace = TRUE,
                     prob = c(0.18, 0.22, 0.22, 0.20, 0.18))
    mult <- c(1.25, 1.12, 1.00, 0.93, 0.85)[income]
    data.frame(exposure = exposure * mult, income = income)
  })
}

#' The mean and share of the rows above the reference (the exposed).
#' @noRd
.vpexposed <- function(x, ref) {
  x <- x[!is.na(x)]
  above <- if (is.finite(ref)) x[x > ref] else x
  list(
    mean = if (length(above)) {
      mean(above)
    } else if (length(x)) {
      mean(x)
    } else {
      NA_real_
    },
    prevalence = if (length(x)) length(above) / length(x) else NA_real_
  )
}

#' The exposure inputs of the pollution check: the demonstration data, a CSV
#'   file, or the numbers given.
#' @noRd
.vpsource <- function(demo, expcsv, poll, ref, em, ep, noscalar) {
  if (isTRUE(demo)) {
    df <- .vpdemo(poll)
    es <- .vpexposed(df$exposure, ref)
    return(list(
      expomean = es$mean, expoprev = es$prevalence, equity_df = df,
      data_source = "demo (synthetic)", dropped = 0L
    ))
  }
  if (!is.null(expcsv)) {
    rd <- .vpread(expcsv, poll)
    if (is.character(rd)) {
      return(rd)
    }
    es <- .vpexposed(rd$df$exposure, ref)
    return(list(
      expomean = es$mean, expoprev = es$prevalence, equity_df = rd$equity_df,
      data_source = expcsv, dropped = rd$dropped
    ))
  }
  if (noscalar) {
    return(paste(
      "no exposure data: pass expcsv, or expomean and expoprev,",
      "or demo = TRUE"
    ))
  }
  list(
    expomean = as.numeric(em), expoprev = as.numeric(ep), equity_df = NULL,
    data_source = "scalar arguments", dropped = 0L
  )
}

#' The reason the inputs of the pollution check cannot be used, NULL when they
#'   can.
#' @noRd
.vpcheck <- function(demo, expcsv, sc, rare, hirate, y) {
  m <- .vpsrc(demo, expcsv, sc)
  if (is.null(m)) .vpnums(sc, rare, hirate, y) else m
}

#' .vpcheck(): one source of exposure, and a population
#' @noRd
.vpsrc <- function(demo, expcsv, sc) {
  noscalar <- is.null(sc$expomean) && is.null(sc$expoprev)
  if (isTRUE(demo) && !is.null(expcsv)) {
    return(paste(
      "expcsv and demo = TRUE are two sources of exposure: pass one",
      "(the demo was used silently before)"
    ))
  }
  if ((isTRUE(demo) || !is.null(expcsv)) && !noscalar) {
    return(paste(
      "expomean / expoprev cannot be combined with expcsv or demo:",
      "pass one source of exposure"
    ))
  }
  if (xor(is.null(sc$expomean), is.null(sc$expoprev))) {
    return(paste0(
      "expomean and expoprev go together: pass both (a zero is a ",
      "measurement), or neither"
    ))
  }
  if (is.null(sc$pop)) {
    return("population must be a single number (persons at risk), not NULL")
  }
  NULL
}

#' .vpcheck(): the numbers and flags, each of the right kind and length
#' @noRd
.vpnums <- function(sc, rare, hirate, y) {
  for (nm in c("expomean", "expoprev", "brate", "basrateper", "pop", "ref")) {
    v <- sc[[nm]]
    if (!is.null(v) && !is.numeric(v) && !(is.logical(v) && all(is.na(v)))) {
      return(sprintf("%s must be a single number, not %s", nm, class(v)[1L]))
    }
  }
  if (is.null(sc$brate)) {
    return("brate must be a single number, not NULL")
  }
  if (!is.logical(rare) || length(rare) != 1L || is.na(rare)) {
    return("rare must be TRUE or FALSE")
  }
  if (!is.logical(hirate) || length(hirate) != 1L || is.na(hirate)) {
    return("hirate must be TRUE or FALSE")
  }
  bp <- sc$basrateper
  if (!is.numeric(bp) || length(bp) != 1L || !is.finite(bp) || bp <= 0) {
    return(paste(
      "basrateper must be a single positive number",
      "(100000 for a rate per 100,000, 1 per person-year)"
    ))
  }
  for (nm in c("ref", "expomean", "expoprev", "brate", "pop")) {
    v <- sc[[nm]]
    if (!is.null(v) && length(v) != 1L) {
      return(sprintf("%s must be a single number, not %d values", nm,
                     length(v)))
    }
  }
  if (!is.character(y) || length(y) != 1L || is.na(y)) {
    return("outcome must be a single outcome name")
  }
  NULL
}

#' The inputs of the pollution check gathered and checked.
#' @noRd
.vpinputs <- function(demo, expcsv, poll, ref, em, ep, rates, rare, hirate,
                      y) {
  if (is.null(ref)) {
    ref <- if (identical(poll, "no2")) 10 else 5.8
  }
  noscalar <- is.null(em) && is.null(ep)
  why <- .vpcheck(
    demo, expcsv,
    list(expomean = em, expoprev = ep, brate = rates$brate,
         basrateper = rates$basrateper, pop = rates$pop, ref = ref),
    rare, hirate, y
  )
  if (!is.null(why)) {
    return(why)
  }
  ref <- suppressWarnings(as.numeric(ref))
  ex <- .vpsource(demo, expcsv, poll, ref, em, ep, noscalar)
  if (is.character(ex)) {
    return(ex)
  }
  list(ex = ex, ref = ref)
}

#' Runs the pollution check and builds its report.
#' @noRd
.vprun <- function(poll, y, reg, years, ex, brate, pop, ref, rare, hirate) {
  asm <- .vpasmlist(poll, ex$expomean, ex$expoprev, brate, pop, ref, rare,
                    hirate)
  report <- list(
    command = "verify-pollution",
    pollutant = poll,
    outcome = y,
    region = reg,
    years = years,
    dsou = ex$data_source,
    inputs = list(
      expomean = ex$expomean,
      expoprev = ex$expoprev,
      baserateper100k = brate,
      population = pop,
      refconc = ref,
      drorows = ex$dropped
    )
  )
  report$assumptions <- asm
  if (!all(vapply(asm, function(a) isTRUE(a$ok), logical(1)))) {
    report$status <- "assumption_failure"
    report$pipeline <- list(skipped = TRUE, reason = "assumptions")
    return(structure(report, exitstat = 1L, class = c("vrfpol", "list")))
  }
  pipe <- .vppipeline(poll, y, ex$expomean, ex$expoprev, brate, pop, ref,
                      rare, hirate, ex$equity_df)
  if (is.character(pipe)) {
    return(.vpfail(pipe))
  }
  report$status <- "ok"
  report$pipeline <- pipe
  structure(report, exitstat = 0L, class = c("vrfpol", "list"))
}

#' The computations of the pollution check, in order; the message of the error
#'   when one fails.
#' @noRd
.vppipeline <- function(poll, y, em, ep, brate, pop, ref, rare, hirate, eq) {
  crf <- .vpcrf(poll, em, y, ref)
  if (inherits(crf, "error")) {
    return(conditionMessage(crf))
  }
  paf <- tryCatch(
    attrfrac(
      crf$rr, ep,
      rrlo = if (is.finite(crf$rrlo)) crf$rrlo,
      rrup = if (is.finite(crf$rrup)) crf$rrup
    ),
    error = function(e) e
  )
  if (inherits(paf, "error")) {
    return(conditionMessage(paf))
  }
  bpp <- brate / 1e5
  del <- max(0, em - ref)
  bpu <- log(crf$rr) / max(del, 1e-9)
  st <- tryCatch(
    list(
      displaced = mortdisp(del, pop * ep, bpp, bpu, rare = rare,
                           hirate = hirate),
      burden = polbrdn(em, ep, bpp, pop,
                       poll = if (poll == "pm25") "PM2.5" else "NO2", y = y,
                       refconc = ref, rare = rare, hirate = hirate),
      equity = .vpequity(eq)
    ),
    error = function(e) e
  )
  if (inherits(st, "error")) {
    return(conditionMessage(st))
  }
  .vpout(crf, paf, st, del, bpu, y)
}

#' .vppipeline(): the concentration-response function, or its error
#' @noRd
.vpcrf <- function(poll, em, y, ref) {
  tryCatch(
    if (poll == "no2") {
      crfno2(em, y = y, refconc = ref)
    } else {
      crfpm25(em, y = y, refconc = ref)
    },
    error = function(e) e
  )
}

#' .vppipeline(): the stages' results as the report's pipeline
#' @noRd
.vpout <- function(crf, paf, st, del, bpu, y) {
  list(
    crf = unclass(crf),
    paf = as.numeric(paf),
    pafci = attr(paf, "ci"),
    displaced = list(
      deadisp = st$displaced,
      expdel = del,
      betaperunit = bpu,
      unit = if (identical(y, "chiast")) "incident cases" else "deaths"
    ),
    burden = unclass(st$burden),
    equity = if (is.null(st$equity)) NULL else unclass(st$equity)
  )
}

#' The equity part of the pollution check, NULL without an income column.
#' @noRd
.vpequity <- function(eq) {
  if (is.null(eq) || !"income" %in% names(eq)) {
    return(NULL)
  }
  complete <- sum(
    is.finite(eq$exposure) & is.finite(suppressWarnings(as.numeric(eq$income)))
  )
  if (complete < 2L || !any(eq$exposure > 0, na.rm = TRUE)) {
    return(NULL)
  }
  expconcidx(eq, "exposure", "income")
}

# ---- reading the exposure file, with base R --------------------------------

#' The exposure CSV read and checked; a message when it cannot be used.
#' @noRd
.vpread <- function(expcsv, poll) {
  df <- .vpload(expcsv)
  if (is.character(df)) {
    return(df)
  }
  raw <- NULL
  naps <- FALSE
  if (!"exposure" %in% names(df) && "value" %in% names(df)) {
    naps <- TRUE
    nv <- .vpnaps(df, poll)
    if (is.character(nv)) {
      return(nv)
    }
    raw <- nv$raw
    df$exposure <- nv$vals
  }
  if (!"exposure" %in% names(df)) {
    return(paste0(
      "CSV missing 'exposure' column (ug/m3); a NAPS pull's 'value' ",
      "column also works."
    ))
  }
  if (!naps && is.character(df$exposure)) {
    raw <- trimws(df$exposure)
  }
  df$exposure <- suppressWarnings(as.numeric(df$exposure))
  why <- .vptext(df$exposure, raw)
  if (!is.null(why)) {
    return(why)
  }
  .vprows(df)
}

#' The exposure CSV loaded; a message when the path is not a usable file.
#' @noRd
.vpload <- function(expcsv) {
  if (!is.character(expcsv) || length(expcsv) != 1L || is.na(expcsv)) {
    return("expcsv must be a single file path")
  }
  if (!file.exists(expcsv)) {
    return(sprintf("exposure CSV not found: %s", expcsv))
  }
  if (dir.exists(expcsv)) {
    return(sprintf("exposure CSV is a directory, not a file: %s", expcsv))
  }
  empty <- sprintf(
    paste0(
      "%s is empty (expected a column exposure in ug/m3, or a NAPS ",
      "pull's value + unit)"
    ),
    expcsv
  )
  if (isTRUE(file.size(expcsv) == 0)) {
    return(empty)
  }
  df <- tryCatch(.vpcsv(expcsv), error = function(e) e)
  if (inherits(df, "error")) {
    return(sprintf("%s: %s", expcsv, conditionMessage(df)))
  }
  if (!nrow(df)) {
    return(empty)
  }
  df
}

#' The exposure column of a file in NAPS form (`value` with `unit`); a message
#'   when it cannot be read.
#' @noRd
.vpnaps <- function(df, poll) {
  if (!"unit" %in% names(df)) {
    return(paste0(
      "a 'value' column needs a 'unit' column (NAPS pulls carry ",
      "one; NO2 is in ppb)"
    ))
  }
  raw <- if (is.character(df$value)) trimws(df$value)
  vals <- suppressWarnings(as.numeric(df$value))
  unit <- tolower(trimws(as.character(df$unit)))
  ugm3 <- c("ug/m3", "\u00b5g/m3", "\u00b5g/m\u00b3", "ug/m\u00b3")
  ppb <- unit %in% c("ppb", "ppbv")
  bad <- !is.na(vals) & !(unit %in% ugm3) & !(ppb & poll == "no2")
  if (any(bad)) {
    return(sprintf(
      "exposure unit %s is not ug/m3 for %s",
      toString(sort(unique(unit[bad]))),
      poll
    ))
  }
  if (any(ppb & !is.na(vals))) {
    # ug/m3 per ppb of NO2 at 25 C and 1 atm (WHO 2021)
    vals[ppb] <- vals[ppb] * 1.88
    message(sprintf(
      "note: %d NO2 row(s) in ppb converted to ug/m3 (x 1.88)",
      sum(ppb & !is.na(vals))
    ))
  }
  list(vals = vals, raw = raw)
}

#' A note on exposure cells that hold text, neither a number nor a mark of a
#'   missing value.
#' @noRd
.vptext <- function(exposure, raw) {
  if (is.null(raw)) {
    return(NULL)
  }
  tok <- is.na(raw) |
    tolower(raw) %in% c("", "na", "n/a", "nan", "null", "-", ".", "missing")
  txt <- is.na(exposure) & !tok
  if (!any(txt)) {
    return(NULL)
  }
  ex <- utils::head(unique(raw[txt]), 3L)
  sprintf(
    paste(
      "%d exposure value(s) are text, not numbers (%s): a censored",
      "value such as \"<2\"",
      "is a measurement, and dropping it would bias the mean and",
      "the prevalence; replace",
      "them (for example with half the detection limit) or remove",
      "those rows"
    ),
    sum(txt),
    toString(sprintf("\"%s\"", ex))
  )
}

#' Drops the rows with no exposure value, with a note.
#' @noRd
.vprows <- function(df) {
  dropped <- sum(is.na(df$exposure))
  if (dropped > 0L) {
    message(sprintf(
      "note: %d row(s) have no exposure value and were dropped",
      dropped
    ))
  }
  if (any(df$exposure < 0, na.rm = TRUE)) {
    return(sprintf(
      "%d exposure value(s) are negative: a concentration cannot be negative",
      sum(df$exposure < 0, na.rm = TRUE)
    ))
  }
  if ("income" %in% names(df) && all(is.na(df$income))) {
    message("note: the income column is empty: no equity analysis")
    df$income <- NULL
  }
  if ("income" %in% names(df) && !is.numeric(df$income)) {
    return(paste0(
      "the income column must be numeric (a quintile or an amount); ",
      "text cannot be ranked"
    ))
  }
  list(df = df, dropped = dropped, equity_df = if ("income" %in% names(df)) df)
}

#' Reads the CSV at `path` after checking what kind of file it is.
#' @noRd
.vpcsv <- function(path) {
  head <- readBin(path, "raw", 6L)
  comp <- length(head) >= 2L &&
    (identical(head[1:2], as.raw(c(0x1f, 0x8b))) ||
      (length(head) >= 3L && identical(head[1:3], charToRaw("BZh"))) ||
      (length(head) >= 6L &&
        identical(head[1:6], as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)))))
  raw <- if (comp) {
    .vpunzip(path, head)
  } else {
    readBin(path, "raw", file.size(path))
  }
  if (!length(raw)) {
    return(data.frame())
  }
  enc <- .vpenc(raw, path)
  dec <- .vpdec(enc$raw, enc$utf16, path)
  txt <- dec$txt
  note <- .vpor(dec$note, enc$note)
  # a copy ends with a newline: an incomplete last line is then a quote
  # left open (.vpquiet())
  changed <- comp || !is.null(enc$utf16) || dec$transformed ||
    !endsWith(txt, "\n")
  rm(raw, enc, dec)
  if (!is.null(note)) {
    message(note)
  }
  Encoding(txt) <- "UTF-8"
  if (!grepl("[^[:space:]]", txt, perl = TRUE)) {
    return(data.frame())
  }
  src <- path
  if (changed) {
    src <- tempfile(fileext = ".csv")
    on.exit(unlink(src), add = TRUE)
    writeLines(txt, src, useBytes = TRUE)
  }
  rm(txt)
  .vpshape(src)
  .vpparse(src)
}

#' A compressed file decompressed, to a bounded size.
#' @noRd
.vpunzip <- function(path, head) {
  lim <- min(2^30, max(64 * 2^20, 200 * file.size(path)), 2^31 - 2)
  bz <- length(head) >= 3L && identical(head[1:3], charToRaw("BZh"))
  if (bz && !.vpbz2end(path)) {
    stop(
      "is truncated: the bzip2 stream has no end-of-stream marker",
      call. = FALSE
    )
  }
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  chunks <- list()
  total <- 0
  withCallingHandlers(
    repeat {
      b <- readBin(con, "raw", 4194304L)
      if (!length(b)) break
      total <- total + length(b)
      if (total > lim) {
        stop("decompresses to more than the size limit", call. = FALSE)
      }
      chunks[[length(chunks) + 1L]] <- b
    },
    warning = function(w) {
      stop(
        sprintf(
          "is truncated or damaged (the decoder said: %s)",
          conditionMessage(w)
        ),
        call. = FALSE
      )
    }
  )
  out <- do.call(c, c(list(raw()), chunks))
  if (identical(head[1:2], as.raw(c(0x1f, 0x8b)))) .vpgzend(path, out)
  out
}

#' Whether the data inflated from a gzip file match the CRC-32 and length
#' that end its last member (RFC 1952). R's reader checks them only when the
#' member ends, not when the file stops inside the data. The last member's
#' data are the last ISIZE bytes inflated.
#' @noRd
.vpgzend <- function(path, out) {
  n <- file.size(path)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, max(0, n - 8))
  if (!.Call(C_clab_r_gzend, out, readBin(con, "raw", 8L))) {
    stop(
      "is truncated or damaged: its data do not match the length and ",
      "CRC-32 at its end",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Whether a bzip2 file ends with its end-of-stream marker.
#' @noRd
.vpbz2end <- function(path) {
  n <- file.size(path)
  if (is.na(n) || n < 14) {
    return(FALSE)
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, n - 11)
  .Call(C_clab_r_bz2end, readBin(con, "raw", 11L))
}

#' The text encoding of the bytes of a file; NUL bytes that are not UTF-16 are
#'   refused.
#' @noRd
.vpenc <- function(raw, path) {
  r <- .Call(C_clab_r_txtenc, raw)
  if (r[1L] < 0L) {
    stop(
      "contains NUL bytes and is not UTF-16 text: not a text CSV",
      call. = FALSE
    )
  }
  if (r[1L] == 0L) {
    return(list(utf16 = NULL, raw = raw, note = NULL))
  }
  enc <- c("UTF-32LE", "UTF-32BE", "UTF-16LE", "UTF-16BE")[r[1L]]
  if (r[2L] > 0L) {
    return(list(utf16 = enc, raw = raw[-seq_len(r[2L])], note = NULL))
  }
  list(
    utf16 = enc,
    raw = raw,
    note = sprintf(
      "note: %s has no byte-order mark; read as %s",
      basename(path),
      enc
    )
  )
}

#' The bytes of a file decoded to UTF-8 text.
#' @noRd
.vpdec <- function(raw, utf16, path) {
  if (!is.null(utf16)) {
    txt <- iconv(list(raw), from = utf16, to = "UTF-8")[[1L]]
    if (is.null(txt) || is.na(txt)) {
      stop(sprintf("cannot be decoded as %s", utf16), call. = FALSE)
    }
    return(list(txt = txt, note = NULL, transformed = FALSE))
  }
  transformed <- FALSE
  if (length(raw) >= 3L && identical(raw[1:3], as.raw(c(0xef, 0xbb, 0xbf)))) {
    raw <- raw[-(1:3)]
    transformed <- TRUE
  }
  txt <- rawToChar(raw)
  if (validUTF8(txt)) {
    return(list(txt = txt, note = NULL, transformed = transformed))
  }
  list(
    txt = iconv(txt, from = "latin1", to = "UTF-8"),
    note = sprintf(
      "note: %s is not UTF-8; read as Latin-1 (ISO-8859-1)",
      basename(path)
    ),
    transformed = TRUE
  )
}

#' Checks the header row of the CSV `src`.
#' @noRd
.vpshape <- function(src) {
  hdr <- names(.vpquiet(
    utils::read.csv(src, nrows = 1L, check.names = FALSE, encoding = "UTF-8")
  ))
  dup <- unique(hdr[duplicated(tolower(trimws(hdr)))])
  if (length(dup)) {
    stop(
      sprintf(
        "has the column name %s more than once: say which column is meant",
        toString(sprintf("'%s'", dup))
      ),
      call. = FALSE
    )
  }
  nf <- utils::count.fields(src, sep = ",", quote = "\"",
                            blank.lines.skip = TRUE, comment.char = "")
  wrong <- which(!is.na(nf) & nf != length(hdr))
  if (length(wrong)) {
    stop(
      sprintf(
        paste(
          "line %d has %d fields where the header has %d:",
          "no rows are used from a file whose rows do not line up"
        ),
        wrong[[1L]],
        nf[[wrong[[1L]]]],
        length(hdr)
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' The CSV `src` read as a data frame.
#' @noRd
.vpparse <- function(src) {
  .vpquiet(utils::read.csv(src, stringsAsFactors = FALSE, encoding = "UTF-8"))
}

#' A read of the CSV in which any warning refuses the file. The file read
#' ends with a newline (.vpcsv()), so "incomplete final line" means the
#' header read reached the end inside a quote.
#' @noRd
.vpquiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (grepl("incomplete final line", m, fixed = TRUE)) {
      m <- "EOF within quoted string"
    }
    stop(
      sprintf(
        paste0(
          "could not be read cleanly (%s): no rows are used from a file ",
          "the parser stopped part-way through"
        ),
        m
      ),
      call. = FALSE
    )
  })
}

# ---- the text report ---------------------------------------------------------

#' The heading lines of the printed pollution report.
#' @noRd
.vphead <- function(report) {
  bar <- strrep("=", 66)
  lines <- c(
    bar,
    sprintf("  %s -- %s -> %s", report$command, toupper(report$pollutant),
            report$outcome)
  )
  if (!is.null(report$region)) {
    lines <- c(lines, sprintf("  region: %s", report$region))
  }
  if (!is.null(report$years)) {
    lines <- c(lines, sprintf("  years:  %s", report$years))
  }
  inp <- report$inputs
  c(
    lines,
    sprintf("  data:   %s", report$dsou),
    bar,
    "",
    "Inputs",
    sprintf("  exposure mean:      %.3f", inp$expomean),
    sprintf("  exposure prevalence:%.3f", inp$expoprev),
    sprintf("  baseline rate/100k: %.2f", inp$baserateper100k),
    sprintf(
      "  population:         %s",
      format(inp$population, big.mark = ",", scientific = FALSE)
    ),
    sprintf("  reference conc:     %s", format(inp$refconc)),
    if (isTRUE(inp$drorows > 0)) {
      sprintf("  rows dropped:       %d (no exposure value)", inp$drorows)
    },
    "",
    "Assumption log"
  )
}

#' The assumptions as printed lines, PASS or FAIL.
#' @noRd
.vpasm <- function(assumptions) {
  vapply(
    assumptions,
    function(a) {
      sprintf(
        "  [%s] %s -- %s",
        if (isTRUE(a$ok)) "PASS" else if (isFALSE(a$ok)) "FAIL" else "NA",
        a$assumption,
        a$note
      )
    },
    ""
  )
}

#' The results of the pollution check as printed lines.
#' @noRd
.vppipe <- function(p) {
  ci <- if (is.null(p$pafci)) {
    ""
  } else {
    sprintf(
      "  (95%% CI %.4f to %.4f, from the risk's limits)",
      p$pafci[["lower"]],
      p$pafci[["upper"]]
    )
  }
  lines <- c(
    "",
    "Concentration-response",
    sprintf("  RR:       %.4f", p$crf$rr),
    sprintf("  source:   %s", p$crf$citation),
    "",
    sprintf("Attributable fraction (PAF): %.4f%s", p$paf, ci),
    "",
    if (identical(p$displaced$unit, "deaths")) {
      "Mortality displaced"
    } else {
      "Cases displaced"
    },
    sprintf(
      "  expected avoided %s: %.1f",
      .vpor(p$displaced$unit, "deaths"),
      p$displaced$deadisp
    ),
    "",
    "Burden of pollution",
    sprintf(
      "  attributable %s:   %.1f",
      .vpor(p$burden$extra$unit, "deaths"),
      p$burden$attrcases
    )
  )
  if (!is.null(p$equity)) {
    lines <- c(
      lines,
      "",
      "Equity analysis",
      sprintf("  concentration index: %.4f", p$equity$concidx)
    )
  }
  lines
}
