# First Nations, Metis and Inuit health and wellbeing: Statistics Canada's
# published tables, read as published (quality flags kept, suppressed cells
# never imputed), compared distinction by distinction, and the Truth and
# Reconciliation Commission's Calls to Action tracked. Data governance is
# part of the method: OCAP(R) (First Nations Information Governance Centre)
# and the CARE Principles (Carroll et al. 2020) say whose data these are.
# The computation is in C++ (src/igh.cpp and the files beside it).

#' Read a Statistics Canada table as published
#'
#' Tidies a Statistics Canada full-table download (the \dQuote{CSV for database
#' loading} form, ): \code{VALUE} as a number, and the
#' quality flags kept beside it in \code{quality}.
#'
#' @details Statistics Canada flags every estimate: \code{E} use with caution,
#' \code{F} too unreliable to be published, \code{x} suppressed to meet the
#' confidentiality requirements of the Statistics Act, \code{..} not available,
#' \code{...} not applicable. An \code{F}, \code{x}, \code{..} or \code{...}
#' cell has no value here (\code{NA}), whatever the file holds, and nothing in
#' the package imputes it or works it back from totals and margins: a suppressed
#' cell is suppressed to protect people, and small First Nations, Metis and
#' Inuit communities most of all. \code{E} values are kept and marked.
#'
#' The census data tables (98-10-xxxx) come in a wide layout: one column per
#' member of a dimension (\code{"Indigenous identity (9):Inuk (Inuit)[6]"}),
#' each followed by its \code{Symbol} column, and no \code{VALUE}. They are made
#' long first: the dimension becomes a column of member labels
#' (\code{"Indigenous identity (9)"}), the values \code{VALUE} and the symbols
#' \code{STATUS}, so every Statistics Canada table reads the same way.
#'
#' A label that names two members of one dimension (the census tables repeat
#' labels under different parents, and the file drops the hierarchy) gets its
#' member number from the coordinate: \code{"Without children [3]"} and
#' \code{"Without children [7]"}. Other labels are left as they are.
#'
#' @param x A data frame read from a Statistics Canada table (English-language
#'   full-table CSV).
#' @return The data frame with \code{VALUE} numeric and a factor column
#'   \code{quality} with levels \code{"ok"}, \code{"E"}, \code{"F"}, \code{"x"},
#'   \code{".."}, \code{"..."}, of class \code{c("sctidy", "data.frame")}.
#' @references Statistics Canada. Symbol legend.
#'   statcan.gc.ca/en/concepts/definitions/guide-symbol
#' @seealso \code{\link{indgap}}, \code{\link{scests}}.
#' @examples
#' # the layout of a downloaded table, with one cell suppressed (made-up
#' # numbers)
#' d <- data.frame(REF_DATE = "2019/2022", GEO = "Canada",
#'                 `Indigenous identity` = c("First Nations", "Inuk (Inuit)",
#'                                           "Non-Indigenous"),
#'                 VALUE = c("45.8", "", "61.5"), STATUS = c("", "x", ""),
#'                 check.names = FALSE)
#' sctidy(d)
#' @export
sctidy <- function(x) .igc("sctidy", x = x)

#' @export
format.sctidy <- function(x, ...) {
  fb <- .igplain(x, c("VALUE", "quality"), "Statistics Canada table, tidied")
  if (!is.null(fb)) {
    return(fb)
  }
  d <- x
  show <- d[setdiff(names(d), "STATUS")]
  more <- if (nrow(d) > 10L) sprintf("(%d more rows)", nrow(d) - 10L)
  q <- table(as.character(d$quality))
  .igb(
    "Statistics Canada table, tidied",
    facts = list(
      Rows = nrow(d),
      Columns = ncol(d),
      Quality = toString(sprintf("%s %d", names(q), as.integer(q)))
    ),
    tables = list(show[seq_len(min(10L, nrow(d))), , drop = FALSE]),
    reading = c(
      "VALUE is numeric; a suppressed cell (F, x, .., ...) is NA, with its",
      "flag kept in quality."
    ),
    notes = c(
      more,
      "Symbols: E use with caution; F too unreliable to be published;",
      "x suppressed (confidentiality); .. not available; ... not applicable."
    )
  )
}

#' Recode another agency's quality flags to the package's
#'
#' The comparisons of this family read Statistics Canada's quality flags:
#' \code{E} (use with caution), and \code{F}, \code{x}, \code{..}, \code{...}
#' (no value: too unreliable, suppressed, not available, not applicable). Other
#' agencies flag otherwise: the Australian Bureau of Statistics \code{np} (not
#' published) and \code{*} or \code{**} (relative standard error over 25% or
#' 50%), Stats NZ \code{S} (suppressed), the US Census Bureau \code{(X)}. This
#' maps them to the package's vocabulary, so suppressed cells are never compared
#' and uncertain ones are marked.
#'
#' @param x The flags, one per row.
#' @param suppr Flags meaning the value must not be used (it becomes
#'   \code{"x"}).
#' @param caution Flags meaning the value may be used with caution (it becomes
#'   \code{"E"}).
#' @return A character vector of \code{"ok"}, \code{"E"} and \code{"x"}; empty
#'   and missing flags are \code{"ok"}. \code{x} must be character (or a
#'   factor); a flag in neither \code{suppr} nor \code{caution} (other than
#'   \code{"ok"}) becomes \code{"ok"} and is named in a warning.
#' @seealso \code{\link{limfrom}}, \code{\link{indgap}}.
#' @examples
#' # flags as an Australian table might give them
#' qltyrcd(c("", "np", "*", "**", NA), suppr = c("np", "**"), caution = "*")
#' @export
qltyrcd <- function(x, suppr = c("F", "x", "..", "..."), caution = "E") {
  .igc("qltyrcd", x = x, suppr = suppr, caution = caution)
}

#' Confidence limits from a standard error, a coefficient of variation or a
#' margin of error
#'
#' The comparisons of this family take an estimate and its published confidence
#' limits. Agencies that publish a standard error, a relative standard error
#' (coefficient of variation) or a margin of error instead give the same
#' information; this turns it into limits at the level wanted.
#'
#' @details With an estimate \eqn{e} and \eqn{z = \Phi^{-1}(1 - \alpha/2)} for
#' the level wanted, a standard error \eqn{s} gives \eqn{e \pm z s}; a
#' coefficient of variation \eqn{c} (in percent) gives \eqn{s = e c / 100}; a
#' margin of error \eqn{m} at level \code{moelev} (the American Community Survey
#' publishes 90% margins) gives \eqn{s = m / z_{moe}}. The limits are symmetric
#' on the estimate's scale.
#'
#' @param est The estimates the limits are centred on.
#' @param se,cv,moe Give one: standard errors, coefficients of variation
#'   (relative standard errors, in percent), or margins of error.
#' @param moelev The confidence level of the margins of error (default 0.90).
#' @param conflvl The level of the limits returned (default 0.95).
#' @return A data frame: \code{estimate}, \code{lower}, \code{upper}, \code{se}.
#' @seealso \code{\link{qltyrcd}}, \code{\link{indgap}}.
#' @examples
#' limfrom(c(24.1, 31.5), cv = c(8, 12))
#' limfrom(1520, moe = 210)                 # a 90% margin of error
#' @export
limfrom <- function(
  est,
  se = NULL,
  cv = NULL,
  moe = NULL,
  moelev = 0.90,
  conflvl = 0.95
) {
  .igc(
    "limfrom",
    est = est,
    se = se,
    cv = cv,
    moe = moe,
    moelev = moelev,
    conflvl = conflvl
  )
}

#' @export
format.limfrom <- function(x, digits = 4L, ...) {
  title <- "Intervals from published precision"
  fb <- .igplain(x, c("estimate", "lower", "upper", "se"), title)
  if (!is.null(fb)) {
    return(fb)
  }
  tb <- data.frame(
    estimate = x$estimate,
    se = x$se,
    interval = .igci(x$lower, x$upper, digits)
  )
  .igb(
    title,
    facts = list(
      Estimates = nrow(x),
      Interval = "normal, estimate +/- z x se (default 95%)"
    ),
    tables = list(tb),
    reading = paste(
      "Each interval is rebuilt from the published standard error, CV or",
      "margin of error, assuming a normal sampling distribution."
    ),
    digits = digits
  )
}

#' Compare published group estimates with a reference group
#'
#' For published estimates with their confidence limits (percentages, rates,
#' years of life expectancy: any scale, with the limits on the same scale), the
#' ratio and the difference of each group to a reference group, with confidence
#' intervals worked from the published ones: First Nations, Metis and Inuit
#' against non-Indigenous people, Aboriginal and Torres Strait Islander people
#' against non-Indigenous Australians, Maori against non-Maori, or any groups
#' with published intervals.
#'
#' @details Each published interval gives a standard error, \eqn{s = (u - l) /
#' (2 z)}, on the scale of the estimate for the difference and on the log scale,
#' \eqn{\sigma = (\log u - \log l) / (2 z)}, for the ratio, with \eqn{z =
#' \Phi^{-1}(1 - \alpha/2)}. This takes the interval as symmetric on that scale;
#' for asymmetric (gamma, logit, Wilson) limits it averages the two sides. With
#' a standard error, a coefficient of variation or a margin of error instead of
#' limits, make limits with \code{\link{limfrom}} first. The groups are treated
#' as independent samples, so the standard error of a difference is
#' \eqn{\sqrt{s_g^2 + s_r^2}} and of a log ratio \eqn{\sqrt{\sigma_g^2 +
#' \sigma_r^2}}. Disjoint domains of one survey are nearly independent but not
#' exactly (they share the design's strata and clusters): the intervals are
#' approximations, and the exact ones need the survey's replicate weights (for
#' Statistics Canada's surveys, in its Research Data Centres). The reference
#' must not contain the groups: compare with non-Indigenous people, not with a
#' total, whose covariance with each group is large and positive (the function
#' warns when the reference looks like a total). A lower limit of 0 (a limit cut
#' at zero) leaves the ratio's interval \code{NA}.
#'
#' Distinctions-based: each group is compared on its own. A row that pools
#' groups is kept but marked \code{poogrp = TRUE}, since a pooled figure can
#' hide what differs between the groups it pools: name such rows in
#' \code{pooled} (by default labels beginning \code{"Total"}, \dQuote{Indigenous
#' identity}, \dQuote{Aboriginal identity} or \dQuote{Indigenous population} are
#' taken as pooled). A row flagged \code{E} (in the group or the reference) is
#' compared and marked; one flagged \code{F}, \code{x}, \code{..} or \code{...},
#' with no estimate, or whose published interval does not contain its own
#' estimate, is not compared. Other agencies' flags: \code{\link{qltyrcd}}.
#'
#' @param data A data frame, one row per group (filtered to one period, age
#'   group, sex and place: a group appearing twice is refused).
#' @param grp,est,lower,upper Column names: the group label, the estimate, and
#'   its published confidence limits, on one scale.
#' @param ref The label of the reference group (default
#'   \code{"Non-Indigenous"}).
#' @param qual Optional column of quality flags in the package's vocabulary
#'   (\code{"E"}, \code{"F"}, \code{"x"}, \code{".."}, \code{"..."}; see
#'   \code{\link{qltyrcd}} for other agencies' flags).
#' @param conflvl The confidence level of the published intervals and of the
#'   ones returned (default 0.95).
#' @param pooled Optional labels of the rows that pool other groups (for example
#'   \code{"Aboriginal and Torres Strait Islander"} beside its parts, or a
#'   \dQuote{Total Indigenous} row); by default the labels described in Details.
#'   Each label must be a group of \code{data}: one that matches no row (a
#'   typo) is refused rather than ignored.
#' @return A data frame of class \code{indgap}, one row per group
#'   but the reference. The columns are \code{group}, \code{estimate},
#'   \code{refest}, \code{ratio}, \code{ratiolo}, \code{ratioup},
#'   \code{difference}, \code{diflo}, \code{difup}, \code{usewithcau} and
#'   \code{poogrp}; attributes \code{reference} and \code{conflvl}.
#'
#' @references Statistics Canada (2025). Health indicators for First Nations
#' adults living off reserve, Metis and Inuit adults, 2015 to 2022. The Daily,
#' 2025-02-13.\cr
#' www150.statcan.gc.ca/n1/daily-quotidien/250213/dq250213f-eng.htm
#'
#' Carroll SR, Garba I, Figueroa-Rodriguez OL, et al. (2020). The CARE
#' Principles for Indigenous Data Governance. Data Science Journal 19:43.
#' \doi{10.5334/dsj-2020-043}
#' @seealso \code{\link{limfrom}}, \code{\link{qltyrcd}},
#'   \code{\link{dastd}}.
#' @examples
#' # made-up percentages for three groups and their pooled total, with
#' # published limits and one estimate flagged E
#' d <- data.frame(
#'   group = c("First Nations", "Metis", "Inuit", "Total Indigenous",
#'             "Non-Indigenous"),
#'   pct = c(45.8, 52.8, 47.8, 49, 61.5),
#'   lo = c(43.9, 50.6, 43.1, 47.5, 61),
#'   hi = c(47.7, 55, 52.5, 50.5, 62),
#'   q = c("", "", "E", "", "")
#' )
#' indgap(d, "group", "pct", "lo", "hi", qual = "q")
#' # your own table: made-up numbers in another country's layout, with relative
#' # standard errors instead of limits and that agency's flags
#' own <- data.frame(
#'   group = c(
#'     "Aboriginal",
#'     "Torres Strait Islander",
#'     "Aboriginal and Torres Strait Islander",
#'     "Non-Indigenous Australians"
#'   ),
#'   pct = c(31.2, 27.5, 30.8, 18.9),
#'   rse = c(4.1, 22.0, 3.9, 1.2),
#'   flag = c("", "*", "", "")
#' )
#' own <- cbind(own, limfrom(own$pct, cv = own$rse)[c("lower", "upper")])
#' own$q <- qltyrcd(own$flag, suppr = c("np", "**"), caution = "*")
#' indgap(
#'   own,
#'   "group",
#'   "pct",
#'   "lower",
#'   "upper",
#'   ref = "Non-Indigenous Australians",
#'   qual = "q",
#'   pooled = "Aboriginal and Torres Strait Islander"
#' )
#' @export
indgap <- function(
  data,
  grp,
  est,
  lower,
  upper,
  ref = "Non-Indigenous",
  qual = NULL,
  conflvl = 0.95,
  pooled = NULL
) {
  .igc(
    "indgap",
    data = data,
    grp = grp,
    est = est,
    lower = lower,
    upper = upper,
    ref = ref,
    qual = qual,
    conflvl = conflvl,
    pooled = pooled
  )
}

#' @export
format.indgap <- function(x, ...) {
  need <- c("group", "estimate", "refest", "ratio", "ratiolo", "difference")
  fb <- .igplain(x, need, "Gap against a reference group")
  if (!is.null(fb)) {
    return(fb)
  }
  lab <- .igpool(.igstar(x$group, x$usewithcau), x$poogrp)
  show <- data.frame(
    group = lab,
    estimate = x$estimate,
    refest = x$refest,
    ratio = x$ratio,
    rci = .igci(x$ratiolo, x$ratioup),
    difference = x$difference,
    dci = .igci(x$diflo, x$difup)
  )
  names(show)[c(5L, 7L)] <- c("ratio CI", "difference CI")
  .igb(
    "Gap against a reference group",
    facts = list(
      Reference = attr(x, "reference"),
      Groups = nrow(x),
      Intervals = .iglvl(attr(x, "conflvl"))
    ),
    tables = list(show),
    reading = .ignull(
      lab,
      x$ratiolo,
      x$ratioup,
      1,
      "Above the reference",
      "Below the reference"
    ),
    notes = c(
      "ratio: the group's estimate over the reference's; difference: the",
      "group's less the reference's, in the estimate's own units.",
      if (any(.igyes(x$poogrp))) {
        paste(
          "(pooled): a pan-Indigenous total, shown beside, not a",
          "distinct group (poogrp)."
        )
      },
      "Report a gap beside each group's own level, not alone.",
      .igcau(x$usewithcau)
    )
  )
}

#' Generations of residential school attendance in a family
#'
#' Codes the intergenerational exposure used in the epidemiology of Indian
#' Residential Schools: how many generations before the respondent attended, 0
#' (neither parents nor grandparents), 1 (parents or grandparents) or 2 (both a
#' parent and a grandparent), as in McQuaid et al. (2017).
#'
#' @details Each answer is \code{TRUE} (attended), \code{FALSE} (did not) or
#' \code{NA} (not known). An unknown answer may hide an attendance, so the count
#' is unknown (\code{NA}) whenever an unknown answer could change it: always
#' when the known answer is \code{FALSE}, and, with \code{partial = "na"} (the
#' default), also when the known answer is \code{TRUE} (the count is 1 or 2).
#' \code{partial = "lowerbnd"} codes that case 1, a lower bound: it keeps the
#' person but can classify them a level too low, which biases a trend toward no
#' effect.
#'
#' @param parent,gparent Logical vectors: did a parent, did a grandparent
#'   attend.
#' @param partial What a known attendance with the other generation unknown
#'   counts as: \code{"na"} (not known, the default) or \code{"lowerbnd"} (1).
#' @return An integer vector of 0, 1, 2 or \code{NA}.
#' @references McQuaid RJ, Bombay A, McInnis OA, Humeny C, Matheson K, Anisman H
#' (2017). Suicide ideation and attempts among First Nations peoples living
#' on-reserve in Canada: the intergenerational and cumulative effects of Indian
#' Residential Schools. Canadian Journal of Psychiatry 62(6):422-430.
#' \doi{10.1177/0706743717702075}
#'
#' Bombay A, Matheson K, Anisman H (2014). The intergenerational effects of
#' Indian Residential Schools: implications for the concept of historical
#' trauma. Transcultural Psychiatry 51(3):320-338.
#' \doi{10.1177/1363461513503380}
#' @seealso \code{\link{exptrnd}}.
#' @examples
#' genexpo(c(TRUE, FALSE, TRUE, NA, FALSE), c(TRUE, FALSE, FALSE, TRUE, NA))
#' genexpo(c(TRUE, FALSE, TRUE, NA, FALSE), c(TRUE, FALSE, FALSE, TRUE, NA),
#'                     partial = "lowerbnd")
#' @export
genexpo <- function(parent, gparent, partial = c("na", "lowerbnd")) {
  .igc("genexpo", parent = parent, gparent = gparent, partial = partial)
}

#' A dose-response across generations of exposure
#'
#' For counts of an outcome (suicidal thoughts, distress) among people with 0, 1
#' or 2 generations of residential school attendance in their family (or any
#' ordered exposure), the odds ratio of each level against the lowest level,
#' with Woolf intervals, and the Cochran-Armitage test for a trend across the
#' levels.
#'
#' @param cases Number with the outcome at each exposure level.
#' @param n Number of people at each level.
#' @param dose The exposure levels, in order (default \code{0, 1, 2, ...}); the
#'   trend test scores them.
#' @param conflvl Confidence level of the odds ratios' intervals (default 0.95).
#' @return A data frame with \code{dose}, \code{cases}, \code{n},
#'   \code{proportion}, \code{orratio}, \code{lower}, \code{upper}; the trend
#'   test is the attribute \code{trend} (an \code{htest}).
#' @details With \eqn{a_k} cases and \eqn{b_k = n_k - a_k} non-cases at level
#' \eqn{k}, the odds ratio against the lowest level is
#' \eqn{(a_k/b_k)/(a_1/b_1)}, with Woolf's interval \eqn{\exp\{\log OR \pm
#' z\sqrt{1/a_k + 1/b_k + 1/a_1 + 1/b_1}\}}. When any cell is empty, one half is
#' added to every cell of every level (Haldane-Anscombe), so every odds ratio is
#' finite. The trend test is the Cochran-Armitage test (\code{prop.trend.test})
#' with the levels as scores. Counts must be unweighted. For a complex survey,
#' take the weighted analysis from the survey itself (\code{\link{svymx}}):
#' this function is for published or aggregated counts.
#' @references McQuaid RJ, et al. (2017), see \code{\link{genexpo}}.
#' @examples
#' # made-up counts by 0, 1 and 2 generations of attendance, to show the
#' # computation
#' et <- exptrnd(cases = c(60, 95, 70), n = c(500, 520, 300))
#' et
#' attr(et, "trend")
#' @export
exptrnd <- function(cases, n, dose = seq_along(cases) - 1L, conflvl = 0.95) {
  .igc("exptrnd", cases = cases, n = n, dose = dose, conflvl = conflvl)
}

#' @export
format.exptrnd <- function(x, ...) {
  need <- c("dose", "cases", "n", "proportion", "orratio", "lower", "upper")
  fb <- .igplain(x, need, "Odds by level of exposure, with the trend")
  if (!is.null(fb)) {
    return(fb)
  }
  tr <- attr(x, "trend")
  show <- data.frame(
    dose = x$dose,
    cases = x$cases,
    n = x$n,
    proportion = x$proportion,
    orratio = x$orratio,
    interval = .igci(x$lower, x$upper)
  )
  k <- nrow(x)
  .igb(
    "Odds by level of exposure, with the trend",
    facts = list(
      Levels = k,
      People = sum(x$n),
      Cases = sum(x$cases),
      Trend = if (inherits(tr, "htest")) {
        sprintf(
          "chi-square %.2f on %s df, p %s",
          unname(tr$statistic),
          format(unname(tr$parameter)),
          .igp(tr$p.value)
        )
      }
    ),
    tables = list(show),
    reading = .etread(x, tr),
    notes = c(
      "orratio: the odds of a case at each level over the odds at the lowest",
      "(Woolf's interval); the trend test asks whether they rise with level.",
      if (any(x$cases == 0 | x$cases == x$n)) {
        "An empty cell: one half was added to every cell (Haldane-Anscombe)."
      }
    )
  )
}

#' Calls to Action completed, year by year
#'
#' The number of the Truth and Reconciliation Commission's 94 Calls to Action
#' that the Yellowhead Institute's annual accountability reports judged
#' complete, by the year they were completed, to the end of 2022: the only
#' counts here are those the reports state.
#'
#' @return A data frame with \code{period}, \code{yearend}, \code{completed},
#'   \code{cumulative} and \code{source}.
#' @details 2015-2019: 8 (13 in all by December 2022, less those of 2020-2022);
#' 2020: none; 2021: 3, all in June; 2022: 2 (Jewell and Mosby 2021, 2022).
#' Other trackers count otherwise: the federal government reports by theme
#' rather than call by call; CBC's Beyond 94 and Indigenous Watchdog grade each
#' call on their own scales and are copyrighted, so their tables are not
#' reproduced. The reports are licensed CC BY-NC 2.5 Canada; the counts are
#' facts taken from them, with the reports cited.
#' @references Jewell E, Mosby I (2021). Calls to Action Accountability: A 2021
#' Status Update on Reconciliation. Yellowhead Institute, Toronto Metropolitan
#' University.
#'
#' Jewell E, Mosby I (2022). Calls to Action Accountability: A 2022 Status
#' Update on Reconciliation. Yellowhead Institute, Toronto Metropolitan
#' University.
#'
#' Truth and Reconciliation Commission of Canada (2015). Truth and
#' Reconciliation Commission of Canada: Calls to Action. Winnipeg.
#' @seealso \code{\link{trcprog}}.
#' @examples
#' trccmpl()
#' @export
trccmpl <- function() .igc("trccmpl")

#' How fast the Calls to Action are being completed
#'
#' From cumulative counts of completed Calls to Action over time, the completion
#' rate (calls a year), its exact Poisson interval, and the year all 94 would be
#' complete at that rate.
#'
#' @details Completions are taken as a Poisson process with a constant rate over
#' the span from the start (\code{start}, the Commission's final report in
#' December 2015) to the last count: the rate is \eqn{c / T} with the exact
#' (Garwood) interval for \eqn{c} completed calls, and the projection is the
#' year \eqn{94 - c} more would take at that rate (at the interval's ends, its
#' range). A constant rate is the simplest model, and the one behind the
#' Yellowhead Institute's projection (\dQuote{42 years, or until 2065}, 2022);
#' with the counts to the end of 2022 and December 2015 as the start it gives
#' 2067. The pace has not been constant, so the projection is a description of
#' the pace so far, not a forecast. Any count of commitments completed over time
#' can be modelled the same way: set \code{total} and \code{start}.
#'
#' @param done Cumulative number completed at each date (default:
#'   \code{\link{trccmpl}}).
#' @param yearend The year each count was made to (end of year).
#' @param total Number of Calls to Action (94).
#' @param start When counting began, as a year with a fraction (2015.96:
#'   December 15, 2015).
#' @param conflvl Confidence level of the completion rate (default 0.95).
#' @return A one-row data frame of class \code{trcprog} with the
#'   columns \code{completed}, \code{remaining}, \code{years}, \code{rate},
#'   \code{ratelo}, \code{rateup}, \code{cmplby}, \code{cmplbyearly} and
#'   \code{cmplbylate};
#' @seealso \code{\link{trccmpl}}.
#' @examples
#' p <- trcprog()
#' p
#'
#' # your own count: Calls completed by the end of each year, made up
#' trcprog(done = c(2, 6, 9, 13), yearend = 2016:2019)
#' @export
trcprog <- function(
  done = trccmpl()$cumulative,
  yearend = trccmpl()$yearend,
  total = 94,
  start = 2015.96,
  conflvl = 0.95
) {
  .igc(
    "trcprog",
    done = done,
    yearend = yearend,
    total = total,
    start = start,
    conflvl = conflvl
  )
}

#' @export
format.trcprog <- function(x, ...) {
  need <- c("completed", "remaining", "years", "rate", "cmplby")
  if (!all(need %in% names(x)) || nrow(x) != 1L) {
    return(c("Truth and Reconciliation Calls to Action: progress", "", .igt(x)))
  }
  tot <- attr(x, "total")
  if (is.null(tot)) {
    tot <- x$completed + x$remaining
  }
  yr <- function(v) if (is.finite(v)) format(round(v)) else "never"
  .igb(
    "Truth and Reconciliation Calls to Action: progress",
    facts = list(
      Completed = sprintf("%d of %d", as.integer(x$completed), as.integer(tot)),
      Remaining = x$remaining,
      `Years counted` = x$years,
      `Calls a year` = sprintf("%.2f %s", x$rate, .igci(x$ratelo, x$rateup)),
      `All complete by` = sprintf(
        "%s (%s to %s)",
        yr(x$cmplby),
        yr(x$cmplbyearly),
        yr(x$cmplbylate)
      )
    ),
    reading = sprintf(
      "At this pace the remaining %d Calls are complete around %s.",
      as.integer(x$remaining),
      yr(x$cmplby)
    ),
    notes = c(
      "Completions taken as a Poisson process at a constant rate; the range",
      "comes from the rate's interval. Trackers count completion differently."
    )
  )
}

# Runs one of these functions in the compiled core.
.igc <- function(fn, ...) .Call(C_clab_r_igh, fn, list(...))

# The printed block of a result: a title, facts, tables, a reading, notes.
.igb <- function(
  title,
  facts = list(),
  tables = list(),
  reading = NULL,
  notes = NULL,
  digits = 4L
) {
  facts <- facts[!vapply(facts, is.null, NA)]
  out <- c(title, strrep("-", nchar(title, type = "width")))
  if (length(facts)) {
    k <- paste0(names(facts), ":")
    v <- vapply(
      facts,
      function(f) {
        paste(format(f, digits = digits), collapse = " ")
      },
      ""
    )
    out <- c(out, paste(formatC(k, width = -max(nchar(k, type = "width"))), v))
  }
  for (i in seq_along(tables)) {
    nm <- names(tables)[i]
    out <- c(
      out,
      "",
      if (!is.null(nm) && nzchar(nm)) nm,
      .igt(tables[[i]], digits)
    )
  }
  if (length(reading)) {
    out <- c(out, "", reading)
  }
  if (length(notes)) {
    out <- c(out, "", notes)
  }
  out
}

# A data frame as lines of text.
.igt <- function(d, digits = 4L) {
  n <- nrow(d)
  d <- structure(
    as.list(d),
    names = names(d),
    row.names = seq_len(n),
    class = "data.frame"
  )
  if (!n) {
    return("(no rows)")
  }
  f <- format(d, digits = digits)
  m <- rbind(names(f), as.matrix(f))
  m[is.na(m)] <- "NA"
  for (j in seq_len(ncol(m))) {
    m[, j] <- formatC(m[, j], width = max(nchar(m[, j], type = "width")))
  }
  apply(m, 1L, paste, collapse = "  ")
}

# A table with only the columns that are there (a subset of a result).
.igplain <- function(x, need, title) {
  if (all(need %in% names(x)) && nrow(x)) {
    return(NULL)
  }
  c(title, "", .igt(x))
}

.igci <- function(lo, hi, d = 3L) {
  ok <- is.finite(lo) & is.finite(hi)
  f <- function(v) trimws(formatC(v, digits = d, format = "g"))
  ifelse(ok, paste0("[", f(lo), ", ", f(hi), "]"), "")
}

# "A", "A and B", "A, B and C"
.iglist <- function(lab) {
  n <- length(lab)
  if (n < 2L) {
    return(paste(lab))
  }
  paste(toString(lab[-n]), "and", lab[n])
}

# TRUE where a logical flag is TRUE, FALSE where it is FALSE or missing
.igyes <- function(v) !is.na(v) & v

.igstar <- function(lab, flag) paste0(lab, ifelse(.igyes(flag), " *", ""))

.igpool <- function(lab, flag) {
  paste0(lab, ifelse(.igyes(flag), " (pooled)", ""))
}

.igcau <- function(flag) {
  if (!any(.igyes(flag))) {
    return(NULL)
  }
  paste(
    "* use with caution: the source flags the estimate (E, high",
    "sampling variability)."
  )
}

.iglvl <- function(l) if (is.numeric(l)) paste0(format(100 * l), "%")

# Which rows have an interval above, below or across the null value.
.ignull <- function(lab, lo, hi, null, up, down) {
  known <- is.finite(lo) & is.finite(hi)
  above <- known & lo > null
  below <- known & hi < null
  mid <- known & !above & !below
  nv <- format(null)
  c(
    if (any(above)) {
      sprintf("%s (interval above %s): %s.", up, nv, .iglist(lab[above]))
    },
    if (any(below)) {
      sprintf("%s (interval below %s): %s.", down, nv, .iglist(lab[below]))
    },
    if (any(mid)) {
      sprintf(
        "Not distinguishable from %s (interval includes it): %s.",
        nv,
        .iglist(lab[mid])
      )
    },
    if (!any(known)) "No row has an interval to judge it by."
  )
}

# The highest and the lowest of the groups, and whether their intervals
# part.
.igrange <- function(lab, est, what = "estimate", lo = NULL, hi = NULL) {
  ok <- which(is.finite(est))
  if (!length(ok)) {
    return(sprintf("No group has a %s to compare.", what))
  }
  if (length(ok) == 1L) {
    return(sprintf("%s: %s %s.", lab[ok], what, format(signif(est[ok], 4))))
  }
  i <- ok[which.max(est[ok])]
  j <- ok[which.min(est[ok])]
  if (est[i] == est[j]) {
    return(sprintf(
      "Every group has the same %s, %s.",
      what,
      format(signif(est[i], 4))
    ))
  }
  top <- sprintf(
    "Highest %s: %s (%s); lowest: %s (%s).",
    what,
    lab[i],
    format(signif(est[i], 4)),
    lab[j],
    format(signif(est[j], 4))
  )
  if (is.null(lo) || !is.finite(lo[i]) || !is.finite(hi[j])) {
    return(top)
  }
  c(
    top,
    if (lo[i] > hi[j]) {
      "Their intervals do not overlap: the gap is more than sampling error."
    } else {
      paste(
        "Their intervals overlap, which alone does not mean they are the",
        "same: the difference has its own, narrower interval."
      )
    }
  )
}

.igp <- function(p) {
  if (!length(p) || is.na(p)) {
    return("not available")
  }
  if (p < 0.001) "< 0.001" else format(signif(p, 3))
}

# format.exptrnd(): the odds at the highest level and the trend's verdict
.etread <- function(x, tr) {
  k <- nrow(x)
  c(
    if (is.finite(x$orratio[k])) {
      sprintf(
        paste(
          "At the highest level (%s) the odds of a case are %s times those",
          "at the lowest (%s)."
        ),
        format(x$dose[k]),
        format(signif(x$orratio[k], 3)),
        format(x$dose[1L])
      )
    },
    if (inherits(tr, "htest")) {
      if (tr$p.value < 0.05) {
        sprintf(
          "The odds %s with exposure: the trend test rejects a flat trend.",
          if (isTRUE(x$orratio[k] < 1)) "fall" else "rise"
        )
      } else {
        "The trend test finds no steady change in the odds with exposure."
      }
    }
  )
}
