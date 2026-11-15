# First Nations, Inuit and Metis frameworks for wellness, put to work on
# published tables: the First Nations Mental Wellness Continuum's four
# directions, Indigenous Services Canada's Community Well-Being index,
# connection to land and community as association across places, healing as
# change over time, a community-weighted index of protective and risk
# factors, and an audit of how deficit-framed a set of results is. Each says
# what it computes and what it cannot. The computation is in src/igf.cpp.

#' The four directions of the First Nations Mental Wellness Continuum, from
#' published indicators
#'
#' Lays out a group's published indicators under the four outcomes of the First
#' Nations Mental Wellness Continuum Framework (hope, belonging, meaning,
#' purpose), each with its interval, and says which directions the tables have
#' no measure of.
#'
#' @details The Framework (Health Canada and Assembly of First Nations 2015)
#' describes mental wellness as hope for the future, a sense of belonging and
#' connectedness to family, community and culture, a sense of meaning and of how
#' one's life fits into creation, and purpose in daily life. It has no numeric
#' scoring, and no survey table measures the four directions as the Framework
#' defines them. This function takes the indicators a user names as proxies for
#' each (\code{hope = "Always or often has a hopeful outlook"}, ...), reports
#' each direction's level (the mean of its proxies, their standard errors
#' combined as independent; the proxies share respondents, so if they are
#' positively correlated that standard error is too small), marks every row as a
#' proxy, and reports a direction given no proxy as missing rather than filling
#' it from another. One indicator named for two directions (Statistics Canada
#' asks about \dQuote{meaning and purpose} together) is used for both and marked
#' \code{shared}.
#'
#' The result is read as a \code{\link{welprof}}. For the balance across the
#' directions, call \code{\link{holbal}} on the rows with a proxy
#' (\code{x[x$domain != "no published proxy", ]}, with \code{"domain"} as its
#' domain column): a missing direction otherwise leaves the balance \code{NA},
#' and a shared proxy counts once per direction it serves, so give such
#' directions half weight. Other frameworks of wellness (Te Whare Tapa Wha, the
#' Social and Emotional Wellbeing domains, a Nation's own) take the same route
#' through \code{\link{welprof}} and \code{\link{holbal}}, which accept any set
#' of domains.
#'
#' @param data A data frame, one row per group and indicator.
#' @param grp,ind,est,lower,upper Column names: the group, the indicator,
#'   the estimate (percent) and its published limits.
#' @param hope,belong,meaning,purp The indicators taken as proxies for each
#'   direction, or \code{NULL} where the tables have none.
#' @param qual Optional column of quality flags: \code{"E"} rows are marked use
#'   with caution, and suppressed ones (\code{"F"}, \code{"x"}, \code{".."})
#'   get no estimate.
#' @param conflvl Confidence level of the published intervals (default 0.95).
#' @return A data frame with the classes \code{mhcont} and
#'   \code{welprof}. Its columns are \code{group},
#'   \code{indicator} (the direction), \code{side} (\code{"wellness"}),
#'   \code{domain} (the proxies), \code{estimate}, \code{lower}, \code{upper},
#'   \code{usewithcau}, \code{shared}, \code{proxy} (always \code{TRUE}).
#' @references Health Canada, Assembly of First Nations (2015). First Nations
#' Mental Wellness Continuum Framework. Ottawa: Health Canada. Summary
#' report:\cr sac-isc.gc.ca/eng/1576093687903/1576093725971
#' @seealso \code{\link{welprof}}, \code{\link{holbal}}.
#' @examples
#' # your own table: made-up percentages for two communities, with limits;
#' # one item stands for both meaning and purpose, and none for belonging
#' own <- data.frame(
#'   community = rep(c("Lakeside", "Riverbend"), each = 2),
#'   item = rep(c("Hopeful", "Life has meaning"), 2),
#'   pct = c(64, 70, 58, 75),
#'   lo = c(59, 65, 52, 70),
#'   hi = c(69, 75, 64, 80)
#' )
#' mhcont(own, "community", "item", "pct", "lo", "hi", hope = "Hopeful",
#'        meaning = "Life has meaning", purp = "Life has meaning")
#' @export
mhcont <- function(
  data,
  grp,
  ind,
  est,
  lower,
  upper,
  hope = NULL,
  belong = NULL,
  meaning = NULL,
  purp = NULL,
  qual = NULL,
  conflvl = 0.95
) {
  .igc(
    "mhcont",
    data = data,
    grp = grp,
    ind = ind,
    est = est,
    lower = lower,
    upper = upper,
    hope = hope,
    belong = belong,
    meaning = meaning,
    purp = purp,
    qual = qual,
    conflvl = conflvl
  )
}

#' @export
format.mhcont <- function(x, ...) {
  need <- c("group", "indicator", "domain", "estimate", "lower", "upper")
  fb <- .igplain(
    x,
    need,
    "First Nations Mental Wellness Continuum: the four directions"
  )
  if (!is.null(fb)) {
    return(fb)
  }
  miss <- unique(x$indicator[is.na(x$estimate)])
  ok <- !is.na(x$estimate) & !duplicated(x$indicator)
  .igb(
    "First Nations Mental Wellness Continuum: the four directions",
    facts = list(
      Groups = paste(unique(x$group), collapse = "; "),
      Directions = "hope, belonging, meaning, purpose",
      Estimates = "percent, from the published tables"
    ),
    tables = .mhtabs(x),
    reading = if (length(miss)) {
      c(
        paste("No published proxy:", toString(miss)),
        "-- reported as missing, never filled from another direction."
      )
    } else {
      "Every direction has a published proxy for every group."
    },
    notes = c(
      if (any(ok)) paste0("Proxy for ", x$indicator[ok], ": ", x$domain[ok]),
      "Every level is a proxy: no survey table measures the directions as",
      "the Framework defines them (Health Canada and AFN 2015).",
      if (any(.igyes(x$shared))) {
        "yes (shared): one indicator serves two directions, counted for each."
      },
      .igcau(x$usewithcau)
    )
  )
}

#' The Community Well-Being index of a community
#'
#' Indigenous Services Canada's Community Well-Being (CWB) index: income,
#' education, housing and labour force activity, each scored from 0 to 100, and
#' their mean, from census profile figures.
#'
#' @details The components, as the index's methodology sets them out (Penney,
#' O'Sullivan and Senecal 2012; O'Sullivan and McHardy 2013; Indigenous Services
#' Canada's reports on trends, 1981 to 2016 and 2021): \describe{
#' \item{income}{\eqn{100 (\log y - \log a) / (\log b - \log a)} for income per
#' capita \eqn{y} (total income over total population) between the bounds
#' \eqn{a} and \eqn{b} (2,650 and 75,000 dollars in the 2016 and 2021 series;
#' 2,000 and 40,000 in the earlier ones), kept between 0 and 100 here.}
#' \item{education}{two thirds of the percentage aged 20 and over with at least
#' a high school certificate, plus one third of the percentage aged 25 and over
#' with a university degree at the bachelor's level or above.}
#' \item{housing}{half the percentage of people in dwellings with no more than
#' one person per room, plus half the percentage in dwellings that need no major
#' repairs.} \item{labour}{half the participation rate, plus half the employment
#' rate of those participating (ages 20 to 64).} } The index is the mean of the
#' four. Indigenous Services Canada computes it from unrounded census microdata,
#' scores only communities of 65 people or more and suppresses the components of
#' communities under 250 people or 40 households: scores from published, rounded
#' profile figures come close to the published index but need not equal it, and
#' a small community's should be treated as the Department would. The full
#' current methodology is available from the Department on request. The index
#' measures socio-economic conditions as the census counts them; it says nothing
#' of culture, language or relationship to land, and is not a measure of
#' wellness as First Nations and Inuit define it. The CWB is Indigenous Services
#' Canada's index, with bounds in the Canadian dollars of its series: for
#' another country, set \code{bounds} in its currency and read the result as a
#' CWB-style score, not comparable with Canada's.
#'
#' @param incpercap Income per capita, in the dollars of the series' bounds.
#' @param highsch20p,uni25p Percentages for education.
#' @param notcro,nomajrep Percentages for housing.
#' @param partic,empl Percentages for labour force activity.
#' @param bounds The income bounds (default the 2016 and 2021 series: 2650 and
#'   75000).
#' @return A data frame, one row per community: \code{income}, \code{education},
#'   \code{housing}, \code{labour}, \code{cwb}.
#' @references Penney C, O'Sullivan E, Senecal S (2012). The Community
#' Well-Being Index (CWB): Examining Well-Being in Inuit Communities, 1981-2006.
#' Aboriginal Affairs and Northern Development Canada.
#'
#' O'Sullivan E, McHardy M (2013). The Community Well-Being (CWB) Index:
#' Disparity in Well-Being Between First Nations and Other Canadian Communities
#' Over Time. In White JP, Peters J, Beavon D, Spence N (eds), Aboriginal
#' Well-Being: Canada's Continuing Challenge. Thompson Educational.
#'
#' Indigenous Services Canada (2019). Report on trends in First Nations
#' communities, 1981 to 2016 (Community Well-Being index).\cr
#' sac-isc.gc.ca/eng/1345816651029/1557323327644
#' @examples
#' # made-up profile figures for two communities, to show the computation
#' cwbscr(incpercap = c(18000, 31000), highsch20p = c(58, 81),
#'           uni25p = c(6, 19), notcro = c(71, 95),
#'           nomajrep = c(62, 90),
#'           partic = c(55, 68), empl = c(78, 92))
#' @export
cwbscr <- function(
  incpercap,
  highsch20p,
  uni25p,
  notcro,
  nomajrep,
  partic,
  empl,
  bounds = c(2650, 75000)
) {
  .igc(
    "cwbscr",
    incpercap = incpercap,
    highsch20p = highsch20p,
    uni25p = uni25p,
    notcro = notcro,
    nomajrep = nomajrep,
    partic = partic,
    empl = empl,
    bounds = bounds
  )
}

#' @export
format.cwbscr <- function(x, digits = 4L, ...) {
  title <- "Community Well-Being scores"
  fb <- .igplain(x, c("income", "education", "housing", "labour", "cwb"), title)
  if (!is.null(fb)) {
    return(fb)
  }
  ok <- !is.na(x$cwb)
  .igb(
    title,
    facts = list(
      Communities = nrow(x),
      `Mean CWB` = if (any(ok)) mean(x$cwb[ok]),
      `CWB range` = if (any(ok)) .igci(min(x$cwb[ok]), max(x$cwb[ok]))
    ),
    tables = list(data.frame(community = seq_len(nrow(x)), as.list(x))),
    reading = paste(
      "Each component runs from 0 to 100, and the CWB score is their mean:",
      "income, education, housing and labour force activity."
    ),
    digits = digits
  )
}

#' Connection across places: does one indicator rise where another does?
#'
#' For two indicators published for the same places (regions, provinces, Inuit
#' regions), their rank correlation across places, an exact permutation p-value,
#' and an interval that carries the sampling error of the published estimates:
#' for connection to land and wellness, community belonging and mental health,
#' education and health.
#'
#' @details Mashford-Pringle and Stewart (2019) describe the land as \dQuote{an
#' instrumental part of all aspects of Indigenous life including health and
#' well-being}; Mashford-Pringle (2016) describes Indigenous health and
#' education as interconnected. Across places, that connection has a testable
#' shadow: where more people are on the land (or belong, or finish school), do
#' more report wellness? This function measures it by Spearman's \eqn{\rho}. Its
#' p-value is exact when there are at most 8 places (over all \eqn{n!}
#' orderings) and a Monte Carlo estimate otherwise (\eqn{(b + 1)/(R + 1)} for
#' \eqn{b} of \eqn{R} = \code{reps} random orderings at least as extreme). Each
#' estimate carries sampling error, so \code{rholo} and \code{rhoup} are the
#' quantiles of \eqn{\rho} over \code{reps} draws of both indicators from normal
#' distributions with their published standard errors. That is a sensitivity
#' band, not a confidence interval for a correlation across places: it shows how
#' much the published sampling error alone moves \eqn{\rho} for these places, it
#' is pulled toward 0 (noise is added to estimates that are already noisy), and
#' it ignores the sampling of the places themselves. Quality flags are not
#' carried: filter suppressed rows out first.
#'
#' This is an association across places (an ecological one). It cannot say that
#' connection to land makes a person well, which needs data on persons; a
#' handful of places gives little power; and a correlation is not a measure of
#' the relationship the authors describe. It is a description of how the
#' published indicators move together.
#'
#' @param data A data frame, one row per place.
#' @param place Column of place names.
#' @param x,y Columns of the two indicators.
#' @param xlo,xup,ylo,yup Optional columns of their published limits (for
#'   the interval).
#' @param reps Draws for the interval and, past 8 places, for the p-value
#'   (default 10,000).
#' @param seed Seed for the draws (default 1).
#' @param conflvl Confidence level of rho's drawn interval (default 0.95).
#' @return A list of class \code{connassoc}: \code{data} (the places
#'   used), \code{rho}, \code{p_value}, \code{exact}, \code{rholo},
#'   \code{rhoup}, \code{n}.
#' @section Acknowledgement: This estimator puts to quantitative use the work of
#' Professor Angela Mashford-Pringle (Timiskaming First Nation) and Professor
#' Suzanne Stewart (Yellowknives Dene First Nation), both of the University of
#' Toronto.
#'
#' See the acknowledgements in the section on Indigenous health of the
#' package's vignette (\code{vignette("crimlab")}).
#' @references Mashford-Pringle A, Stewart SL (2019). Akiikaa (it is the land):
#' exploring land-based experiences with university students in Ontario. Global
#' Health Promotion 26(3_suppl):64-72. \doi{10.1177/1757975919828722}
#'
#' Mashford-Pringle A (2016). Aboriginal health and Aboriginal education are
#' interconnected. Research in the Sociology of Health Care 34:85-96.
#' \doi{10.1108/s0275-495920160000034005}
#' @seealso \code{\link{dualcont}}.
#' @examples
#' # your own places: made-up shares in six regions, without limits
#' own <- data.frame(
#'   region = c("North", "South", "East", "West", "Coast", "Interior"),
#'   language = c(41, 18, 27, 33, 22, 38),
#'   wellbeing = c(68, 55, 60, 66, 57, 63)
#' )
#' connassoc(own, "region", "language", "wellbeing", reps = 2000)
#' @export
connassoc <- function(
  data,
  place,
  x,
  y,
  xlo = NULL,
  xup = NULL,
  ylo = NULL,
  yup = NULL,
  reps = 10000,
  seed = 1,
  conflvl = 0.95
) {
  .igc(
    "connassoc",
    data = data,
    place = place,
    x = x,
    y = y,
    xlo = xlo,
    xup = xup,
    ylo = ylo,
    yup = yup,
    reps = reps,
    seed = seed,
    conflvl = conflvl
  )
}

#' @export
format.connassoc <- function(x, ...) {
  lab <- x$labels
  rho <- x$rho
  size <- if (!is.finite(rho)) {
    NULL
  } else if (abs(rho) >= 0.7) {
    "strongly"
  } else if (abs(rho) >= 0.4) {
    "moderately"
  } else {
    "weakly"
  }
  .igb(
    "Connection across places",
    facts = list(
      Places = x$n,
      Indicators = sprintf("%s and %s", lab[[1L]], lab[[2L]]),
      `Rank correlation` = sprintf("rho = %.2f", rho),
      `p-value` = sprintf(
        "%s (%s)",
        .igp(x$p_value),
        if (x$exact) {
          "exact, every permutation"
        } else {
          "permutation"
        }
      ),
      Interval = if (is.finite(x$rholo)) {
        sprintf(
          "%s %s, with the published sampling error",
          .iglvl(x$conflvl),
          .igci(x$rholo, x$rhoup, 2L)
        )
      }
    ),
    reading = if (is.null(size)) {
      "The rank correlation is undefined (an indicator does not vary)."
    } else {
      sprintf(
        "Across these places the two rank %s %s.",
        size,
        if (rho >= 0) "together" else "in opposite orders"
      )
    },
    notes = c(
      "An association across places only:",
      "it does not show that one causes the other for any person."
    )
  )
}

#' Healing as a process: change in published indicators between two periods
#'
#' For each group and indicator, the change between two survey periods with its
#' confidence interval, wellness indicators first: change over time read as
#' Indigenous scholars describe mental health, as a process of healing, not as a
#' fixed deficit.
#'
#' @details Stewart (2008) describes Indigenous mental health as wellness and as
#' a process of healing; a process is seen in change. The change \eqn{p_2 - p_1}
#' has the interval \eqn{p_2 - p_1 \pm z\sqrt{s_1^2 + s_2^2}}, with the standard
#' errors from the published intervals and the two periods' samples independent
#' (separate cycles of a survey are). A change is only as comparable as the two
#' periods' questions and designs: Statistics Canada flags breaks in a series,
#' and a change across one is not a change in the people.
#'
#' @param data A data frame, one row per group, indicator and period.
#' @param grp,ind,prd,est,lower,upper Column names: the change is in the
#'   estimate's own units (percentage points for percentages).
#' @param from,to The two periods.
#' @param well,distress Optional indicators of wellness (listed first) and of
#'   distress; with \code{distress} given, any other indicator is marked
#'   \code{"other"}.
#' @inheritParams mhcont
#' @param conflvl Confidence level of each change's interval (default 0.95).
#' @return A data frame of class \code{healchg}: \code{group},
#'   \code{indicator}, \code{side}, \code{from}, \code{to}, \code{change},
#'   \code{lower}, \code{upper}, \code{usewithcau}.
#' @section Acknowledgement: This estimator puts to quantitative use the
#' framework of Professor Suzanne Stewart (Yellowknives Dene First Nation,
#' University of Toronto).
#'
#' See the acknowledgements in the section on Indigenous health of the
#' package's vignette (\code{vignette("crimlab")}).
#' @references Stewart SL (2008), see \code{\link{holbal}}.
#' @seealso \code{\link{welprof}}.
#' @examples
#' # your own table: made-up percentages for one community in two surveys
#' own <- data.frame(
#'   community = "Lakeside",
#'   item = rep(c("Good mental health", "Anxiety"), 2),
#'   survey = rep(c("2018", "2023"), each = 2),
#'   pct = c(55, 14, 61, 11),
#'   lo = c(50, 11, 56, 8),
#'   hi = c(60, 17, 66, 14)
#' )
#' healchg(own, "community", "item", "survey", "pct", "lo", "hi",
#'         from = "2018", to = "2023", well = "Good mental health",
#'         distress = "Anxiety")
#' @export
healchg <- function(
  data,
  grp,
  ind,
  prd,
  est,
  lower,
  upper,
  from,
  to,
  well = NULL,
  qual = NULL,
  conflvl = 0.95,
  distress = NULL
) {
  .igc(
    "healchg",
    data = data,
    grp = grp,
    ind = ind,
    prd = prd,
    est = est,
    lower = lower,
    upper = upper,
    from = from,
    to = to,
    well = well,
    qual = qual,
    conflvl = conflvl,
    distress = distress
  )
}

#' @export
format.healchg <- function(x, ...) {
  need <- c("group", "indicator", "side", "from", "to", "change", "lower")
  fb <- .igplain(x, need, "Healing as change between two periods")
  if (!is.null(fb)) {
    return(fb)
  }
  ind <- unique(x$indicator)
  pr <- attr(x, "periods")
  up <- x$lower > 0
  dn <- x$upper < 0
  cnt <- function(v, s) sum(v & x$side %in% s, na.rm = TRUE)
  .igb(
    "Healing as change between two periods",
    facts = list(
      Periods = if (!is.null(pr)) sprintf("%s to %s", pr[["from"]], pr[["to"]]),
      Groups = paste(unique(x$group), collapse = "; "),
      Indicators = length(ind),
      Units = "the estimate's own (percentage points for percentages)"
    ),
    tables = .hctabs(x, ind),
    reading = c(
      sprintf(
        "%d of %d changes have an interval clear of zero.",
        sum(up | dn, na.rm = TRUE),
        sum(!is.na(x$change))
      ),
      if (!all(is.na(x$side))) {
        sprintf(
          paste(
            "Wellness: %d rose and %d fell by more than chance; distress:",
            "%d rose and %d fell."
          ),
          cnt(up, "wellness"),
          cnt(dn, "wellness"),
          cnt(up, "distress"),
          cnt(dn, "distress")
        )
      }
    ),
    notes = c(
      "Wellness indicators come first: a change is read as part of healing.",
      "A change is only as comparable as the two periods' questions and",
      "designs: across a break in the series it is not a change in people.",
      .igcau(x$usewithcau)
    )
  )
}

#' A community-weighted index of protective and risk factors
#'
#' Combines published indicators into one index for each group or place, with
#' weights the community sets: protective factors count for, risk factors
#' against. The weights are never estimated.
#'
#' @details Inuit Tapiriit Kanatami's National Inuit Suicide Prevention Strategy
#' (2016) names risk factors (among them crowded housing, food insecurity,
#' childhood adversity, mental distress) and protective factors (among them
#' cultural continuity, social equity, family strength, healthy child
#' development). How much each matters is a judgement for the people the index
#' describes, so the function takes the weights as given: positive for
#' protective, negative for risk. The index is \eqn{\sum_i w_i p_i / \sum_i
#' |w_i|}, with standard error \eqn{\sqrt{\sum_i w_i^2 s_i^2} / \sum_i |w_i|}
#' from the published intervals (independent estimates). It is a summary of the
#' indicators under the weights given; with other weights it ranks places
#' otherwise, and it is never a ranking of Nations or communities on their
#' culture. For Inuit regions, use Inuit data only.
#'
#' @param data A data frame, one row per group (or place) and indicator.
#' @param grp,ind,est,lower,upper Column names; every indicator must be on one
#'   scale (percent, for example: a rate per 1,000 beside a percentage makes the
#'   weights meaningless).
#' @param w Named weights, one per indicator used: positive for a protective
#'   factor, negative for a risk factor.
#' @inheritParams mhcont
#' @param conflvl Confidence level of each group's index limits (default 0.95).
#' @return A data frame of class \code{commidx}: \code{group},
#'   \code{index}, \code{lower}, \code{upper}, \code{usewithcau}; attributes
#'   \code{weights} and \code{conflvl}.
#' @references Inuit Tapiriit Kanatami (2016). National Inuit Suicide Prevention
#' Strategy. Ottawa. itk.ca/national-inuit-suicide-prevention-strategy
#'
#' Carroll SR, et al. (2020), see \code{\link{holbal}}.
#' @examples
#' # your own table and the weights a community chose: made-up percentages
#' own <- data.frame(
#'   community = rep(c("Lakeside", "Riverbend"), each = 3),
#'   item = rep(c("Speaks the language", "Has a doctor", "Food insecure"), 2),
#'   pct = c(38, 72, 21, 25, 64, 30),
#'   lo = c(33, 67, 17, 20, 58, 25),
#'   hi = c(43, 77, 25, 30, 70, 35)
#' )
#' commidx(own, "community", "item", "pct", "lo", "hi",
#'         w = c("Speaks the language" = 2, "Has a doctor" = 1,
#'               "Food insecure" = -1))
#' @export
commidx <- function(
  data,
  grp,
  ind,
  est,
  lower,
  upper,
  w = NULL,
  qual = NULL,
  conflvl = 0.95
) {
  .igc(
    "commidx",
    data = data,
    grp = grp,
    ind = ind,
    est = est,
    lower = lower,
    upper = upper,
    w = w,
    qual = qual,
    conflvl = conflvl
  )
}

#' @export
format.commidx <- function(x, ...) {
  fb <- .igplain(
    x,
    c("group", "index", "lower", "upper"),
    "Community-weighted index of protective and risk factors"
  )
  if (!is.null(fb)) {
    return(fb)
  }
  w <- attr(x, "weights")
  lab <- .igstar(x$group, x$usewithcau)
  show <- data.frame(
    group = lab,
    index = x$index,
    interval = .igci(x$lower, x$upper)
  )
  .igb(
    "Community-weighted index of protective and risk factors",
    facts = list(
      Groups = nrow(x),
      Indicators = length(w),
      Index = "sum of w x estimate over sum of |w|, on the indicators' scale"
    ),
    tables = list(show),
    reading = .igrange(lab, x$index, "index", x$lower, x$upper),
    notes = c(
      "A higher index: more of the protective factors, less of the risks,",
      "under the weights given.",
      if (length(w)) "Weights (positive protective, negative risk):",
      if (length(w)) sprintf("%+g  %s", unname(w), names(w)),
      "The weights are a judgement for the people the index describes;",
      "with other weights it ranks places otherwise. It is never a ranking",
      "of Nations or communities on their culture.",
      "Risk and protective factors: Inuit Tapiriit Kanatami (2016).",
      .igcau(x$usewithcau)
    )
  )
}

#' How deficit-framed is a set of results?
#'
#' Counts how many of the indicators in a report measure distress, deficit or
#' disparity and how many measure wellness or strength, and names any group
#' reported on deficits alone.
#'
#' @details Walter (2016) shows that statistics about Indigenous peoples cluster
#' around five Ds: disparity, deprivation, disadvantage, dysfunction and
#' difference. Results built only of those describe Indigenous peoples as their
#' problems. This audit does not measure anything about people; it measures the
#' report: the share of its indicators that are deficit-framed, overall and per
#' group, so a report can be balanced before it is published. Pair each gap with
#' a strength-framed level (\code{\link{welprof}}, \code{\link{dualcont}}).
#'
#' @param x A \code{\link{welprof}} result, whose rows are counted by side.
#' @return A list of class \code{frmaud}: \code{bygrp} (a data frame of
#'   \code{group}, \code{wellness}, \code{distress}, \code{shrdef}),
#'   \code{shrdef} overall and \code{defonly} (the groups with no wellness
#'   indicator).
#' @references Walter M (2016). Data politics and Indigenous representation in
#' Australian statistics. In Kukutai T, Taylor J (eds), Indigenous Data
#' Sovereignty: Toward an Agenda, 79-97. ANU Press.
#' \doi{10.22459/caepr38.11.2016.05}
#'
#' Walter M, Andersen C (2013). Indigenous Statistics: A Quantitative Research
#' Methodology. Left Coast Press. \doi{10.4324/9781315426570}
#' @examples
#' # your own table: made-up percentages for two communities, with limits
#' own <- data.frame(
#'   community = rep(c("Lakeside", "Riverbend"), each = 4),
#'   item = rep(c("Hopeful", "Satisfied", "Rarely hopeful", "Low satisfaction"),
#'              2),
#'   pct = c(64, 71, 12, 18, 58, 66, 16, 22),
#'   lo = c(59, 66, 9, 14, 52, 60, 12, 17),
#'   hi = c(69, 76, 15, 22, 64, 72, 20, 27)
#' )
#' frmaud(welprof(own, "community", "item", "pct", "lo", "hi",
#'                well = c("Hopeful", "Satisfied"),
#'                distress = c("Rarely hopeful", "Low satisfaction")))
#' @export
frmaud <- function(x) .igc("frmaud", x = x)

#' @export
format.frmaud <- function(x, ...) {
  bg <- x$bygrp
  .igb(
    "Framing audit",
    facts = list(
      Groups = nrow(bg),
      `Wellness indicators` = sum(bg$wellness),
      `Distress indicators` = sum(bg$distress),
      `Deficit-framed` = sprintf("%.0f%%", 100 * x$shrdef)
    ),
    tables = list(`By group (shrdef: the deficit-framed share)` = bg),
    reading = c(
      sprintf(
        "%.0f%% of the indicators reported are deficit-framed.",
        100 * x$shrdef
      ),
      if (length(x$defonly)) {
        paste("Reported on deficits alone:", paste(x$defonly, collapse = "; "))
      },
      if (x$shrdef > 0.5) {
        c(
          "Most of this report measures disparity or deficit;",
          "pair each with a strength-framed level."
        )
      }
    ),
    notes = c(
      "Walter (2016): statistics about Indigenous peoples cluster around",
      "disparity, deprivation, disadvantage, dysfunction and difference.",
      "This audits the report, not the people it describes."
    )
  )
}

# format.mhcont(): one table per group, a direction to a row
.mhtabs <- function(x) {
  tabs <- lapply(unique(x$group), function(g) {
    r <- x[x$group == g, , drop = FALSE]
    proxy <- ifelse(.igyes(r$shared), "yes (shared)", "yes")
    proxy[is.na(r$estimate)] <- "none published"
    data.frame(
      direction = .igstar(r$indicator, r$usewithcau),
      estimate = r$estimate,
      interval = ifelse(is.na(r$estimate), "", .igci(r$lower, r$upper)),
      proxy = proxy
    )
  })
  names(tabs) <- unique(x$group)
  tabs
}

# format.healchg(): one table per indicator, named with its side
.hctabs <- function(x, ind) {
  tabs <- lapply(ind, function(k) {
    r <- x[x$indicator == k, , drop = FALSE]
    data.frame(
      group = .igstar(r$group, r$usewithcau),
      from = r$from,
      to = r$to,
      change = r$change,
      interval = .igci(r$lower, r$upper)
    )
  })
  side <- x$side[match(ind, x$indicator)]
  names(tabs) <- ifelse(is.na(side), ind, sprintf("%s (%s)", ind, side))
  tabs
}
