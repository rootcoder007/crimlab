# First Nations, Metis and Inuit mental health and wellbeing from Statistics
# Canada's published tables: the estimates and their limits taken out of a
# table's statistics column, rated scales (perceived mental health, life
# satisfaction) compared as whole distributions, survey percentages
# standardised by age, wellness beside distress and its balance across
# domains. The computation is in src/igm.cpp.

#' The estimates of a Statistics Canada table, with their confidence limits
#'
#' Statistics Canada's tables give each estimate on several rows: the estimate
#' itself (a percent, a rate, a count, a life expectancy), its low and high
#' confidence limits, or its coefficient of variation or standard error, all in
#' one statistics column (\code{Statistics}, \code{Characteristics},
#' \code{Statistics (3)}, ...). This puts them side by side, one row per cell of
#' the table and measure, ready for \code{\link{indgap}}, \code{\link{ordcmp}},
#' \code{\link{stdprev}} and the rest. Every table reads: one without a
#' statistics column (a census count) gives each value as its own estimate, with
#' no limits.
#'
#' @details Each label of the statistics column is read as an estimate or as one
#' of its measures of precision, in every form the tables use: \itemize{ \item
#' \code{"Low 95\% confidence interval"}; \item
#' \code{"High 95\% confidence interval, age-standardized rate"}; \item
#' \code{"Confidence interval lower bound of count"}; \item
#' \code{"95\% confidence interval (upper limit)"};
#' \item \code{"95\% confidence interval lower bound, 2020"}; \item
#' \code{"Coefficient of variation for percent"}; \item \code{"Standard error"}.
#' } A limit that names its measure is matched to it (ignoring case and
#' spacing); one that names none, or a name that fits several measures
#' (\code{"potential years of life lost"} in a table with a number, a crude rate
#' and an age-standardised rate of them), is matched to the measure with the
#' same unit (\code{UOM}) whose estimates it contains in every cell, and left
#' unmatched (said in the attribute \code{unmatched}) if no single measure fits.
#' A coefficient of variation \eqn{c} (in percent) gives log-normal limits
#' \eqn{e \exp(\mp z \sqrt{\log(1 + (c/100)^2)})}, which stay positive however
#' large \eqn{c} is; a standard error \eqn{s} the limits \eqn{e \mp z s}; the
#' column \code{cisou} says which (\code{"published"}, \code{"cv"} or
#' \code{"se"}). Published limits that do not contain their own estimate (it
#' happens: a 100% with limits 85.5 to 98.6 in table 41-10-0012) are kept as
#' published and marked \code{"inconsistent"}; every comparison of the family
#' warns and leaves such a row out. A limit is matched to a measure when it fits
#' 95% of the cells. Labels are read in English: use the English-language
#' tables.
#'
#' Every other column but the metadata (\code{VECTOR}, \code{COORDINATE},
#' \code{UOM}, ...) identifies the cell, \code{DGUID} among them (one place name
#' can name several geographies), with the spaces some labels end in trimmed.
#' The quality flag is the estimate's: an \code{F} or \code{x} estimate is
#' \code{NA} (see \code{\link{sctidy}}), and so are its limits.
#'
#' @param x A table read by \code{\link{sctidy}}, or a data frame as
#'   downloaded.
#' @return A data frame, one row per cell and measure: the identifying columns,
#'   \code{measure}, \code{estimate}, \code{lower}, \code{upper}, \code{cisou}
#'   and \code{quality}; the attribute \code{conflvl} is the level of the
#'   limits, and \code{unmatched} names any limits that could not be matched to
#'   one measure.
#' @references Statistics Canada (2023). Self-perceived mental health of First
#'   Nations people living off reserve, Metis and Inuit by age group and gender.
#'   Table 41-10-0067-01. \doi{10.25318/4110006701-eng}
#' @seealso \code{\link{sctidy}}, \code{\link{indgap}},
#'   \code{\link{ordcmp}}.
#' @examples
#' # the layout of a downloaded table: percent and number of persons, with
#' # limits for the percent
#' d <- data.frame(
#'   REF_DATE = "2022",
#'   GEO = "Canada",
#'   `Indigenous identity` = rep(c("First Nations", "Metis"), each = 4),
#'   Statistics = rep(
#'     c(
#'       "Number of persons",
#'       "Percent",
#'       "Low 95% confidence interval",
#'       "High 95% confidence interval"
#'     ),
#'     2
#'   ),
#'   UOM = rep(c("Persons", "Percent", "Percent", "Percent"), 2),
#'   VALUE = c(415000, 47.1, 45.2, 49.0, 320000, 50.3, 48.1, 52.5),
#'   STATUS = "",
#'   check.names = FALSE
#' )
#' scests(d)
#' # a rate with its coefficient of variation instead of limits
#' r <- data.frame(
#'   GEO = "Canada",
#'   Characteristics = c("Rate", "Coefficient of variation for rate"),
#'   VALUE = c(12.5, 8),
#'   STATUS = ""
#' )
#' scests(r)
#' @export
scests <- function(x) .igc("scests", x = x)

#' Compare rated scales (perceived mental health, life satisfaction) between
#' groups
#'
#' For a rated scale published as the percentage in each category with its
#' confidence limits (perceived mental health from poor to excellent, life
#' satisfaction from 0 to 10, agreement with a statement), the probability that
#' a person in each group rates higher than a person in the reference group,
#' with its confidence interval: a comparison of the whole distribution, not of
#' one category cut from it.
#'
#' @details With shares \eqn{p_i} in the group and \eqn{q_j} in the reference
#' over categories ordered from lowest to highest, the probability of a higher
#' rating, ties counted half, is \eqn{A = \sum_i p_i (Q_{i-1} + q_i / 2)} with
#' \eqn{Q_i} the reference's cumulative shares and \eqn{Q_0 = 0}: the
#' Mann-Whitney probability, equal to \eqn{\alpha / (1 + \alpha)} for Agresti's
#' (1980) generalized odds ratio \eqn{\alpha} with ties split equally, and to
#' the mean ridit of the group with the reference as the standard (Bross 1958).
#' 0.5 is no difference; \eqn{2A - 1} is Cliff's \eqn{\delta}, the share rating
#' higher less the share rating lower.
#'
#' Its variance is that of the two-sample U statistic, \eqn{\sigma_1^2 / n_1 +
#' \sigma_0^2 / n_0}, with \eqn{\sigma_1^2} the variance over the group of each
#' category's placement in the reference and \eqn{\sigma_0^2} the converse
#' (Brunner and Munzel 2000). The sample sizes are the effective ones the
#' published limits imply: for each category, \eqn{n = p(1 - p) / s^2} with
#' \eqn{s} the standard error from the limits' width, taking the median over the
#' categories. That carries the design effect of the survey's weights and
#' clusters, which the counts of persons (population totals) do not. The
#' interval is built on the logit scale, so it stays inside 0 and 1.
#'
#' The categories given in \code{lvls} are taken as the whole scale: rows for
#' totals, for combined categories (\dQuote{Excellent or very good}) and for
#' \dQuote{not specified} are left out, and the shares are rescaled to sum to
#' one over the levels. A group with a suppressed or missing category is not
#' compared, since its distribution is unknown; one with an \code{E} category
#' (in the group or the reference) is compared and marked, and one whose
#' published interval misses its own estimate is not compared.
#'
#' A higher rating means better only on a scale that runs that way. Give
#' \code{lvls} from the least to the most favourable so that \code{probbet}
#' reads as its name says.
#'
#' @param data A data frame, one row per group and category (as
#'   \code{\link{scests}} gives, filtered to one geography, age group and
#'   gender).
#' @param grp,cat,est,lower,upper Column names: the group, the category, the
#'   share in it (percent or proportion: only relative shares matter) and its
#'   published limits on the same scale.
#' @param lvls The categories, from the lowest (least favourable) to the
#'   highest.
#' @param ref The label of the reference group (required: any group of
#'   \code{data}).
#' @param qual Optional column of quality flags (see \code{\link{qltyrcd}}).
#' @param conflvl Confidence level of the published limits and of the result
#'   (default 0.95).
#' @return A data frame of class \code{ordcmp}: \code{group},
#'   \code{reference}, \code{probbet}, \code{lower}, \code{upper}, \code{delta},
#'   \code{neff}, \code{neffref}, \code{usewithcau}, \code{poogrp}; the
#'   attribute \code{shares} is the matrix of shares (groups by levels), which
#'
#' @references Agresti A (1980). Generalized odds ratios for ordinal data.
#' Biometrics 36(1):59-67. \doi{10.2307/2530495}
#'
#' Bross IDJ (1958). How to use ridit analysis. Biometrics 14(1):18-38.
#' \doi{10.2307/2527727}
#'
#' Brunner E, Munzel U (2000). The nonparametric Behrens-Fisher problem:
#' asymptotic theory and a small-sample approximation. Biometrical Journal
#' 42(1):17-25.
#' \doi{10.1002/(SICI)1521-4036(200001)42:1<17::AID-BIMJ17>3.0.CO;2-U}
#' @param pooled Optional labels of rows that pool other groups (see
#'   \code{\link{indgap}}).
#' @seealso \code{\link{scests}}, \code{\link{indgap}}.
#' @examples
#' # your own table: a made-up four-level rating for two groups, with limits
#' own <- data.frame(
#'   group = rep(c("Nation A", "Everyone else"), each = 4),
#'   rating = rep(c("Poor", "Fair", "Good", "Excellent"), 2),
#'   pct = c(10, 25, 40, 25, 6, 20, 44, 30),
#'   lo = c(7, 21, 35, 21, 5, 18, 42, 28),
#'   hi = c(13, 29, 45, 29, 7, 22, 46, 32)
#' )
#' ordcmp(own, "group", "rating", "pct", "lo", "hi",
#'        lvls = c("Poor", "Fair", "Good", "Excellent"), ref = "Everyone else")
#' @export
ordcmp <- function(
  data,
  grp,
  cat,
  est,
  lower,
  upper,
  lvls,
  ref,
  qual = NULL,
  conflvl = 0.95,
  pooled = NULL
) {
  .igc(
    "ordcmp",
    data = data,
    grp = grp,
    cat = cat,
    est = est,
    lower = lower,
    upper = upper,
    lvls = lvls,
    ref = ref,
    qual = qual,
    conflvl = conflvl,
    pooled = pooled
  )
}

#' @export
format.ordcmp <- function(x, ...) {
  need <- c("group", "reference", "probbet", "lower", "upper", "delta")
  fb <- .igplain(x, need, "Ordinal comparison against a reference group")
  if (!is.null(fb)) {
    return(fb)
  }
  lab <- .igpool(.igstar(x$group, x$usewithcau), x$poogrp)
  show <- data.frame(
    group = lab,
    probbet = x$probbet,
    interval = .igci(x$lower, x$upper),
    delta = x$delta,
    neff = round(x$neff)
  )
  sh <- attr(x, "shares")
  .igb(
    "Ordinal comparison against a reference group",
    facts = list(
      Reference = x$reference[1L],
      `Reference neff` = round(x$neffref[1L]),
      Levels = if (is.matrix(sh)) ncol(sh),
      Measure = "chance a person in the group rates higher, ties half"
    ),
    tables = list(show),
    reading = .ignull(
      lab,
      x$lower,
      x$upper,
      0.5,
      sprintf("Rate higher than %s", x$reference[1L]),
      sprintf("Rate lower than %s", x$reference[1L])
    ),
    notes = c(
      "probbet 0.5 is no difference; above 0.5 the group rates higher.",
      "delta (Cliff's) = share rating higher less share rating lower.",
      "neff: the effective sample size the published limits imply (it",
      "carries the survey's design effect).",
      if (any(.igyes(x$poogrp))) {
        paste(
          "(pooled): a pan-Indigenous total, shown beside, not a",
          "distinct group (poogrp)."
        )
      },
      .igcau(x$usewithcau)
    )
  )
}

#' Age-standardise survey percentages
#'
#' The percentage (or proportion, or rate) a group would have with the age
#' structure of a standard population, from the estimate and published limits in
#' each age group: for comparing First Nations, Metis and Inuit prevalences
#' (suicidal thoughts, mood and anxiety disorders, heavy drinking) with each
#' other and with an older population, for example, when the age pattern of the
#' outcome is steep.
#'
#' @details With estimates \eqn{p_i}, standard errors \eqn{s_i = (u_i - l_i) /
#' (2z)} and standard weights \eqn{w_i} (the standard's shares), the
#' standardised estimate is \eqn{\sum w_i p_i} with standard error
#' \eqn{\sqrt{\sum w_i^2 s_i^2}}, the age groups being disjoint domains of the
#' survey (treated as independent: see \code{\link{indgap}}). The interval is
#' normal, cut at 0 and at \code{scale}. A group with a suppressed or missing
#' age group, or one whose published interval misses its own estimate, has no
#' standardised value: the age group is not imputed. The standard should be the
#' age structure of the population compared, or a common one: for First Nations,
#' Metis and Inuit comparisons, the Indigenous population's own
#' (\code{\link{cnsstd}}; Kumar and Tjepkema 2019).
#'
#' @param est,lower,upper Numeric vectors (not column names): the estimate in
#'   each age group (and group) and its published limits.
#' @param std The standard population, one count (or share) per age group.
#' @param grp Optional labels, one per row: one result per group. Rows are
#'   matched to \code{std} by position, so order each group's rows by age in the
#'   order of \code{std} (sort by an age order, not alphabetically).
#' @param qual Optional quality flags, one per row (a vector; see
#'   \code{\link{qltyrcd}}).
#' @param conflvl Confidence level of the prevalence intervals (default 0.95).
#' @param scale The upper bound of the scale: 100 for percentages (the default),
#'   1 for proportions, \code{Inf} for rates.
#' @return A data frame of class \code{stdprev}: \code{group},
#'   \code{estimate}, \code{lower}, \code{upper}, \code{se}, \code{usewithcau};
#'   attribute \code{conflvl}.
#' @references Klein RJ, Schoenborn CA (2001). Age adjustment using the 2000
#' projected U.S. population. Healthy People 2010 Statistical Notes 20. National
#' Center for Health Statistics.
#'
#' Kumar MB, Tjepkema M (2019). Suicide among First Nations people, Metis and
#' Inuit (2011-2016): findings from the 2011 Canadian Census Health and
#' Environment Cohort (CanCHEC). National Household Survey: Aboriginal Peoples,
#' Statistics Canada catalogue no. 99-011-X2019001.
#' @seealso \code{\link{dastd}} for counts of events,
#'   \code{\link{cnsstd}}, \code{\link{indgap}}.
#' @examples
#' # made-up percentages in three age groups for two groups, standardised to a
#' # population of 500, 300 and 200 in those age groups
#' stdprev(
#'   est = c(12, 9, 6, 8, 6, 4),
#'   lower = c(9, 6, 3, 7, 5, 3),
#'   upper = c(15, 12, 9, 9, 7, 5),
#'   std = c(500, 300, 200),
#'   grp = rep(c("Nation A", "Everyone else"), each = 3)
#' )
#' @export
stdprev <- function(
  est,
  lower,
  upper,
  std,
  grp = NULL,
  qual = NULL,
  conflvl = 0.95,
  scale = 100
) {
  .igc(
    "stdprev",
    est = est,
    lower = lower,
    upper = upper,
    std = std,
    grp = grp,
    qual = qual,
    conflvl = conflvl,
    scale = scale
  )
}

#' @export
format.stdprev <- function(x, ...) {
  fb <- .igplain(
    x,
    c("group", "estimate", "lower", "upper"),
    "Age-standardised prevalence"
  )
  if (!is.null(fb)) {
    return(fb)
  }
  lab <- .igstar(x$group, x$usewithcau)
  show <- data.frame(
    group = lab,
    estimate = x$estimate,
    interval = .igci(x$lower, x$upper),
    se = x$se
  )
  .igb(
    "Age-standardised prevalence",
    facts = list(
      Groups = nrow(x),
      Method = "direct standardisation to the standard population given"
    ),
    tables = list(show),
    reading = .igrange(lab, x$estimate, "prevalence", x$lower, x$upper),
    notes = c(
      "Each group's prevalence as if it had the standard population's ages,",
      "so the groups compare without their age structures.",
      "Compare only rates standardised to the same population.",
      .igcau(x$usewithcau)
    )
  )
}

#' Wellness beside distress: a strengths-first profile of published indicators
#'
#' Puts a group's wellness indicators (life satisfaction, perceived mental
#' health, a hopeful outlook, a sense of belonging, a strong cultural identity)
#' first and its distress indicators (suicidal thoughts, mood and anxiety
#' disorders, fair or poor mental health) after them, each with its published
#' interval, for each First Nations, Metis and Inuit group: so that a table of
#' mental health is never read, or plotted, as a list of deficits alone.
#'
#' @details Indigenous scholars of mental health define it as wellness and as a
#' process of healing, grounded in community, cultural identity, a holistic view
#' and interdependence (Stewart 2008), and ask that strengths be reported
#' together with disparities (Stewart 2020). This function puts that into how
#' published estimates are reported; it adds no estimate of its own. The
#' optional \code{domain} groups indicators by the holistic domains of the
#' Medicine Wheel teachings that Gabriel, Mirza and Stewart (2022) describe
#' (physical, emotional, mental, spiritual). That grouping is a way to lay out
#' the table, not a validated scale, and published tables have few indicators of
#' the spiritual domain: the profile says what is there and does not stand one
#' indicator in for another.
#'
#' Aggregate tables show how indicators vary together across groups; they cannot
#' show that culture, land or belonging protects mental health within a group,
#' which would be an ecological inference. The differences these tables show
#' between groups are the effects of colonial policy and of structural
#' conditions, not of the people in the groups (Gabriel et al. 2022); their
#' interpretation belongs with First Nations, Metis and Inuit communities.
#'
#' @param data A data frame, one row per group and indicator (as
#'   \code{\link{scests}} gives).
#' @param grp,ind,est,lower,upper Column names: the group, the indicator,
#'   the percentage and its published limits.
#' @param well,distress The indicators that measure wellness or strength, and
#'   those that measure distress; rows for any other indicator (totals,
#'   \dQuote{not specified}) are left out.
#' @param domain Optional named character vector giving each indicator's
#'   holistic domain, for example \code{c("Life satisfaction 8 to 10" =
#'   "mental", ...)}.
#' @param qual Optional column of quality flags, one per row of `data`.
#' @return A data frame of class \code{welprof}, wellness rows
#'   first, with the columns \code{group}, \code{indicator}, \code{side}
#'   (\code{"wellness"} or \code{"distress"}), \code{domain}, \code{estimate},
#'   \code{lower}, \code{upper} and \code{usewithcau};
#' @section Acknowledgement: This estimator puts to quantitative use the
#' framework of Professor Suzanne Stewart (Yellowknives Dene First Nation,
#' University of Toronto).
#'
#' See the acknowledgements in the section on Indigenous health of the
#' package's vignette (\code{vignette("crimlab")}).
#' @references Stewart SL (2008). Promoting Indigenous mental health: cultural
#' perspectives on healing from Native counsellors in Canada. International
#' Journal of Health Promotion and Education 46(2):49-56.
#' \doi{10.1080/14635240.2008.10708129}
#'
#' Stewart S (2020). Addressing strengths and disparities in Indigenous health.
#' International Journal of Indigenous Health 15(1):4.
#' \doi{10.32799/ijih.v15i1.34976}
#'
#' Gabriel MD, Mirza S, Stewart SL (2022). Exploring mental health and holistic
#' healing through the life stories of Indigenous youth who have experienced
#' homelessness. International Journal of Environmental Research and Public
#' Health 19(20):13402. \doi{10.3390/ijerph192013402}
#' @seealso \code{\link{scests}}, \code{\link{indgap}},
#'   \code{\link{ordcmp}}.
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
#' welprof(own, "community", "item", "pct", "lo", "hi",
#'         well = c("Hopeful", "Satisfied"),
#'         distress = c("Rarely hopeful", "Low satisfaction"))
#' @export
welprof <- function(
  data,
  grp,
  ind,
  est,
  lower,
  upper,
  well,
  distress,
  domain = NULL,
  qual = NULL
) {
  .igc(
    "welprof",
    data = data,
    grp = grp,
    ind = ind,
    est = est,
    lower = lower,
    upper = upper,
    well = well,
    distress = distress,
    domain = domain,
    qual = qual
  )
}

#' @export
format.welprof <- function(x, ...) {
  need <- c("group", "indicator", "side", "estimate", "lower", "upper")
  fb <- .igplain(x, need, "Wellness profile")
  if (!is.null(fb)) {
    return(fb)
  }
  ind <- unique(x$indicator)
  .igb(
    "Wellness profile",
    facts = list(
      Groups = paste(unique(x$group), collapse = "; "),
      `Wellness indicators` = length(unique(x$indicator[x$side == "wellness"])),
      `Distress indicators` = length(unique(x$indicator[x$side == "distress"]))
    ),
    tables = .wptabs(x, ind),
    reading = .wpread(x, ind),
    notes = c(
      "Wellness indicators come first, so each group's strengths are read",
      "before its distress (Stewart 2008).",
      .igcau(x$usewithcau)
    )
  )
}

#' Holistic balance: how evenly wellness is held across the domains of a whole
#' person
#'
#' For each group, from published percentages of wellness in each holistic
#' domain (physical, emotional, mental, spiritual, or the domains a community
#' names), the overall level of wellness and how balanced it is across the
#' domains, with Monte Carlo intervals from the published ones.
#'
#' @details Indigenous understandings of mental health are holistic: wellness is
#' a balance among the domains of a person, and illness a loss of that balance
#' (Stewart 2008; Gabriel, Mirza and Stewart 2022). A single indicator, or an
#' average of indicators, cannot show balance: a group with 80% in two domains
#' and 20% in two others has the same average as one with 50% in all four. This
#' method, which is ours and inspired by that work rather than taken from it,
#' measures balance as the ratio of the weighted geometric to the weighted
#' arithmetic mean of the domain levels \eqn{p_k}: \deqn{B = \exp\left(\sum_k
#' w_k \log p_k\right) / \sum_k w_k p_k,} one less Atkinson's (1970) index with
#' inequality aversion 1. \eqn{B = 1} exactly when every domain with a positive
#' weight is at the same level, and it falls toward 0 as one domain falls behind
#' the others; it does not depend on the overall level, which is reported beside
#' it (\code{level}, \eqn{\sum_k w_k p_k}) with the weakest domain.
#'
#' The weights \code{w} say how much each domain counts. They are equal by
#' default; a community that holds its own view of the domains sets its own (the
#' CARE principle of authority to control, Carroll et al. 2020), and the result
#' says which weights it used.
#'
#' Several indicators in one domain are averaged, their standard errors combined
#' as independent (they come from the same respondents, so that standard error
#' is too small). The interval for \eqn{B} comes from \code{reps} draws of every
#' domain level from a normal distribution on the logit scale with the standard
#' error of its published interval, the domains independent. \eqn{B} is at most
#' 1 and reaches 1 only at equal levels, so noise in the levels tends to pull
#' the draws below \eqn{\hat B} (\eqn{B} is not concave; the direction of the
#' bias is seen in the draws, not guaranteed): the percentile interval is
#' shifted back by the draws' mean bias (\eqn{\bar B^* - \hat B}), kept inside
#' \eqn{(0, 1]} and widened if need be to contain \eqn{\hat B}. The level's
#' interval is normal. A group with a suppressed or missing domain, or one whose
#' published interval misses its own estimate, has no balance: the domain is not
#' imputed.
#'
#' What the domains mean, and which indicators stand for each, is for the
#' communities to say; the published tables have few measures of spiritual
#' wellness, and a domain with no fitting indicator should be left out rather
#' than filled with another.
#'
#' @param grp,domain,est,lower,upper Column names: the group, the domain of the
#'   indicator (any set of domains: Medicine Wheel, Te Whare Tapa Wha, a
#'   Nation's own), the percentage of wellness (0 to 100, strictly inside;
#'   multiply proportions by 100) and its published limits.
#' @param w Optional named weights, one per domain (they are rescaled to sum to
#'   one).
#' @param qual Optional column of quality flags for the estimates.
#' @param reps Monte Carlo draws for the intervals (default 10,000).
#' @param seed Seed for the draws, so the intervals are reproducible (default
#'   1).
#' @param conflvl Confidence level of the balance intervals (default 0.95).
#' @return A data frame of class \code{holbal}: \code{group},
#'   \code{level}, \code{lvllo}, \code{lvlup}, \code{balance}, \code{ballo},
#'   \code{balup}, \code{weadom}, \code{usewithcau}; the attributes
#'   \code{domains} (the matrix of domain levels, groups by domains) and
#'   \code{weights}.
#' @section Acknowledgement: This estimator puts to quantitative use the
#' framework of Professor Suzanne Stewart (Yellowknives Dene First Nation,
#' University of Toronto).
#'
#' See the acknowledgements in the section on Indigenous health of the
#' package's vignette (\code{vignette("crimlab")}).
#' @references Atkinson AB (1970). On the measurement of inequality. Journal of
#' Economic Theory 2(3):244-263. \doi{10.1016/0022-0531(70)90039-6}
#'
#' Carroll SR, Garba I, Figueroa-Rodriguez OL, et al. (2020). The CARE
#' Principles for Indigenous Data Governance. Data Science Journal 19:43.
#' \doi{10.5334/dsj-2020-043}
#'
#' Gabriel MD, Mirza S, Stewart SL (2022). Exploring mental health and holistic
#' healing through the life stories of Indigenous youth who have experienced
#' homelessness. International Journal of Environmental Research and Public
#' Health 19(20):13402. \doi{10.3390/ijerph192013402}
#'
#' Stewart SL (2008). Promoting Indigenous mental health: cultural perspectives
#' on healing from Native counsellors in Canada. International Journal of Health
#' Promotion and Education 46(2):49-56. \doi{10.1080/14635240.2008.10708129}
#' @seealso \code{\link{welprof}}, \code{\link{dualcont}}.
#' @examples
#' # your own table: made-up percentages of wellness in four domains a Nation
#' # set for itself, for two communities
#' own <- data.frame(
#'   community = rep(c("Lakeside", "Riverbend"), each = 4),
#'   domain = rep(c("body", "mind", "spirit", "family"), 2),
#'   pct = c(62, 58, 70, 75, 55, 61, 66, 80),
#'   lo = c(57, 53, 65, 70, 49, 55, 60, 75),
#'   hi = c(67, 63, 75, 80, 61, 67, 72, 85)
#' )
#' holbal(own, "community", "domain", "pct", "lo", "hi", reps = 2000)
#' @export
#' @inheritParams welprof
holbal <- function(
  data,
  grp,
  domain,
  est,
  lower,
  upper,
  w = NULL,
  qual = NULL,
  reps = 10000,
  seed = 1,
  conflvl = 0.95
) {
  .igc(
    "holbal",
    data = data,
    grp = grp,
    domain = domain,
    est = est,
    lower = lower,
    upper = upper,
    w = w,
    qual = qual,
    reps = reps,
    seed = seed,
    conflvl = conflvl
  )
}

#' @export
format.holbal <- function(x, ...) {
  need <- c("group", "level", "lvllo", "lvlup", "balance", "ballo", "balup")
  fb <- .igplain(x, need, "Holistic balance across the domains")
  if (!is.null(fb)) {
    return(fb)
  }
  lab <- .igstar(x$group, x$usewithcau)
  w <- attr(x, "weights")
  tb <- sort(table(x$weadom), decreasing = TRUE)
  .igb(
    "Holistic balance across the domains",
    facts = list(
      Groups = nrow(x),
      Weights = .hbwts(w),
      Balance = "weighted geometric over arithmetic mean of the domain levels"
    ),
    tables = .hbtabs(x, lab),
    reading = c(
      .igrange(lab, x$balance, "balance", x$ballo, x$balup),
      if (length(tb)) {
        sprintf(
          "The weakest domain is %s for %d of %d groups.",
          names(tb)[1L],
          as.integer(tb[1L]),
          nrow(x)
        )
      }
    ),
    notes = c(
      "balance is 1 when every domain is at the same level and falls toward",
      "0 as one falls behind; level is the weighted mean, reported beside it.",
      "A method inspired by Stewart (2008), not taken from it.",
      "The weights say how much each domain counts; a community may set",
      "its own (CARE: authority to control, Carroll et al. 2020).",
      .igcau(x$usewithcau)
    )
  )
}

#' Wellness and distress as two continua
#'
#' For each group, one wellness indicator against one distress indicator, each
#' with its published interval, and the ratio of the two: the dual-continua view
#' of mental health, in which wellness is not the absence of distress and is
#' measured in its own right.
#'
#' @details Keyes (2002, 2005) showed that flourishing and mental illness are
#' correlated but distinct axes: a person, or a population, can carry distress
#' and also wellness. That is how Indigenous scholars of mental health describe
#' it too, as wellness and as a process of healing rather than as illness alone
#' (Stewart 2008). Plotting the groups on both axes keeps a strength visible
#' beside the burden it sits next to. The \code{strratio} is \eqn{R = W / D},
#' wellness over distress, with the log-normal interval \eqn{\exp\{\log R \pm
#' z\sqrt{\sigma_W^2 + \sigma_D^2}\}} from the two published intervals: above 1,
#' more people report the wellness than the distress. The two estimates come
#' from the same respondents; their covariance is not published, so the interval
#' treats them as independent and is too narrow when they are negatively
#' correlated, as answers about wellness and distress usually are. Any two
#' indicators on one scale can be set against each other.
#'
#' @param x A wellness profile from \code{\link{welprof}}.
#' @param well,distress One wellness and one distress indicator of \code{x}.
#' @param conflvl Confidence level of the published intervals and of the ratio
#'   (default 0.95).
#' @return A data frame of class \code{dualcont} with \code{group},
#'   \code{wellness}, \code{wello}, \code{welup}, \code{distress},
#'   \code{distrlo}, \code{distrup}, \code{strratio},
#'   \code{ratiolo}, \code{ratioup} and \code{usewithcau}.
#' @section Acknowledgement: This estimator puts to quantitative use the
#' framework of Professor Suzanne Stewart (Yellowknives Dene First Nation,
#' University of Toronto).
#'
#' See the acknowledgements in the section on Indigenous health of the
#' package's vignette (\code{vignette("crimlab")}).
#' @references Keyes CLM (2002). The mental health continuum: from languishing
#' to flourishing in life. Journal of Health and Social Behavior 43(2):207-222.
#' \doi{10.2307/3090197}
#'
#' Keyes CLM (2005). Mental illness and/or mental health? Investigating axioms
#' of the complete state model of health. Journal of Consulting and Clinical
#' Psychology 73(3):539-548. \doi{10.1037/0022-006X.73.3.539}
#'
#' Stewart SL (2008), see \code{\link{holbal}}.
#' @seealso \code{\link{welprof}}.
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
#' w <- welprof(own, "community", "item", "pct", "lo", "hi",
#'              well = "Hopeful", distress = "Rarely hopeful")
#' dualcont(w, "Hopeful", "Rarely hopeful")
#' @export
dualcont <- function(x, well, distress, conflvl = 0.95) {
  .igc("dualcont", x = x, well = well, distress = distress, conflvl = conflvl)
}

#' @export
format.dualcont <- function(x, ...) {
  need <- c("group", "wellness", "wello", "distress", "strratio", "ratiolo")
  fb <- .igplain(x, need, "Dual continuum: wellness beside distress")
  if (!is.null(fb)) {
    return(fb)
  }
  lab <- attr(x, "labels")
  show <- data.frame(
    group = .igstar(x$group, x$usewithcau),
    wellness = x$wellness,
    wi = .igci(x$wello, x$welup),
    distress = x$distress,
    di = .igci(x$distrlo, x$distrup),
    ratio = x$strratio,
    ri = .igci(x$ratiolo, x$ratioup)
  )
  names(show)[c(3L, 5L, 7L)] <- c("interval", "interval ", "ratio interval")
  .igb(
    "Dual continuum: wellness beside distress",
    facts = list(
      Wellness = if (!is.null(lab)) lab[["wellness"]],
      Distress = if (!is.null(lab)) lab[["distress"]],
      Groups = nrow(x),
      Ratio = "wellness over distress, log-normal interval"
    ),
    tables = list(show),
    reading = c(
      "A ratio above 1: more people report the wellness than the distress.",
      sprintf(
        "%d of %d groups have the whole interval above 1.",
        sum(x$ratiolo > 1, na.rm = TRUE),
        nrow(x)
      )
    ),
    notes = c(
      "Both estimates come from the same respondents and their covariance",
      "is not published: the interval treats them as independent.",
      .igcau(x$usewithcau)
    )
  )
}

#' format.welprof(): one table per indicator, named with its side
#' @noRd
.wptabs <- function(x, ind) {
  tabs <- lapply(ind, function(k) {
    r <- x[x$indicator == k, , drop = FALSE]
    out <- data.frame(
      group = .igstar(r$group, r$usewithcau),
      estimate = r$estimate,
      interval = .igci(r$lower, r$upper)
    )
    if (!all(is.na(r$domain))) {
      out$domain <- ifelse(is.na(r$domain), "", as.character(r$domain))
    }
    out
  })
  side <- x$side[match(ind, x$indicator)]
  names(tabs) <- sprintf("%s (%s)", ind, side)
  tabs
}

#' format.welprof(): on which indicators the groups' intervals part
#' @noRd
.wpread <- function(x, ind) {
  apart <- vapply(
    ind,
    function(k) {
      g <- x[x$indicator == k & is.finite(x$estimate), , drop = FALSE]
      nrow(g) > 1L &&
        isTRUE(g$lower[which.max(g$estimate)] > g$upper[which.min(g$estimate)])
    },
    NA
  )
  if (any(apart)) {
    sprintf(
      "The groups' intervals part on %d of %d indicators: %s.",
      sum(apart),
      length(ind),
      paste(ind[apart], collapse = "; ")
    )
  } else {
    sprintf(
      "On none of the %d indicators do the groups' intervals part: %s",
      length(ind),
      "the data cannot rank the groups on any of them."
    )
  }
}

#' format.holbal(): the balance table, with the domain levels when kept
#' @noRd
.hbtabs <- function(x, lab) {
  show <- data.frame(
    group = lab,
    level = x$level,
    lci = .igci(x$lvllo, x$lvlup),
    balance = x$balance,
    bci = .igci(x$ballo, x$balup),
    weakest = x$weadom
  )
  names(show)[c(3L, 5L)] <- c("level CI", "balance CI")
  tabs <- list(show)
  m <- attr(x, "domains")
  if (is.matrix(m) && nrow(m) == nrow(x)) {
    tabs[["Domain levels (percent)"]] <- cbind(
      data.frame(group = x$group),
      as.data.frame(unclass(m), row.names = NULL, optional = TRUE)
    )
  }
  tabs
}

#' format.holbal(): the weights, "equal" when they are
#' @noRd
.hbwts <- function(w) {
  if (length(w)) {
    if (length(unique(w)) == 1L) {
      sprintf("equal, over %d domains", length(w))
    } else {
      toString(sprintf("%s %s", names(w), format(signif(w, 3))))
    }
  }
}
