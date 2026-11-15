# The 2021 Census Indigenous Population Profile (Statistics Canada
# 98-510-X2021001 to -016): read as published, shares with intervals that
# carry the long-form sample and random rounding, an Indigenous standard
# population for age standardisation, and a Community Well-Being score from
# public figures; and the Truth and Reconciliation Commission's Calls to
# Action 19 (health gaps) and 55 (annual reporting), measured from the
# published tables. The computation is in src/igc.cpp.

#' Read the 2021 Census Indigenous Population Profile
#'
#' Reads a table of the Indigenous Population Profile, 2021 Census of Population
#' (Statistics Canada catalogue 98-510-X2021001 to 98-510-X2021016), as
#' published, into one row per geography, gender, age group,
#' characteristic and group: First Nations, Metis, Inuit, the other identity
#' groups and non-Indigenous people (or, in the household tables, Indigenous and
#' non-Indigenous households), each in its own rows.
#'
#' @details The profile gives each characteristic in a hierarchy shown only by
#' indentation (\dQuote{Total - Highest certificate, diploma or degree},
#' \dQuote{High (secondary) school diploma or equivalency certificate}, ...):
#' the level, the parent and the topic's total (\code{rootid}) are worked out
#' from it, so a share can be taken of its parent. A topic marked "25%
#' sample data" comes from the long-form sample (\code{sample = "25\%"}); the
#' others from the full count (\code{"100\%"}). Counts are randomly rounded to
#' multiples of 5 (of 10 below 10) by Statistics Canada, and a cell suppressed
#' to protect confidentiality (\code{x}), too unreliable (\code{F}), not
#' available (\code{..}) or not applicable (\code{...}) has no value here;
#' \code{E} values are kept and marked. The global non-response rates of the
#' geography (\code{tnrsf}, \code{tnrlf}, in percent) are kept: Statistics
#' Canada suppresses an area whose long-form rate is 50% or more, and a high
#' rate means a larger error than the sample alone.
#'
#' @param x A profile table as published (a data frame read from its CSV).
#' @return A data frame of class \code{cnsprof}: \code{geo},
#'   \code{dguid}, \code{geolvl}, \code{residence} (on or off reserve, in the
#'   tables that split it), \code{gender}, \code{age}, \code{chaid},
#'   \code{characteristic}, \code{level}, \code{parid}, \code{rootid},
#'   \code{sample}, \code{group}, \code{value}, \code{quality}, \code{tnrsf},
#'   \code{tnrlf}.
#' @references Statistics Canada (2023). Indigenous Population Profile, 2021
#' Census of Population. Statistics Canada Catalogue no. 98-510-X2021001.
#' Released June 21, 2023. Published at
#' www12.statcan.gc.ca/census-recensement/2021/dp-pd/ipp-ppa/index.cfm?Lang=E
#' (the site answers scripts with 403); the Wayback Machine holds a copy of that
#' address captured on 8 March 2026 (web.archive.org, snapshot 20260308015010).
#'
#' Statistics Canada (2022). Guide to the Census of Population, 2021. Catalogue
#' no. 98-304-X.
#' @seealso \code{\link{profshr}}, \code{\link{cnsstd}},
#'   \code{\link{cnscwb}}.
#' @examples
#' # made-up numbers in the layout of a population profile table, for one
#' # place and two groups (the published tables hold many more rows)
#' nm <- c(
#'   "Total - Age groups of the population - 100% data",
#'   paste("Number of total income recipients aged 15 years and over in",
#'         "private households in 2020 - 25% sample data"),
#'   "Average total income in 2020 among recipients ($)",
#'   paste("Total - Highest certificate, diploma or degree for the population",
#'         "aged 25 years and over in private households - 25% sample data"),
#'   "  No certificate, diploma or degree",
#'   "  Bachelor's degree or higher",
#'   "Participation rate",
#'   "Unemployment rate"
#' )
#' x <- data.frame(
#'   GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
#'   AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3),
#'                 rep("25 to 64 years", 2)),
#'   RESIDENCE = "Total - Residence on or off reserve",
#'   CHARACTERISTIC_ID = 1:8, CHARACTERISTIC_NAME = nm,
#'   `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14),
#'   SYMBOL_IND_IDENT = "",
#'   NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6),
#'   SYMBOL_NON_IND = "",
#'   check.names = FALSE
#' )
#' p <- cnsprof(x)
#' p[p$chaid == 5, c("group", "characteristic", "level", "parid", "value")]
#' @export
cnsprof <- function(x) .igc("cnsprof", x = x)

#' @export
format.cnsprof <- function(x, ...) {
  need <- c("geo", "dguid", "group", "chaid", "characteristic", "value")
  fb <- .igplain(x, need, "Indigenous Population Profile, tidied")
  if (!is.null(fb)) {
    return(fb)
  }
  keep <- c(
    "geo",
    "dguid",
    "geolvl",
    "residence",
    "gender",
    "age",
    "group",
    "chaid",
    "characteristic",
    "level",
    "parid",
    "rootid",
    "sample",
    "value",
    "quality",
    "tnrsf",
    "tnrlf"
  )
  n <- min(nrow(x), 20L)
  .igb(
    "Indigenous Population Profile, tidied",
    facts = list(
      Rows = nrow(x),
      Geographies = .igvals(x$geo),
      Groups = .igvals(x$group),
      Characteristics = length(unique(x$chaid)),
      Samples = .igvals(x$sample)
    ),
    tables = list(x[seq_len(n), intersect(keep, names(x)), drop = FALSE]),
    reading = .cpread(x),
    notes = c(
      if (nrow(x) > n) sprintf("(%d more rows)", nrow(x) - n),
      "One row per geography, group and characteristic: the value with its",
      "level and parent in the profile's hierarchy, ready for profshr().",
      "Counts are randomly rounded to multiples of 5 (10 below 10).",
      "Statistics Canada suppresses an area whose long-form non-response",
      "is 50% or more; a high rate means a larger error than the sample."
    )
  )
}

#' The share of a population with a characteristic, from the census profile
#'
#' The share of a group with a characteristic of the 2021 Census Indigenous
#' Population Profile (foster children among children, speakers of an Indigenous
#' language, households needing major repairs, ...), with an approximate
#' confidence interval that carries the long-form sample and the random rounding
#' of the counts.
#'
#' @details The share is the numerator's count over the denominator's (by
#' default the numerator's parent in the profile's hierarchy), for each
#' geography, gender, age group and group. Statistics Canada does not publish
#' standard errors for the profile, so the interval is an approximation built
#' from what is known about how the counts were made: \itemize{ \item a
#' long-form ("25% sample") count is estimated from a sample of about one
#' household in four, so the share has variance \eqn{d (1 - f) p (1 - p) / (f
#' N)} for a sampling fraction \eqn{f} (\code{smpfrac}, 0.25) of a population of
#' \eqn{N} (the denominator) and a design effect \eqn{d} (\code{deff}, 1 by
#' default: the clustering of people in households raises it, and a value of 1.5
#' to 2 is prudent for characteristics shared within a household). On reserves
#' (rows whose \code{residence} is \code{"On reserve"}) every household got the
#' long form, so \eqn{f = 1} there and the sampling term vanishes; \item every
#' count is randomly rounded, up or down without bias, to a multiple of \eqn{b}
#' (\code{roubase}, 5). Averaged over the remainders, that adds variance \eqn{c
#' = (b^2 - 1)/6} to each count (4 for \eqn{b = 5}), so \eqn{c (1 + p^2) / N^2}
#' to the share (the only error of a full-count share). } In all, \deqn{V =
#' \frac{d (1 - f) p (1 - p)}{f N} + \frac{c (1 + p^2)}{N^2}.}{V = d (1 - f) p
#' (1 - p) / (f N) + c (1 + p^2) / N^2.} The interval is Wilson's at the
#' effective sample size \eqn{n_e = \max(p(1-p), 1/N) / V} (the floor \eqn{1/N}
#' keeps a share of 0 or 100% from an infinite \eqn{n_e}), so it stays inside 0
#' and 100. It does not carry non-response (see \code{tnrlf}) or the edit and
#' imputation error, so it is narrower than the full error; treat it as a lower
#' bound on the uncertainty.
#'
#' @param x A \code{\link{cnsprof}}.
#' @param num The characteristic counted: its \code{chaid}, or its name
#'   exactly.
#' @param den The characteristic it is a share of (id or name); default its
#'   parent.
#' @param smpfrac The long-form sampling fraction (default 0.25); full-count
#'   characteristics use 1 whatever is given.
#' @param deff Design effect (default 1).
#' @param roubase The base of the random rounding (default 5, Statistics
#'   Canada's); 1 for counts that are not rounded.
#' @param conflvl Confidence level of the share's limits (default 0.95).
#' @return A data frame of class \code{profshr} with \code{geo},
#'   \code{dguid}, \code{residence}, \code{gender}, \code{age}, \code{group},
#'   \code{numerator}, \code{denominator}, \code{estimate} (percent),
#'   \code{lower}, \code{upper}, \code{quality}; ready for \code{\link{indgap}}.
#' @references Statistics Canada (2022). Guide to the Census of Population,
#' 2021, chapter 9 (data quality) and chapter 10 (confidentiality and random
#' rounding). Catalogue no. 98-304-X.
#'
#' Wilson EB (1927). Probable inference, the law of succession, and statistical
#' inference. Journal of the American Statistical Association 22(158):209-212.
#' \doi{10.1080/01621459.1927.10502953}
#' @seealso \code{\link{cnsprof}}.
#' @examples
#' # made-up numbers in the layout of a population profile table, for one
#' # place and two groups (the published tables hold many more rows)
#' nm <- c(
#'   "Total - Age groups of the population - 100% data",
#'   paste("Number of total income recipients aged 15 years and over in",
#'         "private households in 2020 - 25% sample data"),
#'   "Average total income in 2020 among recipients ($)",
#'   paste("Total - Highest certificate, diploma or degree for the population",
#'         "aged 25 years and over in private households - 25% sample data"),
#'   "  No certificate, diploma or degree",
#'   "  Bachelor's degree or higher",
#'   "Participation rate",
#'   "Unemployment rate"
#' )
#' x <- data.frame(
#'   GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
#'   AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3),
#'                 rep("25 to 64 years", 2)),
#'   RESIDENCE = "Total - Residence on or off reserve",
#'   CHARACTERISTIC_ID = 1:8, CHARACTERISTIC_NAME = nm,
#'   `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14),
#'   SYMBOL_IND_IDENT = "",
#'   NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6),
#'   SYMBOL_NON_IND = "",
#'   check.names = FALSE
#' )
#' p <- cnsprof(x)
#' # the share with no certificate, of the population aged 25 and over
#' profshr(p, num = "No certificate, diploma or degree")
#' @export
profshr <- function(
  x,
  num,
  den = NULL,
  smpfrac = 0.25,
  deff = 1,
  conflvl = 0.95,
  roubase = 5
) {
  .igc(
    "profshr",
    x = x,
    num = num,
    den = den,
    smpfrac = smpfrac,
    deff = deff,
    conflvl = conflvl,
    roubase = roubase
  )
}

#' @export
format.profshr <- function(x, ...) {
  need <- c("geo", "group", "numerator", "denominator", "estimate", "lower")
  fb <- .igplain(x, need, "Share of a population with a characteristic")
  if (!is.null(fb)) {
    return(fb)
  }
  ch <- attr(x, "characteristics")
  n <- min(nrow(x), 20L)
  h <- x[seq_len(n), , drop = FALSE]
  show <- data.frame(
    geo = h$geo,
    dguid = h$dguid,
    residence = h$residence,
    gender = h$gender,
    age = h$age,
    group = h$group,
    numerator = h$numerator,
    denominator = h$denominator,
    estimate = h$estimate,
    interval = .igci(h$lower, h$upper),
    quality = h$quality
  )
  .igb(
    "Share of a population with a characteristic",
    facts = list(
      Numerator = if (!is.null(ch)) ch[["numerator"]],
      Denominator = if (!is.null(ch)) ch[["denominator"]],
      Rows = nrow(x),
      Estimate = "percent"
    ),
    tables = list(show),
    reading = .igrange(.iglab(x), x$estimate, "share", x$lower, x$upper),
    notes = c(
      if (nrow(x) > n) sprintf("(%d more rows)", nrow(x) - n),
      "The interval carries the long-form sample and the random rounding,",
      "not non-response: read it as a lower bound on the uncertainty."
    )
  )
}

#' An Indigenous standard population from the census
#'
#' The age structure of a group in a geography of the 2021 Census Indigenous
#' Population Profile, as the standard of \code{\link{dastd}} or
#' \code{\link{stdprev}}.
#'
#' @details Comparisons of First Nations, Metis and Inuit health are
#' standardised to the Indigenous population's own, younger, age structure
#' (Kumar and Tjepkema 2019): a standard as old as Canada's gives weight to ages
#' where few Indigenous people are. This takes the population in private
#' households of each age group (the profile's first characteristic) for the
#' group and geography given. Age groups that do not overlap are chosen by
#' default (0 to 4, 5 to 9, 10 to 14, 15 to 19, 20 to 24, 25 to 34, 35 to 44, 45
#' to 54, 55 to 64, 65 to 74, 75 and over); collapse them to the age groups of
#' the rates being standardised.
#'
#' @param grp The group (default \code{"Indigenous identity"}).
#' @param geo The geography: its name (default \code{"Canada"}) or, where names
#'   repeat, its DGUID.
#' @param ages The age groups, as the profile names them.
#' @param resid On or off reserve, in the tables that split it (default the
#'   total).
#' @return A data frame: \code{age}, \code{population}, \code{weight} (shares
#'   summing to 1).
#' @references Kumar MB, Tjepkema M (2019). Suicide among First Nations people,
#'   Metis and Inuit (2011-2016): findings from the 2011 Canadian Census Health
#'   and Environment Cohort (CanCHEC). Statistics Canada catalogue no.
#'   99-011-X2019001.
#' @seealso \code{\link{dastd}}, \code{\link{stdprev}}.
#' @examples
#' # made-up population counts by age group for two groups, in the layout of a
#' # profile table
#' ages <- c("0 to 14 years", "15 to 64 years", "65 years and over")
#' x <- data.frame(
#'   GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
#'   AGE_GROUP = c("Total - Age groups", ages),
#'   RESIDENCE = "Total - Residence on or off reserve",
#'   CHARACTERISTIC_ID = 1,
#'   CHARACTERISTIC_NAME = "Total - Age groups of the population - 100% data",
#'   `IND_IDENT (2)` = c(1200, 400, 700, 100), SYMBOL_IND_IDENT = "",
#'   NON_IND = c(5000, 800, 3300, 900), SYMBOL_NON_IND = "",
#'   check.names = FALSE
#' )
#' cnsstd(cnsprof(x), geo = "Lakeside", ages = ages)
#' @export
#' @inheritParams profshr
cnsstd <- function(
  x,
  grp = "Indigenous identity",
  geo = "Canada",
  ages = c(
    "0 to 4 years",
    "5 to 9 years",
    "10 to 14 years",
    "15 to 19 years",
    "20 to 24 years",
    "25 to 34 years",
    "35 to 44 years",
    "45 to 54 years",
    "55 to 64 years",
    "65 to 74 years",
    "75 years and over"
  ),
  resid = "Total - Residence on or off reserve"
) {
  .igc("cnsstd", x = x, grp = grp, geo = geo, ages = ages, resid = resid)
}

#' A Community Well-Being score from the census profile
#'
#' Indigenous Services Canada's Community Well-Being (CWB) index
#' (\code{\link{cwbscr}}) computed from the public 2021 Census Indigenous
#' Population Profile, for each geography and group.
#'
#' @details The Department computes the CWB from census microdata; the public
#' profile gives close but not the same inputs, and this function says where it
#' departs: \describe{ \item{income}{per capita: the number of total income
#' recipients aged 15 and over times their average total income, over the
#' population in private households (all ages), as the CWB does.}
#' \item{education}{the shares with at least a high school diploma (all but
#' those with no certificate, diploma or degree) and with a bachelor's degree or
#' higher, both among those aged 25 and over: the profile has no age group of 20
#' and over, which the CWB uses for the first.} \item{labour force}{the
#' participation rate and the employment rate of participants (100 less the
#' unemployment rate) of those aged 25 to 64: the CWB uses 20 to 64.}
#' \item{housing}{the shares of households with one person or fewer per room and
#' needing no major repairs, from the household profile (\code{dwell}), which
#' splits Indigenous from non-Indigenous households only, so First Nations,
#' Metis and Inuit people all get the pooled Indigenous households' housing: the
#' CWB counts persons, not households. Without \code{dwell} the housing
#' component and the score are \code{NA}.} } A group of fewer than \code{minpop}
#' people in a geography gets no score. The Department scores communities of 65
#' people or more; the default here is 250, because the public profile's counts
#' are rounded to multiples of 5 and its long-form counts rest on a few dozen
#' households at that size. Geographies are matched by DGUID, so places that
#' share a name stay apart. The inputs are those of the 2021 profile
#' (98-510-X2021001 to -016): other editions name their characteristics
#' differently. Read the result as a CWB-style description of socio-economic
#' conditions, close to but not the Department's index, and never as a measure
#' of the wellness of First Nations, Metis or Inuit people as they define it
#' (see \code{\link{cwbscr}}).
#'
#' @param x A \code{\link{cnsprof}} of a population table.
#' @param dwell Optional \code{\link{cnsprof}} of the matching household table.
#' @param groups The groups to score (default First Nations, Metis, Inuit and
#'   non-Indigenous people).
#' @param bounds Income bounds of \code{\link{cwbscr}}.
#' @param minpop The smallest population scored (default 250; the Department
#'   uses 65).
#' @return A data frame of class \code{cnscwb}: \code{geo},
#'   \code{dguid}, \code{group}; the components \code{income}, \code{education},
#'   \code{housing}, \code{labour} and the score \code{cwb} (0 to 100, \code{NA}
#'   where not scored); and the inputs \code{incpercap} (dollars),
#'   \code{highsch25p}, \code{uni25p}, \code{population}, \code{participation},
#'   \code{employment}, \code{nomajrep}, \code{notcro} (percentages).
#' @references Indigenous Services Canada (2019). Report on trends in First
#'   Nations communities, 1981 to 2016 (Community Well-Being index).\cr
#'   sac-isc.gc.ca/eng/1345816651029/1557323327644
#' @seealso \code{\link{cwbscr}}.
#' @examples
#' # made-up numbers in the layout of a population profile table, for one
#' # place and two groups (the published tables hold many more rows)
#' nm <- c(
#'   "Total - Age groups of the population - 100% data",
#'   paste("Number of total income recipients aged 15 years and over in",
#'         "private households in 2020 - 25% sample data"),
#'   "Average total income in 2020 among recipients ($)",
#'   paste("Total - Highest certificate, diploma or degree for the population",
#'         "aged 25 years and over in private households - 25% sample data"),
#'   "  No certificate, diploma or degree",
#'   "  Bachelor's degree or higher",
#'   "Participation rate",
#'   "Unemployment rate"
#' )
#' x <- data.frame(
#'   GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
#'   AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3),
#'                 rep("25 to 64 years", 2)),
#'   RESIDENCE = "Total - Residence on or off reserve",
#'   CHARACTERISTIC_ID = 1:8, CHARACTERISTIC_NAME = nm,
#'   `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14),
#'   SYMBOL_IND_IDENT = "",
#'   NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6),
#'   SYMBOL_NON_IND = "",
#'   check.names = FALSE
#' )
#' p <- cnsprof(x)
#' cnscwb(p, groups = c("Indigenous identity", "Non-Indigenous"), minpop = 0)
#' @export
#' @inheritParams cnsstd
cnscwb <- function(
  x,
  dwell = NULL,
  groups = c("First Nations", "M\u00e9tis", "Inuit", "Non-Indigenous"),
  bounds = c(2650, 75000),
  resid = "Total - Residence on or off reserve",
  minpop = 250
) {
  .igc(
    "cnscwb",
    x = x,
    dwell = dwell,
    groups = groups,
    bounds = bounds,
    resid = resid,
    minpop = minpop
  )
}

#' @export
format.cnscwb <- function(x, ...) {
  need <- c("geo", "group", "income", "education", "housing", "labour", "cwb")
  fb <- .igplain(x, need, "Community Well-Being index, from the census")
  if (!is.null(fb)) {
    return(fb)
  }
  n <- min(nrow(x), 20L)
  .igb(
    "Community Well-Being index, from the census",
    facts = list(
      Rows = nrow(x),
      `Not scored` = if (anyNA(x$cwb)) sum(is.na(x$cwb)),
      Scale = "0 to 100 for each component and the score"
    ),
    tables = list(x[seq_len(n), need, drop = FALSE]),
    reading = .igrange(paste(x$geo, x$group, sep = ", "), x$cwb, "score"),
    notes = c(
      if (nrow(x) > n) sprintf("(%d more rows)", nrow(x) - n),
      "cwb is the mean of income, education, housing and labour force",
      "activity (Indigenous Services Canada's Community Well-Being index).",
      "CWB-style, close to but not the Department's index: the public",
      "profile's inputs differ (ages 25 and over; housing from households,",
      "pooled for First Nations, M\u00e9tis and Inuit people).",
      "Never a measure of wellness as First Nations, M\u00e9tis or Inuit",
      "people define it. A group below minpop people gets no score.",
      if (anyNA(x$housing)) "Housing needs the household table (dwell)."
    )
  )
}

#' The indicators of Calls to Action 19 and 55, and the tables that measure them
#'
#' Call to Action 19 asks the federal government to \dQuote{establish measurable
#' goals to identify and close the gaps in health outcomes between Aboriginal
#' and non-Aboriginal communities, and to publish annual progress reports and
#' assess long-term trends}, focusing on \dQuote{infant mortality, maternal
#' health, suicide, mental health, addictions, life expectancy, birth rates,
#' infant and child health issues, chronic diseases, illness and injury
#' incidence, and the availability of appropriate health services}. Call to
#' Action 55 asks for annual data on, among others, \dQuote{the number of
#' Aboriginal children--including Metis and Inuit children--in care, compared
#' with non-Aboriginal children} and \dQuote{the educational and income
#' attainments of Aboriginal peoples in Canada compared with non-Aboriginal
#' people}. This lists, for each indicator, the Statistics Canada tables
#' (or the census profile) that measure it, the measure, and which
#' direction is better.
#'
#' @details The list is a map of what the published tables can measure, not a
#' judgement of what the Calls require: maternal health is measured here only
#' through birth outcomes in Inuit regions, and several indicators (addictions,
#' the availability of services) only through proxies. The direction is the
#' direction of better health for the people measured.
#'
#' @return A data frame with \code{call}, \code{indicator} (the Call's words),
#'   \code{table}, \code{measure}, \code{better} (\code{"higher"} or
#'   \code{"lower"}) and \code{note}.
#' @references Truth and Reconciliation Commission of Canada (2015). Truth and
#' Reconciliation Commission of Canada: Calls to Action. Winnipeg. Government of
#' Canada catalogue IR4-8/2015E.
#' publications.gc.ca/collections/collection_2015/trc/IR4-8-2015-eng.pdf
#' @seealso \code{\link{trchlthgap}}, \code{\link{trcprog}}.
#' @examples
#' ti <- trcind()
#' ti[ti$call == "19", c("indicator", "table", "measure", "better")]
#' @export
trcind <- function() .igc("trcind")

#' The Call to Action 19 health gaps, in one place
#'
#' For each health indicator of Call to Action 19 that the published tables
#' measure, each group's estimate against non-Indigenous people's, with its
#' interval, turned into one scale, the disadvantage ratio: the factor by which
#' the group's outcome is worse (above
#' 1) or better (below 1).
#'
#' @details The gaps come from \code{\link{indgap}} on the latest period of each
#' table, with \code{ref} (non-Indigenous people) as the reference. To put
#' years, rates and percentages on one scale, each gap is given as a ratio
#' oriented by the indicator's better direction (\code{\link{trcind}}): the
#' group's estimate over the reference's where lower is better (infant
#' mortality, chronic conditions), the reference's over the group's where higher
#' is better (life expectancy, good mental health, a regular provider): \deqn{D
#' = \hat\theta_g / \hat\theta_r \ \textrm{(lower is better)}, \qquad D =
#' \hat\theta_r / \hat\theta_g \ \textrm{(higher is better)}.}{D = g / r (lower
#' is better), D = r / g (higher is better).} It is a ratio of the published
#' estimates, not a relative risk: for life expectancy, D = 1.1 means the
#' reference lives 10% longer, not that the group's risk is 10% higher. Above 1
#' the group fares worse; the interval is the published one carried through
#' (\code{indgap}'s log-scale interval, inverted where needed). The difference
#' on the indicator's own scale is kept beside it.
#'
#' This is a description of gaps that Call 19 asks to be closed, measured as the
#' surveys measure them; it is not a ranking of peoples. Report it with the
#' trends (\code{\link{healchg}}) and the strengths (\code{\link{welprof}}) that
#' sit beside it.
#'
#' Each table is in the long layout of \code{\link{scests}}: columns
#' \code{group}, \code{indicator} or \code{measure},
#' \code{period}, \code{estimate}, \code{lower}, \code{upper}, \code{quality},
#' and \code{age} for the health and suicidal-ideation tables (the all-ages
#' total is used). The latest period with estimates is taken. Table 41-10-0021
#' publishes Indigenous groups only, so it adds a gap only when given a table
#' with \code{ref} rows.
#'
#' @param health A health indicators table (table 13-10-0924), or \code{NULL}
#'   (the default).
#' @param lifeexp Life expectancy (table 17-10-0160), or \code{NULL} (the
#'   default).
#' @param infmort Infant mortality (table 17-10-0161), or \code{NULL} (the
#'   default).
#' @param suiide Suicidal thoughts (table 41-10-0021), or \code{NULL} (the
#'   default). The name is kept as it was first published;
#'   it is the suicidal-ideation table.
#' @param ref The reference group's label (default \code{"Non-Indigenous"}).
#' @param agetot The label of the all-ages row in the \code{age} column, one
#'   string for both tables or two (\code{health}, then \code{suiide}).
#'   Default \code{NULL}: Statistics Canada's
#'   \code{"Total, 18 years and over"} and
#'   \code{"Total, 15 years and over"}. A table of your own with no
#'   \code{age} column is used whole.
#' @return A data frame of class \code{trchlthgap}: \code{indicator},
#'   \code{measure}, \code{period}, \code{group}, \code{estimate},
#'   \code{refest}, \code{difference} (group minus reference, on the indicator's
#'   scale), \code{disratio}, \code{lower}, \code{upper}, \code{worse} (the
#'   interval above 1), \code{usewithcau}.
#' @references Truth and Reconciliation Commission of Canada (2015), see
#'   \code{\link{trcind}}.
#' @seealso \code{\link{trcind}}, \code{\link{indgap}}.
#' @examples
#' # your own table: made-up life expectancies with their limits; the tables
#' # you do not have are NULL
#' le <- data.frame(
#'   period = "2017/2019",
#'   group = c("First Nations", "Inuit", "Non-Indigenous"),
#'   estimate = c(74.2, 70.9, 81.8),
#'   lower = c(73.6, 69.1, 81.6),
#'   upper = c(74.8, 72.7, 82.0),
#'   quality = "ok"
#' )
#' trchlthgap(health = NULL, lifeexp = le, infmort = NULL, suiide = NULL)
#' @export
trchlthgap <- function(
  health = NULL,
  lifeexp = NULL,
  infmort = NULL,
  suiide = NULL,
  ref = "Non-Indigenous",
  agetot = NULL
) {
  .igc(
    "trchlthgap",
    health = health,
    lifeexp = lifeexp,
    infmort = infmort,
    suiide = suiide,
    ref = ref,
    agetot = agetot
  )
}

#' @export
format.trchlthgap <- function(x, ...) {
  need <- c(
    "indicator",
    "measure",
    "period",
    "group",
    "estimate",
    "refest",
    "disratio",
    "lower",
    "upper",
    "worse"
  )
  fb <- .igplain(x, need, "Truth and Reconciliation Call 19: health gaps")
  if (!is.null(fb)) {
    return(fb)
  }
  .igb(
    "Truth and Reconciliation Call 19: health gaps",
    facts = list(
      Measures = length(unique(x$measure)),
      Groups = paste(unique(x$group), collapse = "; "),
      `Clearly worse` = sprintf(
        "%d of %d gaps (interval above 1)",
        sum(.igyes(x$worse)),
        nrow(x)
      )
    ),
    tables = .thtabs(x),
    reading = c(
      "disratio is the disparity ratio, turned so that above 1 is worse for",
      "the group than for the reference (refest)."
    ),
    notes = c(
      "Call 19 asks for measurable goals to close the gaps in health",
      "outcomes; Call 55 for annual reports on progress.",
      .igcau(x$usewithcau)
    )
  )
}

#' A simulated policy panel for the examples
#'
#' Twenty districts in four regions over 36 months; a reform adopted region by
#' region in three stages (months 13, 19 and 25), with the fourth region never
#' adopting it; crime and arrest rates and counts. The numbers are made up, with
#' a known effect of the reform (a fall of 3 in the crime rate per 100,000), so
#' the examples and tests can check that each estimator recovers it.
#'
#' The panel is simulated: it describes no real place, district or reform,
#' and its numbers must not be read as data. The same seed gives the same
#' panel in every front end of the package (the draws are those of R's
#' default generator after \code{set.seed(seed)}), and the policy functions'
#' examples use it.
#'
#' @param seed The seed (default 20240101).
#' @return A data frame: \code{district}, \code{region}, \code{month},
#' \code{refmonth} (\code{NA} for never), \code{reformed} (0 or 1),
#' \code{population}, \code{logpop}, \code{crimes}, \code{crmrate},
#' \code{arrrate}.
#' @examples
#' # the panel is simulated: made-up numbers with a known effect of -3
#' p <- polex()
#' head(p)
#' # the same seed gives the same panel
#' identical(polex(7), polex(7))
#' @export
polex <- function(seed = 20240101) .igc("polex", seed = seed)

# format.cnsprof(): the distinct values of a column, or their count and
# range when there are more than four
.igvals <- function(v) {
  u <- unique(as.character(v[!is.na(v)]))
  if (length(u) <= 4L) {
    paste(u, collapse = "; ")
  } else {
    sprintf("%d, from %s to %s", length(u), u[1L], u[length(u)])
  }
}

# format.cnsprof(): what the profile holds and how much of it is flagged
.cpread <- function(x) {
  q <- as.character(x$quality)
  bad <- !is.na(q) & q != "ok"
  c(
    sprintf(
      "%d values of %d characteristics in %d geographies.",
      nrow(x),
      length(unique(x$chaid)),
      length(unique(x$geo))
    ),
    if (any(bad)) {
      tb <- table(q[bad])
      sprintf(
        "%d %s a quality flag (%s): read those with caution.",
        sum(bad),
        if (sum(bad) == 1L) "carries" else "carry",
        paste(names(tb), as.integer(tb), collapse = ", ")
      )
    } else {
      "None carries a quality flag."
    }
  )
}

# a row's label: its place, residence, age and group, totals left out
.iglab <- function(x) {
  parts <- lapply(c("geo", "residence", "age", "gender", "group"), function(k) {
    v <- as.character(x[[k]])
    ifelse(is.na(v) | grepl("^Total", v), NA, v)
  })
  apply(do.call(cbind, parts), 1L, function(r) {
    toString(r[!is.na(r)])
  })
}

# format.trchlthgap(): one table per indicator, measure and period
.thtabs <- function(x) {
  key <- sprintf("%s: %s (%s)", x$indicator, x$measure, x$period)
  tabs <- lapply(unique(key), function(k) {
    g <- x[key == k, , drop = FALSE]
    data.frame(
      group = .igstar(g$group, g$usewithcau),
      estimate = g$estimate,
      refest = g$refest,
      disratio = g$disratio,
      interval = .igci(g$lower, g$upper),
      worse = ifelse(.igyes(g$worse), "yes", "")
    )
  })
  names(tabs) <- unique(key)
  tabs
}
