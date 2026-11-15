#' Concentration-response functions for PM2.5 and NO2
#'
#' A concentration-response function (CRF) turns an air-pollution
#' concentration into the relative risk (RR) of a health outcome compared
#' with a counterfactual \dQuote{clean air} level: an RR of 1.05 means 5%
#' more deaths (or cases) than at the reference concentration. These are
#' the first step of any burden estimate; \code{\link{polbrdn}} calls them.
#'
#' They belong to environmental criminology in both its established senses:
#' the study of environmental crime and harm (green criminology), where the
#' burden a polluting activity places on an area or group, and how unevenly
#' it falls, is the measure of the harm; and the study of the environment
#' as a driver of offending and victimisation, where exposure series are
#' joined to crime counts by place and day, as in work relating air
#' pollution to crime (Bondy, Roth and Sager 2020; Herrnstadt et al.
#' 2019). The second example below is of that kind.
#'
#' @details
#' \code{crfpm25()} is log-linear for all-cause mortality with the pooled
#' cohort estimate of the WHO 2021 guideline review (Chen and Hoek 2020:
#' RR 1.08, 95% CI 1.06-1.09, per 10 micrograms per cubic metre),
#' \deqn{RR(z) = \exp(\ln(1.08) (z - z_{cf}) / 10)}
#' for \eqn{z > z_{cf}} and 1 otherwise; for the cause-specific outcomes
#' (ischaemic heart disease, stroke) it is the Integrated Exposure-Response
#' curve of Burnett et al. (2014, equation 1),
#' \deqn{RR(z) = 1 + \alpha (1 - e^{-\gamma (z - z_{cf})^{\delta}}),}
#' with the GBD 2013 parameter triples; the IER was fitted per cause and has
#' no all-cause form. \code{crfno2()} is the log-linear function of the WHO
#' 2021 review (Huangfu and Atkinson 2020: RR 1.02, 95% CI 1.01-1.04, per 10
#' micrograms per cubic metre for all-cause mortality),
#' \deqn{RR(z) = \exp(\beta (z - z_{cf}) / 10).}
#' A vector of exposures gives the mean RR and log-RR; the per-unit
#' relative risks are then in \code{extra$rrperunit}.
#'
#' Every function has zero excess risk at and below the reference
#' concentration (\eqn{z \le z_{cf}} gives RR 1). The published 95% limits
#' of the log-linear coefficients are carried through the same curve into
#' \code{rrlo} and \code{rrup}; the IER parameters come without an interval,
#' so those are \code{NA} for \code{"ihd"} and \code{"stroke"}, and also
#' when \code{betaper10} is supplied.
#'
#' Refusals: an unknown \code{y} (the message lists the available ones); an
#' exposure that is missing, infinite, negative, logical or empty; and a
#' \code{refconc} (or \code{betaper10}) that is not a single finite number,
#' or a negative \code{refconc}. A supplied \code{betaper10} may be negative
#' (a protective association) and replaces the table, in which case
#' \code{y} is only a label.
#'
#' @param expo Ambient concentration in micrograms per cubic metre, a single
#'   non-negative number or a vector of them (one per area, day or person).
#' @param y PM2.5: \code{"allcause"}, \code{"ihd"} or \code{"stroke"}
#'   (deaths). NO2: \code{"allcause"} and \code{"respiratory"} (deaths;
#'   Huangfu and Atkinson 2020, RR 1.02 and 1.03 per 10 micrograms per cubic
#'   metre) or \code{"chiast"} (incident childhood asthma, not deaths;
#'   Khreis et al. 2017, OR 1.05 per 4 micrograms per cubic metre, scaled to
#'   per 10).
#' @param refconc Counterfactual concentration below which no excess risk is
#'   assumed. PM2.5 default 5.8, the lower bound of the GBD 2019 theoretical
#'   minimum risk exposure level (the WHO 2021 guideline is 5); NO2 default
#'   10, the WHO 2021 guideline.
#' @param betaper10 Optional log relative risk per 10 micrograms per cubic
#'   metre (a single finite number, e.g. \code{log(1.04)}) overriding the NO2
#'   outcome table; \code{NULL} (the default) uses the table.
#' @param x A result.
#' @param ... Not used.
#' @return A list of class \code{crf} with: \code{rr} and \code{logrr}, the
#'   relative risk and its logarithm at the exposure (their means over a
#'   vector of exposures); \code{rrlo} and \code{rrup}, the RR at the
#'   published 95% limits (\code{NA} where there are none); \code{refconc};
#'   \code{expconc} (the exposure, or the mean of a vector); \code{pollutant}
#'   (\code{"PM2.5"} or \code{"NO2"}); \code{citation} (the source as a
#'   string); and \code{extra}, a list with the coefficient
#'   (\code{betaper10} and \code{betaper10ci} for the log-linear forms,
#'   \code{alpha}, \code{gamma}, \code{delta} for the IER), \code{form}
#'   (PM2.5 only), \code{outcome}, and \code{rrperunit}, the RR at each
#'   exposure (\code{NULL} for one exposure).
#' @references Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the
#'   air: the contemporaneous relationship between air pollution and crime.
#'   Journal of the Association of Environmental and Resource Economists
#'   7(3), 555-585. \doi{10.1086/707127}. Herrnstadt, E., Heyes, A.,
#'   Muehlegger, E. and Saberian, S. (2019). Air pollution as a cause of
#'   violent crime. Working paper. Chen, J. and Hoek, G. (2020). Long-term
#'   exposure to PM and
#'   all-cause and cause-specific mortality: a systematic review and
#'   meta-analysis. Environment International 143, 105974. Huangfu, P. and
#'   Atkinson, R. (2020). Long-term exposure to NO2 and O3 and all-cause and
#'   respiratory mortality: a systematic review and meta-analysis.
#'   Environment International 144, 105998. Khreis, H. et al. (2017).
#'   Exposure to traffic-related air pollution and risk of development of
#'   childhood asthma: a systematic review and meta-analysis. Environment
#'   International 100, 1-31. GBD 2019 Risk Factors Collaborators (2020). The
#'   Lancet 396(10258), 1223-1249. Burnett, R. T. et al. (2014). An
#'   integrated risk function for estimating the global burden of disease
#'   attributable to ambient fine particulate matter exposure. Environmental
#'   Health Perspectives 122(4), 397-403. WHO (2021). Global Air Quality
#'   Guidelines.
#' @examples
#' # All-cause mortality at an annual PM2.5 of 12 micrograms per cubic metre
#' r <- crfpm25(12)
#' c(rr = r$rr, lower = r$rrlo, upper = r$rrup)
#' # the cause-specific IER curve for ischaemic heart disease
#' crfpm25(12, y = "ihd")$rr
#' # NO2 and respiratory deaths; at or below the reference no excess risk
#' crfno2(25, y = "respiratory")$rr
#' crfno2(8)$rr
#' # a vector of exposures (one per area): the mean RR and each area's RR
#' crfpm25(c(5, 12, 20))$extra$rrperunit
#' # a coefficient of the reader's own: RR 1.04 per 10 micrograms
#' crfno2(30, betaper10 = log(1.04))$rr
#' # the environment as a driver of offending: a simulated daily NO2 series
#' # and simulated daily counts of assaults for three areas, joined by area
#' # and day
#' set.seed(1)
#' d <- data.frame(area = rep(c("a", "b", "c"), each = 30), day = rep(1:30, 3))
#' d$no2 <- stats::rlnorm(90, log(20), 0.3)
#' d$assaults <- stats::rpois(90, 3 + 0.05 * d$no2)
#' # the relative risk each area-day's exposure carries, and the mean count
#' # on days above and below its median
#' d$rr <- crfno2(d$no2)$extra$rrperunit
#' tapply(d$assaults, d$rr > stats::median(d$rr), mean)
#' # an outcome without a published function is refused
#' try(crfpm25(12, y = "asthma"))
#' @export
crfpm25 <- function(expo, y = "allcause", refconc = 5.8) {
  f <- "crfpm25"
  .crfout(y, c("allcause", "ihd", "stroke"), f, "")
  refconc <- .polsc(refconc, "refconc", f)
  z <- .polexp(expo, f)
  if (y == "allcause") {
    b <- log(1.08)
    ci <- log(c(1.06, 1.09))
    return(.crfobj(
      z, refconc, 0L, b, ci, "PM2.5",
      paste0(
        "Chen & Hoek (2020) Environ Int 143:105974; WHO (2021) Global AQ ",
        "Guidelines"
      ),
      list(betaper10 = b, betaper10ci = ci, form = "log-linear", outcome = y)
    ))
  }
  p <- list(ihd = c(1.91, 0.14, 0.49), stroke = c(1.46, 0.13, 0.61))[[y]]
  .crfobj(
    z, refconc, 1L, p, NULL, "PM2.5",
    "Burnett et al. (2014) EHP 122(4):397-403",
    list(alpha = p[1L], gamma = p[2L], delta = p[3L], form = "IER", outcome = y)
  )
}

#' @rdname crfpm25
#' @export
crfno2 <- function(expo, y = "allcause", refconc = 10, betaper10 = NULL) {
  f <- "crfno2"
  ci <- c(NA_real_, NA_real_)
  if (is.null(betaper10)) {
    .crfout(
      y, c("allcause", "respiratory", "chiast"), f,
      " or pass betaper10 explicitly"
    )
    b <- list(
      allcause = log(1.02), respiratory = log(1.03),
      chiast = log(1.05) * 10 / 4
    )[[y]]
    ci <- list(
      allcause = log(c(1.01, 1.04)), respiratory = log(c(1.00, 1.05)),
      chiast = log(c(1.02, 1.07)) * 10 / 4
    )[[y]]
  } else {
    b <- .polsc(betaper10, "betaper10", f, nonneg = FALSE)
  }
  refconc <- .polsc(refconc, "refconc", f)
  z <- .polexp(expo, f)
  cite <- if (identical(y, "chiast") && is.null(betaper10)) {
    "Khreis et al. (2017) Environ Int 100:1-31; WHO (2021) Global AQ Guidelines"
  } else {
    paste0(
      "Huangfu & Atkinson (2020) Environ Int 144:105998; WHO (2021) ",
      "Global AQ Guidelines"
    )
  }
  .crfobj(
    z, refconc, 0L, b, ci, "NO2", cite,
    list(betaper10 = b, betaper10ci = ci, outcome = y)
  )
}

.crfout <- function(y, avail, fn, tail) {
  if (!is.character(y) || length(y) != 1L || !y %in% avail) {
    .ckstop(
      fn, "Unknown outcome '%s'. Available: %s%s", paste(y, collapse = ","),
      toString(avail), tail
    )
  }
}

# The curve at the exposures, its limits through the published interval.
.crfobj <- function(z, refconc, form, par, ci, pol, cite, extra) {
  r <- .Call(C_clab_r_crf, z, refconc, form, par)
  lo <- up <- rep(NA_real_, length(z))
  if (!is.null(ci)) {
    lo <- .Call(C_clab_r_crf, z, refconc, 0L, ci[1L])[[1L]]
    up <- .Call(C_clab_r_crf, z, refconc, 0L, ci[2L])[[1L]]
  }
  one <- length(z) == 1L
  m <- function(v) if (one) v[[1L]] else mean(v)
  structure(
    list(
      rr = m(r[[1L]]),
      logrr = m(r[[2L]]),
      rrlo = m(lo),
      rrup = m(up),
      refconc = refconc,
      expconc = m(z),
      pollutant = pol,
      citation = cite,
      extra = c(extra, list(rrperunit = if (one) NULL else r[[1L]]))
    ),
    class = "crf"
  )
}

# One finite number, optionally non-negative or positive.
.polsc <- function(x, what, fn, nonneg = TRUE, positive = FALSE) {
  if (is.logical(x) && !all(is.na(x))) {
    .ckstop(fn, "`%s` must be numeric, not logical", what)
  }
  if (is.factor(x)) {
    .ckstop(
      fn,
      "`%s` must be numeric, not a factor (its level codes are not values)",
      what
    )
  }
  if (is.null(x)) .ckstop(fn, "`%s` must be numeric, got NULL", what)
  if (is.list(x) || !(is.numeric(x) || is.logical(x))) {
    .ckstop(fn, "`%s` must be numeric, got %s", what, class(x)[1L])
  }
  v <- as.numeric(x)
  if (length(v) != 1L || is.na(v) || !is.finite(v)) {
    .ckstop(fn, "`%s` must be a single finite number", what)
  }
  if (positive && v <= 0) .ckstop(fn, "`%s` must be positive", what)
  if (nonneg && v < 0) .ckstop(fn, "`%s` must be non-negative", what)
  v
}

.polexp <- function(expo, fn) {
  if (is.logical(expo)) .ckstop(fn, "`expo` must be numeric, not logical")
  if (is.null(expo)) .ckstop(fn, "`expo` must be numeric, got NULL")
  if (is.factor(expo) || is.list(expo) || !is.numeric(expo)) {
    .ckstop(fn, "`expo` must be numeric, got %s", class(expo)[1L])
  }
  z <- as.numeric(expo)
  if (!length(z) || anyNA(z) || !all(is.finite(z))) {
    .ckstop(
      fn,
      paste(
        "`expo` must be one or more finite concentrations, with no missing",
        "values"
      )
    )
  }
  if (any(z < 0)) .ckstop(fn, "`expo` must be non-negative")
  z
}

#' Population attributable fraction
#'
#' The share of all cases in a population that is due to an exposure: the
#' fraction that would not occur if the exposure were removed and nothing
#' else changed. It combines how strongly the exposure raises risk (the
#' relative risk) with how common it is (the prevalence); multiplied by the
#' number of cases it gives the attributable cases, as \code{\link{polbrdn}}
#' does.
#'
#' @details
#' Levin's formula (Levin 1953; Rothman, Greenland and Lash 2008, chapter
#' 5): \deqn{PAF = p (RR - 1) / (1 + p (RR - 1)),} with \eqn{p} the
#' proportion exposed. With everyone exposed it is \eqn{1 - 1/RR}. It
#' assumes the relative risk is causal and unconfounded and that removing
#' the exposure changes nothing else. A relative risk below 1 gives a
#' negative fraction (the exposure protects).
#'
#' When both confidence limits of the RR are given, the fraction's limits
#' are the formula applied to them: the formula increases with the RR, so
#' the interval of the fraction is the image of the interval of the risk.
#' Refusals: an \code{rr} or limit that is not a single finite positive
#' number (logical values included), a prevalence outside \eqn{[0, 1]},
#' only one of the two limits, or limits that do not enclose \code{rr}.
#'
#' @param rr Relative risk among the exposed, a single finite positive
#'   number.
#' @param expoprev Proportion of the population exposed, in \code{[0, 1]}.
#' @param rrlo,rrup Optional confidence limits of \code{rr}, both or
#'   neither.
#' @param x An \code{attrfrac} result.
#' @param ... Not used.
#' @return A single number of class \code{attrfrac}, the attributable
#'   fraction: in \code{[0, 1)} for \code{rr >= 1}, negative for a
#'   protective exposure; when the limits were given it carries an attribute
#'   \code{ci}, a named numeric vector (\code{lower}, \code{upper}), and it
#'   keeps its inputs as attributes \code{rr} and \code{expoprev}.
#'   \code{as.numeric()} gives the plain number.
#' @seealso \code{\link{crfpm25}} for a relative risk from a concentration,
#'   \code{\link{polbrdn}} for the attributable cases.
#' @references Rothman, K. J., Greenland, S. and Lash, T. L. (2008). Modern
#'   Epidemiology, 3rd ed., chapter 5. Levin, M. L. (1953). The occurrence of
#'   lung cancer in man. Acta Unio Internationalis Contra Cancrum 9,
#'   531-541.
#' @examples
#' attrfrac(1.5, 0.4)
#' attrfrac(2, 1) # everyone exposed: 1 - 1/RR
#' paf <- attrfrac(1.08, 1, rrlo = 1.06, rrup = 1.09)
#' attr(paf, "ci")
#' attrfrac(0.8, 0.5) # protective
#' try(attrfrac(1.5, 40))
#' @export
attrfrac <- function(rr, expoprev, rrlo = NULL, rrup = NULL) {
  f <- "attrfrac"
  rr <- .polsc(rr, "rr", f, nonneg = FALSE, positive = TRUE)
  p <- .polsc(expoprev, "expoprev", f)
  if (p > 1) .ckstop(f, "expoprev must be in [0,1]; got %s", format(p))
  paf <- .Call(C_clab_r_levin, rr, p)
  ci <- NULL
  if (!is.null(rrlo) || !is.null(rrup)) {
    if (is.null(rrlo) || is.null(rrup)) .ckstop(f, "rrlo and rrup go together")
    lo <- .polsc(rrlo, "rrlo", f, nonneg = FALSE, positive = TRUE)
    hi <- .polsc(rrup, "rrup", f, nonneg = FALSE, positive = TRUE)
    if (lo > rr || hi < rr) .ckstop(f, "rr must lie between rrlo and rrup")
    ci <- stats::setNames(
      .Call(C_clab_r_levin, c(lo, hi), p), c("lower", "upper")
    )
  }
  structure(paf, ci = ci, rr = rr, expoprev = p, class = "attrfrac")
}

#' @rdname attrfrac
#' @export
print.attrfrac <- function(x, ...) {
  ci <- attr(x, "ci")
  cat(sprintf(
    "Attributable fraction %.4g (relative risk %s, %s exposed)%s\n",
    as.numeric(x), format(attr(x, "rr")), format(attr(x, "expoprev")),
    if (is.null(ci)) "" else sprintf(", 95%% %.4g to %.4g", ci[[1L]], ci[[2L]])
  ))
  if (as.numeric(x) < 0) cat("Negative: the exposure protects.\n")
  invisible(x)
}

#' Deaths avoided under a counterfactual exposure reduction
#'
#' The expected number of deaths (or cases) a year that a given drop in
#' exposure would avoid in a population, from the baseline rate and a
#' log-linear concentration-response coefficient. This is the health-impact
#' function used in regulatory benefit assessments.
#'
#' @details
#' The BenMAP-CE health-impact function (US EPA 2018; Anenberg et al.
#' 2010): \deqn{\Delta Y = y_0 N (1 - e^{-\beta \Delta x}),} with \eqn{y_0}
#' the baseline rate per person-year, \eqn{N} the population, \eqn{\beta}
#' the log relative risk per unit of exposure and \eqn{\Delta x} the
#' reduction. A negative \code{expdel} (an increase) gives a negative
#' number, the extra cases. A coefficient published per 10 units must be
#' divided by 10 first (for example \code{log(1.08) / 10} for the PM2.5
#' all-cause RR of 1.08 per 10).
#'
#' The baseline rate is per person-year, not per 100,000: 800 deaths per
#' 100,000 a year is 0.008. A rate below 1e-5 or above 0.1 is refused unless
#' \code{rare} or \code{hirate} says it is meant, because it is also what a
#' rate in the wrong unit looks like; a rate above 1, or below 1e-11, is
#' always refused. Every input must be a single finite number, and the
#' population and rate non-negative.
#'
#' @param expdel Counterfactual reduction in exposure, in the unit of
#'   \code{betaperunit} (positive means the exposure goes down).
#' @param pop At-risk population, a non-negative number of persons.
#' @param brate Baseline rate in cases per person-year (0.008 is 800 per
#'   100,000).
#' @param betaperunit Log relative risk per unit of exposure (per unit, not
#'   per 10), a single finite number.
#' @param rare,hirate \code{TRUE} to admit a baseline rate below 1e-5 or
#'   above 0.1 per person-year.
#' @return A single number: the expected cases avoided per year (negative
#'   for an increase in exposure).
#' @seealso \code{\link{polbrdn}}.
#' @references US EPA (2018). Environmental Benefits Mapping and Analysis
#'   Program - Community Edition, user manual. Anenberg, S. C. et al.
#'   (2010). An estimate of the global burden of anthropogenic ozone and fine
#'   particulate matter on premature human mortality using atmospheric
#'   modeling. Environmental Health Perspectives 118(9), 1189-1195.
#' @examples
#' # a 10-unit cut, a million people, 800 deaths per 100,000 a year
#' mortdisp(10, 1e6, 0.008, 0.0039)
#' # the PM2.5 all-cause coefficient is published per 10 micrograms
#' mortdisp(5, 2.8e6, 0.0075, log(1.08) / 10)
#' # a rare outcome must be declared
#' mortdisp(5, 1e6, 5e-6, 0.01, rare = TRUE)
#' # a rate per 100,000 passed as a rate per person is refused
#' try(mortdisp(10, 1e6, 800, 0.0039))
#' @export
mortdisp <- function(
  expdel,
  pop,
  brate,
  betaperunit,
  rare = FALSE,
  hirate = FALSE
) {
  f <- "mortdisp"
  expdel <- .polsc(expdel, "expdel", f, nonneg = FALSE)
  pop <- .polsc(pop, "pop", f)
  brate <- .polrate(.polsc(brate, "brate", f), rare, hirate, f)
  betaperunit <- .polsc(betaperunit, "betaperunit", f, nonneg = FALSE)
  .Call(C_clab_r_hif, expdel, pop, brate, betaperunit)
}

# A baseline rate per person-year, checked against the band.
.polrate <- function(rate, rare, hirate, fn) {
  for (a in list(list("rare", rare), list("hirate", hirate))) {
    if (!is.logical(a[[2L]]) || length(a[[2L]]) != 1L || is.na(a[[2L]])) {
      .ckstop(fn, "`%s` must be TRUE or FALSE", a[[1L]])
    }
  }
  if (rate > 1) {
    .ckstop(
      fn,
      paste(
        "`brate` is a rate per person-year, between 0 and 1 (0.008 is 800",
        "per 100,000); %s looks like a rate per 100,000"
      ),
      format(rate)
    )
  }
  if (rate == 0) return(rate)
  if (rate < 1e-11) {
    .ckstop(
      fn,
      paste(
        "`brate` %s per person-year is below one case per 100 billion",
        "person-years: not a usable rate"
      ),
      format(rate)
    )
  }
  if (rate < 1e-5 && !rare) {
    .ckstop(
      fn,
      paste(
        "`brate` %s per person-year is below one case per 100,000 a year. A",
        "rate per 100,000 divided by 100,000 is the rate here (800 per",
        "100,000 is 0.008); if the outcome really is this rare, say so with",
        "rare = TRUE"
      ),
      format(rate)
    )
  }
  if (rate > 0.1 && !hirate) {
    .ckstop(
      fn,
      paste(
        "`brate` %s per person-year is %s per 100,000 a year, above a whole",
        "population's all-cause death rate (an age-specific rate for the",
        "oldest ages can be this high). If it is a rate per 100,000, divide",
        "it by 100,000; if the outcome really is this common, say so with",
        "hirate = TRUE"
      ),
      format(rate),
      format(rate * 1e5, big.mark = ",", scientific = FALSE)
    )
  }
  rate
}

#' End-to-end pollution burden
#'
#' Chains the concentration-response function, the attributable fraction
#' and the baseline case count into annual attributable cases, as the GBD
#' risk factor estimates do (GBD 2019 Risk Factors Collaborators 2020). It
#' serves a single population (a city, a region) whose mean exposure,
#' baseline rate and size are known; for many areas at once see
#' \code{\link{polbrdarea}}.
#'
#' @details
#' The steps are: the relative risk at \code{expomean} from
#' \code{\link{crfpm25}} or \code{\link{crfno2}} (with the reference
#' concentration given or the pollutant's default); the attributable
#' fraction from \code{\link{attrfrac}} with \code{expoprev}; the baseline
#' cases \eqn{y_0 N}; and the attributable cases, their product. The result
#' is a yearly count of deaths (or incident cases for \code{"chiast"}).
#'
#' Refusals: an unknown pollutant, an outcome the pollutant's CRF does not
#' have, any input that is not a single finite non-negative number, and a
#' baseline rate outside the band described in \code{\link{mortdisp}}
#' unless \code{rare} or \code{hirate} admits it.
#'
#' @param expomean The mean exposure among the exposed, in micrograms per
#'   cubic metre. With \code{expoprev = 1} that is the population mean.
#' @param expoprev Proportion of the population at that level (1 for
#'   ambient air), in \code{[0, 1]}.
#' @param brate Cases per person-year in the unexposed scenario (0.008 is
#'   800 per 100,000).
#' @param pop At-risk population (persons).
#' @param poll \code{"PM2.5"} (also written \code{"PM25"}) or \code{"NO2"},
#'   in either case.
#' @param y Outcome key passed to the concentration-response function.
#' @param refconc Counterfactual concentration; \code{NULL} takes the
#'   pollutant's default (PM2.5 5.8, NO2 10).
#' @param rare,hirate As in \code{\link{mortdisp}}.
#' @param x A result.
#' @param ... Not used.
#' @return A list of class \code{polbrdn} with: \code{paf},
#'   \code{attrcases} and \code{basecases} (per year), \code{population},
#'   \code{brate} (per person-year), \code{expomean}, \code{refconc} (the one
#'   used), \code{pollutant}, \code{citation}, and \code{extra}, a list with
#'   \code{rr}, \code{logrr}, \code{outcome} and \code{unit}
#'   (\code{"deaths"} or \code{"incident cases"}).
#' @references GBD 2019 Risk Factors Collaborators (2020). Global burden of
#'   87 risk factors in 204 countries and territories, 1990-2019. The Lancet
#'   396(10258), 1223-1249.
#' @examples
#' # ambient NO2 of 25 micrograms per cubic metre over a million people
#' b <- polbrdn(25, 1, 0.008, 1e6, poll = "NO2")
#' b$attrcases
#' # PM2.5 and ischaemic heart disease deaths, against the WHO guideline
#' polbrdn(12, 1, 0.0012, 2.8e6, poll = "PM2.5", y = "ihd",
#'         refconc = 5)$attrcases
#' # childhood asthma counts incident cases, not deaths
#' polbrdn(30, 1, 0.012, 5e5, y = "chiast")$extra$unit
#' try(polbrdn(40, 1, 0.008, 1e6, poll = "O3"))
#' @export
polbrdn <- function(
  expomean,
  expoprev,
  brate,
  pop,
  poll = "NO2",
  y = "allcause",
  refconc = NULL,
  rare = FALSE,
  hirate = FALSE
) {
  f <- "polbrdn"
  pol <- .polkey(poll)
  if (is.na(pol)) {
    .ckstop(
      f, "Unknown pollutant '%s'. Use 'PM2.5' or 'NO2'",
      paste(poll, collapse = ",")
    )
  }
  expomean <- .polsc(expomean, "expomean", f)
  expoprev <- .polsc(expoprev, "expoprev", f)
  brate <- .polrate(.polsc(brate, "brate", f), rare, hirate, f)
  pop <- .polsc(pop, "pop", f)
  if (!is.null(refconc)) refconc <- .polsc(refconc, "refconc", f)
  crf <- if (pol == "PM2.5") {
    crfpm25(expomean, y = y, refconc = if (is.null(refconc)) 5.8 else refconc)
  } else {
    crfno2(expomean, y = y, refconc = if (is.null(refconc)) 10 else refconc)
  }
  if (expoprev > 1) {
    .ckstop(f, "expoprev must be in [0,1]; got %s", format(expoprev))
  }
  paf <- .Call(C_clab_r_levin, crf$rr, expoprev)
  base <- brate * pop
  structure(
    list(
      paf = paf,
      attrcases = paf * base,
      basecases = base,
      population = pop,
      brate = brate,
      expomean = expomean,
      refconc = crf$refconc,
      pollutant = pol,
      citation = paste0(crf$citation, "; Rothman et al. (2008) section 5"),
      extra = list(
        rr = crf$rr,
        logrr = crf$logrr,
        outcome = y,
        unit = if (identical(y, "chiast")) "incident cases" else "deaths"
      )
    ),
    class = "polbrdn"
  )
}

.polkey <- function(poll) {
  if (!is.character(poll) || length(poll) != 1L || is.na(poll)) {
    return(NA_character_)
  }
  k <- gsub("[^a-z0-9]", "", tolower(poll))
  if (k == "pm25") return("PM2.5")
  if (k == "no2") return("NO2")
  NA_character_
}

.colnm <- function(x, what, fn, several = FALSE) {
  ok <- is.character(x) && length(x) >= 1L && !anyNA(x) && all(nzchar(x)) &&
    (several || length(x) == 1L)
  if (!ok) {
    .ckstop(
      fn, "`%s` must name %s of `data`", what,
      if (several) "one or more columns" else "one column"
    )
  }
  x
}

#' Per-area pollution burden
#'
#' Applies \code{\link{polbrdn}} to every row of a table with one row per
#' area (forward sortation area, census tract, neighbourhood, ...) and
#' returns the rows sorted by attributable cases, worst first: where in a
#' city or region the pollution burden falls.
#'
#' @details
#' Each area is treated as wholly exposed at its own concentration
#' (\code{expoprev = 1}), with the pollutant's default reference
#' concentration. Every area is checked as \code{\link{polbrdn}} checks a
#' single population, and the first area that fails stops the call with
#' that area's error. Refusals: a column name that is not a single string or
#' not in the table, an empty table, and a non-numeric exposure, population
#' or rate column. The area labels can be of any type.
#'
#' @param areatab A data frame with one row per area: a label, the exposure
#'   in micrograms per cubic metre, the population in persons and the
#'   baseline rate per person-year.
#' @param areacol,expcol,popcol,basratecol Column names (single strings) of
#'   the label, exposure, population and rate; defaults \code{"fsa"},
#'   \code{"exposure"}, \code{"population"} and \code{"brate"}.
#' @param poll,y,rare,hirate Passed to \code{\link{polbrdn}} for every area.
#' @param x A result.
#' @param ... Not used.
#' @return A data frame of class \code{polbrdarea} with one row per area,
#'   sorted by \code{attrcases}, largest first: the four input columns plus
#'   \code{rr}, \code{paf}, \code{attrcases} and \code{basecases}.
#' @seealso \code{\link{expconcidx}} for how the burden is spread across
#'   income.
#' @examples
#' # simulated areas
#' t <- data.frame(fsa = c("A1", "B2"), exposure = c(28, 18),
#'                 population = c(40000, 60000), brate = 0.008)
#' polbrdarea(t)
#' try(polbrdarea(t[, c("fsa", "exposure")]))
#' @export
polbrdarea <- function(
  areatab,
  areacol = "fsa",
  expcol = "exposure",
  popcol = "population",
  basratecol = "brate",
  poll = "NO2",
  y = "allcause",
  rare = FALSE,
  hirate = FALSE
) {
  f <- "polbrdarea"
  for (a in c("areacol", "expcol", "popcol", "basratecol")) {
    .colnm(get(a), a, f)
  }
  need <- c(areacol, expcol, popcol, basratecol)
  miss <- setdiff(need, names(areatab))
  if (length(miss)) .ckstop(f, "areatab missing columns: %s", toString(miss))
  if (!nrow(areatab)) .ckstop(f, "areatab has no rows")
  for (cl in c(expcol, popcol, basratecol)) {
    if (!is.numeric(areatab[[cl]])) {
      .ckstop(
        f,
        paste0(
          "column `%s` must be numeric (a factor or character column is ",
          "not a count)"
        ),
        cl
      )
    }
  }
  rows <- lapply(seq_len(nrow(areatab)), function(i) {
    r <- tryCatch(
      polbrdn(
        expomean = as.numeric(areatab[[expcol]][i]),
        expoprev = 1,
        brate = as.numeric(areatab[[basratecol]][i]),
        pop = as.numeric(areatab[[popcol]][i]),
        poll = poll, y = y, rare = rare, hirate = hirate
      ),
      error = function(e) {
        m <- sub("^(polbrdn|crfpm25|crfno2): ", "", conditionMessage(e))
        m <- sub("\\.$", "", m)
        mp <- c(expomean = "expcol", expo = "expcol", pop = "popcol",
                brate = "basratecol")
        for (k in names(mp)) {
          m <- gsub(
            sprintf("`%s`", k), sprintf("`%s`", mp[[k]]), m,
            fixed = TRUE
          )
        }
        .ckstop(f, "%s", m)
      }
    )
    out <- data.frame(a = areatab[[areacol]][i], stringsAsFactors = FALSE)
    names(out) <- areacol
    out[[expcol]] <- as.numeric(areatab[[expcol]][i])
    out[[popcol]] <- as.numeric(areatab[[popcol]][i])
    out[[basratecol]] <- as.numeric(areatab[[basratecol]][i])
    out$rr <- r$extra$rr
    out$paf <- r$paf
    out$attrcases <- r$attrcases
    out$basecases <- r$basecases
    out
  })
  res <- do.call(rbind, rows)
  res <- res[order(-res$attrcases), , drop = FALSE]
  rownames(res) <- NULL
  class(res) <- c("polbrdarea", "data.frame")
  res
}

#' Concentration index of exposure by income
#'
#' Measures whether exposure to pollution is concentrated among lower-income
#' or higher-income people. The index runs from -1 (all the exposure borne
#' by the poorest) through 0 (exposure shared evenly across income) to +1
#' (all of it borne by the richest); environmental-justice analyses report
#' it next to the burden.
#'
#' @details
#' Wagstaff, Paci and van Doorslaer (1991):
#' \deqn{CI = (2 / \mu) \, \mathrm{cov}(h_i, R_i)}
#' with \eqn{R_i = (rank_i - 0.5) / n} the fractional income rank (tied
#' incomes share their average rank, so the index does not depend on row
#' order) and the population covariance. Negative values mean lower-income
#' units bear more exposure. Both columns must be numeric: a factor or
#' character income is refused rather than ranked by its level codes, and
#' exposures must be non-negative.
#'
#' Rows with a missing exposure or income are dropped. The interpretation
#' calls an index below -0.01 a burden on the poor, above 0.01 a burden on
#' the rich, and anything between an even distribution. Income only needs
#' an order: an amount, a quintile number or a rank all give the same
#' index.
#'
#' @param data A data frame with one row per person or area.
#' @param expo,income Column names (single strings): the exposure (a
#'   non-negative concentration) and the income measure (numeric).
#' @param x A result.
#' @param ... Not used.
#' @return A list of class \code{expconcidx} with \code{concidx} (a number
#'   in \eqn{[-1, 1]}), \code{interpretation} (a sentence), \code{expomean},
#'   \code{citation} and \code{extra}, a list with \code{n} (rows used),
#'   \code{cov_h_R} (the covariance) and \code{curve}, a data frame of the
#'   concentration curve: \code{population}, the cumulative population
#'   share, poorest first, and \code{exposure}, the cumulative share of total
#'   exposure.
#' @references Wagstaff, A., Paci, P. and van Doorslaer, E. (1991). On the
#'   measurement of inequalities in health. Social Science and Medicine
#'   33(5), 545-557.
#' @examples
#' d <- data.frame(exposure = c(30, 25, 20, 15, 10), income = 1:5)
#' ci <- expconcidx(d, "exposure", "income")
#' ci$concidx
#' ci$interpretation
#' try(expconcidx(data.frame(exposure = 1:3, income = c("a", "b", "c")),
#'                "exposure", "income"))
#' @export
expconcidx <- function(data, expo, income) {
  f <- "expconcidx"
  .colnm(expo, "expo", f)
  .colnm(income, "income", f)
  miss <- setdiff(c(expo, income), names(data))
  if (length(miss)) .ckstop(f, "data is missing columns: %s", toString(miss))
  df <- data[, c(expo, income), drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 2L) {
    .ckstop(f, "Need at least 2 rows to compute concentration index")
  }
  h <- .cknum(df[[expo]], expo, f, scalar = FALSE)
  inc <- .cknum(df[[income]], income, f, scalar = FALSE)
  if (any(h < 0)) {
    .ckstop(
      f, "`%s` must be non-negative: an exposure is a concentration", expo
    )
  }
  if (length(unique(inc)) < 2L) {
    .ckstop(
      f,
      paste0(
        "`%s` must take at least two values: with one income level there ",
        "is no income gradient and the index is undefined"
      ),
      income
    )
  }
  if (mean(h) <= 0) {
    .ckstop(f, "Mean exposure is zero; concentration index undefined")
  }
  r <- .Call(C_clab_r_concidx, h, inc)
  ci <- r[[1L]][1L]
  structure(
    list(
      concidx = ci,
      interpretation = .concread(ci),
      expomean = r[[1L]][3L],
      citation = paste0(
        "Wagstaff, Paci & van Doorslaer ",
        "(1991) Soc Sci Med 33(5):545-557"
      ),
      extra = list(
        n = nrow(df), cov_h_R = r[[1L]][2L],
        curve = data.frame(population = r[[2L]], exposure = r[[3L]])
      )
    ),
    class = "expconcidx"
  )
}

.concread <- function(ci) {
  if (ci < -0.01) {
    sprintf(
      paste0(
        "CI = %.4f. Pro-poor exposure burden: lower-income individuals ",
        "bear disproportionately higher pollution."
      ),
      ci
    )
  } else if (ci > 0.01) {
    sprintf(
      paste0(
        "CI = %.4f. Pro-rich exposure burden: higher-income individuals ",
        "bear disproportionately higher pollution."
      ),
      ci
    )
  } else {
    sprintf(
      "CI = %.4f. Exposure distributed approximately evenly across income.",
      ci
    )
  }
}

#' @rdname crfpm25
#' @export
print.crf <- function(x, ...) {
  cat(sprintf(
    "Concentration-response, %s, %s: RR %.6g at %s ug/m3 (reference %s)%s\n",
    x$pollutant, x$extra$outcome, x$rr, format(x$expconc), format(x$refconc),
    if (is.finite(x$rrlo) && is.finite(x$rrup)) {
      sprintf(", 95%% %.6g to %.6g", x$rrlo, x$rrup)
    } else {
      ""
    }
  ))
  if (length(x$extra$rrperunit) > 1L) {
    cat("Several exposures were given: the risk shown is their mean.\n")
  }
  cat("Source:", x$citation, "\n")
  invisible(x)
}

#' @rdname polbrdn
#' @export
print.polbrdn <- function(x, ...) {
  cat(sprintf(
    paste0(
      "Burden of %s, %s: about %.4g of the %.4g %s a year are attributable",
      " to %s above %s ug/m3 (attributable fraction %.2f%%, RR %.6g)\n"
    ),
    x$pollutant, x$extra$outcome, x$attrcases, x$basecases, x$extra$unit,
    x$pollutant, format(x$refconc), 100 * x$paf, x$extra$rr
  ))
  cat("Source:", x$citation, "\n")
  invisible(x)
}

#' @rdname polbrdarea
#' @export
print.polbrdarea <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  att <- sum(d$attrcases, na.rm = TRUE)
  base <- sum(d$basecases, na.rm = TRUE)
  cat(sprintf(
    "Pollution burden by area: %d areas, %.4g attributable of %.4g cases%s\n",
    nrow(d), att, base,
    if (base > 0) sprintf(" (%.3g%%)", 100 * att / base) else ""
  ))
  print(utils::head(d, 40L), row.names = FALSE, digits = 4L)
  invisible(x)
}

#' @rdname expconcidx
#' @export
print.expconcidx <- function(x, ...) {
  cat(sprintf(
    "Concentration index of exposure by income: %.4f (%d rows)\n",
    x$concidx, x$extra$n
  ))
  cat(strwrap(x$interpretation, width = 70L), sep = "\n")
  invisible(x)
}
