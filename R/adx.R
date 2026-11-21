# Alcohol and drug use: a blood alcohol concentration from what was drunk,
# standard drinks, and the scores of three screening questionnaires. The
# command line and the other languages compute the same in C++
# (cli/jval.inc); the cases of cli/test.tsv hold the two together.

#' Estimated blood alcohol concentration (Widmark)
#'
#' The blood alcohol concentration a person would have after a number of
#' standard drinks, by Widmark's formula: the alcohol drunk, spread over
#' the share of the body that takes it up, less what has been eliminated
#' since drinking began. It is the estimate surveys of drinking report
#' (an "eBAC") and toxicologists use for back-calculation; it is not a
#' measurement, and a person's own rate of elimination and body water
#' move it by a third either way.
#'
#' \deqn{BAC = \frac{g \, d}{10 \, r \, W} - \beta t}{
#'   BAC = g d / (10 r W) - beta t}
#' in grams per 100 mL, with `d` drinks of `g` grams of ethanol each, `W`
#' the body weight in kilograms, `r` the Widmark factor, `beta` the
#' elimination per hour and `t` the hours since drinking began; a
#' negative value is 0. The form with fluid ounces and pounds (Matthews
#' and Miller 1979), `5.14 oz / (r lb) - 0.015 t`, is the same formula:
#' 5.14 is 100 times the grams in a fluid ounce of ethanol over the grams
#' in a pound.
#'
#' The default `r`, 0.73 for men and 0.66 for women, are the factors of
#' that form; Widmark's own means were 0.68 and 0.55. With `height`, `r`
#' is the factor of Seidl, Jensen and Alt (2000), which follows the body
#' water of a person of that weight and height:
#' `0.31608 - 0.004821 W + 0.004632 H` for men and
#' `0.31223 - 0.006446 W + 0.004466 H` for women, `H` in centimetres.
#' Give `r` to use another. `beta` 0.015 is the usual mean; reported
#' rates run from 0.010 to 0.025.
#'
#' The Canadian Postsecondary Education Alcohol and Drug Use Survey
#' derives its `ebac_tot` this way (Health Canada 2022, Table 3): the
#' Seidl factors, 13.6 grams a drink, the drinks and the hours of the
#' heaviest drinking day, with `beta` 0.017; a value of 0.5 or more is
#' set aside there as not survivable, and its `ebac_legal` is a value
#' above 0.08. The guide writes the elimination inside the division by
#' ten, `(A / (r W) - 0.017 t) / 10`, which is 0.0017 grams per 100 mL
#' an hour; the values of the public file fall with the hours at that
#' rate, so `beta = 0.0017` is what reproduces them, and `beta = 0.017`
#' is the rate of Seidl, Jensen and Alt.
#'
#' Errors: an argument that is not numeric; a negative `drinks` or
#' `hours`; a `weight`, `height` or `r` that is not positive (a Seidl
#' factor that is not positive too: a weight and height no person has);
#' `beta` not one number
#' from 0; `grams` not one positive number; a `sex` other than `"male"`
#' or `"female"`; arguments of different lengths (a length of 1 is
#' repeated).
#'
#' @param drinks The number of standard drinks.
#' @param weight The body weight, in `unit`.
#' @param hours The hours since drinking began.
#' @param sex `"male"` or `"female"`, one or one per person: which default
#'   `r` applies. Not used when `r` is given.
#' @param height The height in centimetres, for the factor of Seidl,
#'   Jensen and Alt; `NULL` (the default) for the fixed factors.
#' @param r The Widmark factor, or `NULL` for the factor of `sex` (and of
#'   `height`, when given).
#' @param beta The elimination, in grams per 100 mL per hour.
#' @param grams The grams of ethanol in one standard drink: 14 in the
#'   United States, 13.45 in Canada (see [stddrk()]), 13.6 in the survey
#'   above.
#' @param unit `"kg"` (the default) or `"lb"`: the unit of `weight`.
#' @return A numeric vector: the concentration in grams per 100 mL, never
#'   below 0; `NA` where an input is missing.
#' @references
#'   Widmark, E. M. P. (1932). *Die theoretischen Grundlagen und die
#'   praktische Verwendbarkeit der gerichtlich-medizinischen
#'   Alkoholbestimmung*. Urban & Schwarzenberg, Berlin.
#'
#'   Matthews, D. B. and Miller, W. R. (1979). Estimating blood alcohol
#'   concentration: two computer programs and their applications in
#'   therapy and research. *Addictive Behaviors*, 4(1), 55-60.
#'   \doi{10.1016/0306-4603(79)90021-2}
#'
#'   Seidl, S., Jensen, U. and Alt, A. (2000). The calculation of blood
#'   ethanol concentrations in males and females. *International Journal
#'   of Legal Medicine*, 114(1-2), 71-77. \doi{10.1007/s004140000154}
#'
#'   Health Canada (2022). *Canadian Postsecondary Education Alcohol and
#'   Drug Use Survey, 2021-22: Public Use Microdata File user guide*.
#'   Table 3, derived variables `ebac_tot` and `ebac_legal`.
#' @seealso [stddrk()] for standard drinks, [hvydrk()] and [lrdg()] for
#'   heavy drinking and the low-risk guideline; [etohlr()] for ethanol
#'   found after death.
#' @examples
#' # four drinks in two hours, 180 lb
#' ebac(4, 180, 2, unit = "lb")
#' # the same person at 0, 1, 2 and 3 hours, and a woman of the same weight
#' ebac(4, 81.6, 0:3)
#' ebac(4, 81.6, 2, sex = "female")
#' # the survey's way: body water from weight and height, 13.6 g a drink
#' ebac(5, 62, 3, sex = "female", height = 164, beta = 0.017, grams = 13.6)
#' # over 0.08?
#' ebac(c(2, 5), c(70, 70), 1) >= 0.08
#' try(ebac(4, 0, 2))
#' @export
ebac <- function(
  drinks,
  weight,
  hours = 0,
  sex = "male",
  height = NULL,
  r = NULL,
  beta = 0.015,
  grams = 14,
  unit = c("kg", "lb")
) {
  fn <- "ebac"
  unit <- .ckmarg("unit", fn)
  drinks <- .adxnum(drinks, "drinks", fn, "numbers from 0")
  weight <- .adxnum(weight, "weight", fn, "positive numbers")
  hours <- .adxnum(hours, "hours", fn, "numbers from 0")
  kg <- if (unit == "lb") weight * 0.45359237 else weight
  r <- .adxr(sex, kg, height, r, fn)
  one <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!one(beta) || beta < 0) {
    .ckstop(fn, "`beta` must be one number from 0")
  }
  if (!one(grams) || grams <= 0) {
    .ckstop(fn, "`grams` must be one positive number")
  }
  len <- c(length(drinks), length(weight), length(hours), length(r))
  if (!all(len %in% c(1L, max(len)))) {
    .ckstop(fn, "the arguments must have one length, or length 1")
  }
  pmax(0, grams * drinks / (10 * r * kg) - beta * hours)
}

#' Standard drinks in a volume of a beverage
#'
#' The number of standard drinks in `vol` millilitres of a beverage of
#' `abv` per cent alcohol by volume: the grams of ethanol (the volume of
#' ethanol times its density, 0.789 g per mL at 20 degrees Celsius) over
#' the grams a country calls one drink. Countries differ by a factor of
#' two and more (Kalinowski and Humphreys 2016), so a count of drinks
#' means nothing without the standard it was counted in.
#'
#' Errors: `vol` or `abv` not numeric, `vol` negative, `abv` outside 0 to
#' 100; `std` not one of the standards; `grams` not one positive number;
#' `vol` and `abv` of different lengths (a length of 1 is repeated).
#'
#' @param vol The volume, in millilitres.
#' @param abv The alcohol by volume, in per cent (5 for a beer of 5%).
#' @param std The standard: `"ca"` (Canada, 13.45 g, the default), `"us"`
#'   (United States, 14 g), `"uk"` (a unit, 8 g), `"au"` (Australia, 10 g)
#'   or `"who"` (World Health Organization, 10 g).
#' @param grams The grams of ethanol in one drink, in place of `std`.
#' @return A numeric vector: standard drinks.
#' @references
#'   Kalinowski, A. and Humphreys, K. (2016). Governmental standard drink
#'   definitions and low-risk alcohol consumption guidelines in 37
#'   countries. *Addiction*, 111(7), 1293-1298. \doi{10.1111/add.13341}
#' @seealso [ebac()].
#' @examples
#' # a 341 mL bottle of 5% beer, a 142 mL glass of 12% wine, 43 mL of 40%
#' # spirits: one Canadian standard drink each
#' stddrk(c(341, 142, 43), c(5, 12, 40))
#' # a 750 mL bottle of 13% wine, by four standards
#' vapply(c("ca", "us", "uk", "au"), function(s) stddrk(750, 13, s), 0)
#' try(stddrk(341, 500))
#' @export
stddrk <- function(
  vol,
  abv,
  std = c("ca", "us", "uk", "au", "who"),
  grams = NULL
) {
  fn <- "stddrk"
  std <- .ckmarg("std", fn)
  if (!is.numeric(vol) || !length(vol)) {
    .ckstop(fn, "`vol` must be numeric")
  }
  if (!is.numeric(abv) || !length(abv)) {
    .ckstop(fn, "`abv` must be numeric")
  }
  if (any(vol < 0, na.rm = TRUE)) {
    .ckstop(fn, "`vol` must be numbers from 0")
  }
  if (any(abv < 0 | abv > 100, na.rm = TRUE)) {
    .ckstop(fn, "`abv` must be per cent, from 0 to 100")
  }
  if (is.null(grams)) {
    grams <- c(ca = 13.45, us = 14, uk = 8, au = 10, who = 10)[[std]]
  } else if (
    !is.numeric(grams) || length(grams) != 1L || !is.finite(grams) || grams <= 0
  ) {
    .ckstop(fn, "`grams` must be one positive number")
  }
  n <- max(length(vol), length(abv))
  if (!all(c(length(vol), length(abv)) %in% c(1L, n))) {
    .ckstop(fn, "the arguments must have one length, or length 1")
  }
  as.numeric(vol) * as.numeric(abv) / 100 * 0.789 / grams
}

#' Scores of the AUDIT, the DAST-10, the CUDIT-R and the ASSIST
#'
#' Scores four screening questionnaires from their item answers, one row
#' per respondent and one column per item, in the order of the
#' questionnaire. An answer outside the scoring of its item is refused,
#' not scored: a miscoded column must not pass for a low score.
#'
#' `auditc()`: the Alcohol Use Disorders Identification Test (Saunders et
#' al. 1993). Ten items scored 0 to 4 (items 9 and 10: 0, 2 or 4), total
#' 0 to 40, and the zones of the World Health Organization's manual: 0 to
#' 7 low risk, 8 to 15 hazardous, 16 to 19 harmful, 20 and over possible
#' dependence. The first three items are the AUDIT-C (Bush et al. 1998),
#' 0 to 12; with `sex` its screen is positive from 4 for men and from 3
#' for women (Bradley et al. 2003). Three columns are scored as the
#' AUDIT-C alone.
#'
#' `dast10()`: the Drug Abuse Screening Test in its ten-item form (Skinner
#' 1982). Ten yes or no items, 1 for yes; item 3 ("are you always able to
#' stop using drugs when you want to?") counts when the answer is no.
#' Total 0 to 10: 0 none, 1 to 2 low, 3 to 5 moderate, 6 to 8
#' substantial, 9 to 10 severe.
#'
#' `cudit()`: the Cannabis Use Disorders Identification Test, revised
#' (Adamson et al. 2010). Eight items scored 0 to 4 (item 8: 0, 2 or 4),
#' total 0 to 32: 8 and over hazardous use, 12 and over a possible
#' cannabis use disorder.
#'
#' `canast()`: the cannabis score of the Alcohol, Smoking and Substance
#' Involvement Screening Test of the World Health Organization (WHO
#' ASSIST Working Group 2002), as the Canadian surveys of alcohol and
#' drug use derive it (Health Canada 2022, `can_assist`). Six items. The
#' first four are answered 0 never, 1 once or twice, 2 monthly, 3 weekly,
#' 4 daily or almost daily, over the past three months: use, a strong
#' desire to use, problems from use, and failing what was expected. The
#' last two are answered 0 never, 1 yes but not in the past three months,
#' 2 yes in the past three months: the concern of others, and a failed
#' attempt to cut down. The answers are weighted as the test weights
#' them (0, 2, 3, 4, 6; 0, 3, 4, 5, 6; 0, 4, 5, 6, 7; 0, 5, 6, 7, 8; and
#' 0, 3, 6 for each of the last two), total 0 to 39: 0 to 3 low risk, 4
#' to 26 moderate, 27 and over high (Humeniuk et al. 2008).
#'
#' A respondent with a missing item has missing scores. These are
#' screens: a score says who to assess, not who has a disorder.
#'
#' Errors: `x` not a matrix or data frame of numbers with the columns of
#' the questionnaire (`canast()`: six); an answer outside its item's
#' scoring; `sex` not
#' `"male"` or `"female"`, one or one per row.
#'
#'
#' `canast()` takes a `clpre` list from [clprep()] as `x` too: its
#' elements fill the arguments.
#'
#' @param x The answers: a numeric matrix or data frame, a row per
#'   respondent and a column per item.
#' @param sex `auditc()`: `"male"` or `"female"`, one or one per row, for
#'   the AUDIT-C screen; `NULL` (the default) leaves it out.
#' @return A data frame with a row per respondent. `auditc()`: `total`
#'   (`NA` with three items), `auditc`, `zone` (`"low"`, `"hazardous"`,
#'   `"harmful"`, `"dependence"`; `NA` with three items) and `cpos`, the
#'   AUDIT-C screen (`NA` without `sex`). `dast10()`: `total` and `level`
#'   (`"none"`, `"low"`, `"moderate"`, `"substantial"`, `"severe"`).
#'   `cudit()`: `total` and `level` (`"low"`, `"hazardous"`,
#'   `"disorder"`). `canast()`: `total` and `level` (`"low"`,
#'   `"moderate"`, `"high"`).
#' @references
#'   Saunders, J. B., Aasland, O. G., Babor, T. F., de la Fuente, J. R.
#'   and Grant, M. (1993). Development of the Alcohol Use Disorders
#'   Identification Test (AUDIT): WHO collaborative project on early
#'   detection of persons with harmful alcohol consumption-II.
#'   *Addiction*, 88(6), 791-804.
#'   \doi{10.1111/j.1360-0443.1993.tb02093.x}
#'
#'   Bush, K., Kivlahan, D. R., McDonell, M. B., Fihn, S. D. and Bradley,
#'   K. A. (1998). The AUDIT alcohol consumption questions (AUDIT-C): an
#'   effective brief screening test for problem drinking. *Archives of
#'   Internal Medicine*, 158(16), 1789-1795.
#'   \doi{10.1001/archinte.158.16.1789}
#'
#'   Bradley, K. A., Bush, K. R., Epler, A. J., Dobie, D. J., Davis, T.
#'   M., Sporleder, J. L., Maynard, C., Burman, M. L. and Kivlahan, D. R.
#'   (2003). Two brief alcohol-screening tests from the Alcohol Use
#'   Disorders Identification Test (AUDIT). *Archives of Internal
#'   Medicine*, 163(7), 821-829. \doi{10.1001/archinte.163.7.821}
#'
#'   Skinner, H. A. (1982). The drug abuse screening test. *Addictive
#'   Behaviors*, 7(4), 363-371. \doi{10.1016/0306-4603(82)90005-3}
#'
#'   Adamson, S. J., Kay-Lambkin, F. J., Baker, A. L., Lewin, T. J.,
#'   Thornton, L., Kelly, B. J. and Sellman, J. D. (2010). An improved
#'   brief measure of cannabis misuse: the Cannabis Use Disorders
#'   Identification Test-Revised (CUDIT-R). *Drug and Alcohol
#'   Dependence*, 110(1-2), 137-143.
#'   \doi{10.1016/j.drugalcdep.2010.02.017}
#'
#'   WHO ASSIST Working Group (2002). The Alcohol, Smoking and Substance
#'   Involvement Screening Test (ASSIST): development, reliability and
#'   feasibility. *Addiction*, 97(9), 1183-1194.
#'   \doi{10.1046/j.1360-0443.2002.00185.x}
#'
#'   Humeniuk, R., Ali, R., Babor, T. F., Farrell, M., Formigoni, M. L.,
#'   Jittiwutikarn, J., de Lacerda, R. B., Ling, W., Marsden, J.,
#'   Monteiro, M., Nhiwatiwa, S., Pal, H., Poznyak, V. and Simon, S.
#'   (2008). Validation of the Alcohol, Smoking and Substance Involvement
#'   Screening Test (ASSIST). *Addiction*, 103(6), 1039-1047.
#'   \doi{10.1111/j.1360-0443.2007.02114.x}
#'
#'   Health Canada (2022). *Canadian Postsecondary Education Alcohol and
#'   Drug Use Survey, 2021-22: Public Use Microdata File user guide*.
#'   Table 3, derived variable `can_assist`.
#' @seealso [ebac()], [stddrk()], [hvydrk()] and [lrdg()].
#' @examples
#' # made-up answers of three respondents
#' a <- rbind(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'            c(3, 2, 2, 1, 1, 0, 1, 1, 0, 2),
#'            c(4, 4, 4, 3, 3, 2, 3, 3, 4, 4))
#' auditc(a, sex = c("female", "male", "male"))
#' # the first three items alone are the AUDIT-C
#' auditc(a[, 1:3], sex = "female")
#' dast10(rbind(c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
#'              c(1, 1, 0, 1, 1, 0, 1, 1, 0, 1)))
#' cudit(rbind(c(1, 0, 0, 0, 0, 0, 0, 0),
#'             c(4, 3, 2, 1, 2, 1, 2, 4)))
#' # weekly use, a monthly desire, no problems: moderate risk
#' canast(rbind(c(1, 0, 0, 0, 0, 0), c(3, 2, 0, 0, 0, 0),
#'              c(4, 4, 3, 3, 2, 2)))
#' try(auditc(cbind(a[, 1:9], 3)))
#' @export
auditc <- function(x, sex = NULL) {
  fn <- "auditc"
  m <- .adxmat(x, c(10L, 3L), "10 columns (the AUDIT) or 3 (the AUDIT-C)", fn)
  .adxrng(m, c(9L, 10L), fn)
  n <- nrow(m)
  if (!is.null(sex)) {
    if (
      !is.character(sex) ||
        anyNA(sex) ||
        !all(sex %in% c("male", "female")) ||
        !length(sex) %in% c(1L, n)
    ) {
      .ckstop(fn, "`sex` must be \"male\" or \"female\", one or one per row")
    }
  }
  full <- ncol(m) == 10L
  c3 <- as.numeric(rowSums(m[, 1:3, drop = FALSE]))
  total <- if (full) as.numeric(rowSums(m)) else rep(NA_real_, n)
  zone <- c("low", "hazardous", "harmful", "dependence")[
    findInterval(total, c(0, 8, 16, 20))
  ]
  cpos <- if (is.null(sex)) {
    rep(NA, n)
  } else {
    c3 >= ifelse(rep_len(sex, n) == "male", 4, 3)
  }
  data.frame(
    total = total,
    auditc = c3,
    zone = zone,
    cpos = cpos,
    stringsAsFactors = FALSE
  )
}

#' @rdname auditc
#' @export
dast10 <- function(x) {
  fn <- "dast10"
  m <- .adxmat(x, 10L, "10 columns (the DAST-10)", fn)
  bad <- which(!is.na(m) & m != 0 & m != 1, arr.ind = TRUE)
  if (nrow(bad)) {
    .ckstop(fn, "item %d must be 0 (no) or 1 (yes)", min(bad[, "col"]))
  }
  m[, 3L] <- 1 - m[, 3L]
  total <- as.numeric(rowSums(m))
  data.frame(
    total = total,
    level = c("none", "low", "moderate", "substantial", "severe")[
      findInterval(total, c(0, 1, 3, 6, 9))
    ],
    stringsAsFactors = FALSE
  )
}

#' @rdname auditc
#' @export
cudit <- function(x) {
  fn <- "cudit"
  m <- .adxmat(x, 8L, "8 columns (the CUDIT-R)", fn)
  .adxrng(m, 8L, fn)
  total <- as.numeric(rowSums(m))
  data.frame(
    total = total,
    level = c("low", "hazardous", "disorder")[
      findInterval(total, c(0, 8, 12))
    ],
    stringsAsFactors = FALSE
  )
}

#' @rdname auditc
#' @export
canast <- function(x) {
  if (inherits(x, "clpre")) {
    return(.prep(canast, x, match.call(), parent.frame()))
  }
  fn <- "canast"
  m <- .adxmat(x, 6L, "6 columns (the ASSIST items on cannabis)", fn)
  w <- list(
    c(0, 2, 3, 4, 6),
    c(0, 3, 4, 5, 6),
    c(0, 4, 5, 6, 7),
    c(0, 5, 6, 7, 8),
    c(0, 3, 6),
    c(0, 3, 6)
  )
  for (j in 1:6) {
    v <- m[, j]
    top <- length(w[[j]]) - 1L
    if (!all(v[!is.na(v)] %in% 0:top)) {
      .ckstop(fn, "item %d must be answered 0 to %d", j, top)
    }
    m[, j] <- w[[j]][v + 1L]
  }
  total <- as.numeric(rowSums(m))
  data.frame(
    total = total,
    level = c("low", "moderate", "high")[findInterval(total, c(0, 4, 27))],
    stringsAsFactors = FALSE
  )
}

#' Heavy drinking, and drinking above the low-risk guideline
#'
#' Two flags from a number of standard drinks and a sex, as the Canadian
#' surveys of alcohol and drug use derive them (Health Canada 2022).
#' `hvydrk()`: heavy drinking, four or more drinks for women and five or
#' more for men on one occasion, the measure of Wechsler et al. (1995).
#' `lrdg()`: more drinks in a day than the daily limit of Canada's
#' Low-Risk Alcohol Drinking Guidelines for chronic harm, two for women
#' and three for men (Stockwell et al. 2012).
#'
#' The first counts the threshold itself as heavy (`>=`); the second is
#' over the limit only above it (`>`). Both thresholds are arguments, for
#' another definition or another country's guideline. The guidelines of
#' 2011 that `lrdg()` follows were replaced in 2023 by Canada's Guidance
#' on Alcohol and Health, which gives a continuum of risk by drinks a
#' week and no daily limit by sex; the survey variables, and so the
#' defaults here, are those of 2011.
#'
#' Errors: `drinks` not numeric or negative; `sex` not `"male"` or
#' `"female"`; `f` or `m` not one number from 0; `drinks` and `sex` of
#' different lengths (a length of 1 is repeated).
#'
#'
#' A `clpre` list from [clprep()] may be given as `drinks`: its elements
#' fill the arguments.
#'
#' @param drinks The number of standard drinks: on one occasion
#'   (`hvydrk()`), in a day (`lrdg()`).
#' @param sex `"male"` or `"female"`, one or one per value of `drinks`.
#' @param f,m The thresholds for women and for men.
#' @return A logical vector: `TRUE` where the drinking is heavy
#'   (`hvydrk()`) or over the limit (`lrdg()`); `NA` where `drinks` is
#'   missing.
#' @references
#'   Wechsler, H., Dowdall, G. W., Davenport, A. and Rimm, E. B. (1995). A
#'   gender-specific measure of binge drinking among college students.
#'   *American Journal of Public Health*, 85(7), 982-985.
#'   \doi{10.2105/AJPH.85.7.982}
#'
#'   Stockwell, T., Butt, P., Beirness, D., Gliksman, L. and Paradis, C.
#'   (2012). The basis for Canada's new low-risk drinking guidelines: a
#'   relative risk approach to estimating hazardous levels and patterns
#'   of alcohol use. *Drug and Alcohol Review*, 31(2), 126-134.
#'   \doi{10.1111/j.1465-3362.2011.00342.x}
#'
#'   Health Canada (2022). *Canadian Postsecondary Education Alcohol and
#'   Drug Use Survey, 2021-22: Public Use Microdata File user guide*.
#'   Table 3, derived variables `alc12_30d_prev`, `u_d_lrdg` and
#'   `alc13_lrdg`.
#' @seealso [ebac()] and [stddrk()]; [auditc()] for the screening
#'   questionnaires.
#' @examples
#' d <- c(2, 3, 4, 5, NA)
#' hvydrk(d, "female")
#' hvydrk(d, "male")
#' lrdg(d, c("female", "male", "male", "female", "male"))
#' # another definition: six drinks for everyone
#' hvydrk(d, "female", f = 6, m = 6)
#' try(hvydrk(-1, "male"))
#' # a survey's public file, by its preset: the codes 997 to 999 are not
#' # answers, and the sex is coded 1 and 2
#' d <- data.frame(alc13a = c(3, 5, 997, 2), dvdemq01 = c(1, 2, 1, 2))
#' lrdg(clprep(d, "cpads-lrdg"))
#' @export
hvydrk <- function(drinks, sex, f = 4, m = 5) {
  if (inherits(drinks, "clpre")) {
    return(.prep(hvydrk, drinks, match.call(), parent.frame()))
  }
  .adxflag(drinks, sex, f, m, "hvydrk", `>=`)
}

#' @rdname hvydrk
#' @export
lrdg <- function(drinks, sex, f = 2, m = 3) {
  if (inherits(drinks, "clpre")) {
    return(.prep(lrdg, drinks, match.call(), parent.frame()))
  }
  .adxflag(drinks, sex, f, m, "lrdg", `>`)
}

#' drinks against the threshold of each sex by `op`, or the refusal of fn.
#' @noRd
.adxflag <- function(drinks, sex, f, m, fn, op) {
  if (!is.numeric(drinks) || !length(drinks)) {
    .ckstop(fn, "`drinks` must be numeric")
  }
  if (any(drinks < 0, na.rm = TRUE)) {
    .ckstop(fn, "`drinks` must be numbers from 0")
  }
  if (
    !is.character(sex) ||
      !length(sex) ||
      anyNA(sex) ||
      !all(sex %in% c("male", "female"))
  ) {
    .ckstop(fn, "`sex` must be \"male\" or \"female\"")
  }
  one <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 0
  if (!one(f)) {
    .ckstop(fn, "`f` must be one number from 0")
  }
  if (!one(m)) {
    .ckstop(fn, "`m` must be one number from 0")
  }
  n <- max(length(drinks), length(sex))
  if (!all(c(length(drinks), length(sex)) %in% c(1L, n))) {
    .ckstop(fn, "the arguments must have one length, or length 1")
  }
  op(rep_len(as.numeric(drinks), n), ifelse(rep_len(sex, n) == "male", m, f))
}

#' The answers as a numeric matrix with one of `k` numbers of columns, or
#' the refusal of fn.
#' @noRd
.adxmat <- function(x, k, what, fn) {
  if (is.data.frame(x) && all(vapply(x, is.numeric, NA))) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || !ncol(x) %in% k) {
    .ckstop(fn, "`x` must be a numeric matrix or data frame of %s", what)
  }
  matrix(as.numeric(x), nrow(x), ncol(x))
}

#' Items scored 0 to 4, the items `two` 0, 2 or 4, or the refusal of fn
#' (the first item out of its scoring).
#' @noRd
.adxrng <- function(m, two, fn) {
  for (j in seq_len(ncol(m))) {
    v <- m[, j]
    v <- v[!is.na(v)]
    if (j %in% two) {
      if (!all(v %in% c(0, 2, 4))) {
        .ckstop(fn, "item %d must be scored 0, 2 or 4", j)
      }
    } else if (!all(v %in% 0:4)) {
      .ckstop(fn, "item %d must be scored 0 to 4", j)
    }
  }
  invisible(NULL)
}

#' Numbers that are `what` ("numbers from 0": none negative; "positive
#' numbers": none 0 or negative), or the refusal of fn.
#' @noRd
.adxnum <- function(v, name, fn, what) {
  if (!is.numeric(v) || !length(v)) {
    .ckstop(fn, "`%s` must be numeric", name)
  }
  bad <- if (what == "positive numbers") v <= 0 else v < 0
  if (any(bad, na.rm = TRUE)) {
    .ckstop(fn, "`%s` must be %s", name, what)
  }
  as.numeric(v)
}

#' The Widmark factor of ebac(): `r` as given, or that of `sex`, fixed or,
#' with `height`, by Seidl, Jensen and Alt from the weight `kg`.
#' @noRd
.adxr <- function(sex, kg, height, r, fn) {
  if (!is.null(height)) {
    height <- .adxnum(height, "height", fn, "positive numbers")
  }
  if (!is.null(r)) {
    return(.adxnum(r, "r", fn, "positive numbers"))
  }
  if (
    !is.character(sex) ||
      !length(sex) ||
      anyNA(sex) ||
      !all(sex %in% c("male", "female"))
  ) {
    .ckstop(fn, "`sex` must be \"male\" or \"female\"")
  }
  if (is.null(height)) {
    return(ifelse(sex == "male", 0.73, 0.66))
  }
  len <- c(length(sex), length(kg), length(height))
  if (!all(len %in% c(1L, max(len)))) {
    .ckstop(fn, "the arguments must have one length, or length 1")
  }
  r <- ifelse(
    sex == "male",
    0.31608 - 0.004821 * kg + 0.004632 * height,
    0.31223 - 0.006446 * kg + 0.004466 * height
  )
  if (any(r <= 0, na.rm = TRUE)) {
    .ckstop(fn, "the weight and height give a factor that is not positive")
  }
  r
}
