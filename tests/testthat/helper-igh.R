# Made-up tables in the layouts of Statistics Canada's published tables,
# for the tests of the Indigenous health, census and Calls to Action
# functions. The numbers are invented: they follow the layouts (columns,
# labels, the profile's indentation and quality flags), not any published
# value.

# The population table of the Indigenous Population Profile: two places
# (one under 250 people), three residences, four age groups, twelve
# characteristics and the nine identity groups with their SYMBOL_ columns.
.ipp <- function() {
  ch <- data.frame(
    id = c(1, 13, 20, 37, 38, 545, 546, 560, 629, 631, 807, 809, 920, 923),
    name = c(
      "Total - Age groups of the population - 100% data",
      paste(
        "Total - Household and family characteristics for the population",
        "in private households - 25% sample data"
      ),
      "  Total - Children in census families",
      "  Total - Persons not in census families",
      "    Foster children",
      paste(
        "Total - Highest certificate, diploma or degree for the population",
        "aged 25 years and over in private households - 25% sample data"
      ),
      "  No certificate, diploma or degree",
      "  Bachelor's degree or higher",
      "Participation rate",
      "Unemployment rate",
      paste(
        "Number of total income recipients aged 15 years and over in",
        "private households in 2020 - 25% sample data"
      ),
      "Average total income in 2020 among recipients ($)",
      paste(
        "Total - Low-income status in 2020 for the population in private",
        "households - 25% sample data"
      ),
      paste(
        "Prevalence of low income based on the Low-income measure, after",
        "tax (LIM-AT) (%)"
      )
    ),
    # each characteristic's size against the population
    f = c(
      1,
      0.97,
      0.6,
      0.3,
      0.02,
      0.55,
      0.12,
      0.2,
      63,
      9,
      0.75,
      41000,
      0.97,
      18
    ),
    rate = c(rep(FALSE, 8), TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  geo <- data.frame(
    name = c("Canada", "Lakeside"),
    dguid = c("2021A000011124", "2021A0005000001"),
    lvl = c("Country", "Census subdivision"),
    size = c(40000, 150),
    tnr = c(3.1, 7.4)
  )
  ages <- c(
    "Total - Age groups",
    "0 to 14 years",
    "25 years and over",
    "25 to 64 years"
  )
  af <- c(1, 0.25, 0.55, 0.4)
  res <- c("Total - Residence on or off reserve", "On reserve", "Off reserve")
  rf <- c(1, 0.4, 0.6)
  grp <- c(
    "TOTAL_IND_IDENT (1)",
    "IND_IDENT (2)",
    "SINGLE_IND_RESP (3)",
    "FIRST_NATIONS",
    "M\u00c9TIS",
    "INUK",
    "MULTI_IND_RESP (4)",
    "IND_RESP_NIE (5)",
    "NON_IND"
  )
  gf <- c(20, 1, 0.97, 0.58, 0.34, 0.04, 0.016, 0.019, 19)
  g <- expand.grid(
    c = seq_len(nrow(ch)),
    a = seq_along(ages),
    r = seq_along(res),
    p = seq_len(nrow(geo))
  )
  d <- data.frame(
    CENSUS_YEAR = "2021",
    DGUID = geo$dguid[g$p],
    GEO_LEVEL = geo$lvl[g$p],
    GEO_NAME = geo$name[g$p],
    TNR_SF = geo$tnr[g$p] - 1,
    TNR_LF = geo$tnr[g$p],
    GENDER = "Total - Gender",
    AGE_GROUP = ages[g$a],
    RESIDENCE = res[g$r],
    CHARACTERISTIC_ID = as.character(ch$id[g$c]),
    CHARACTERISTIC_NAME = ch$name[g$c],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(grp)) {
    v <- ifelse(
      ch$rate[g$c],
      ch$f[g$c] + k / 10,
      5 * round(geo$size[g$p] * gf[k] * af[g$a] * rf[g$r] * ch$f[g$c] / 5)
    )
    d[[grp[k]]] <- as.character(v)
    d[[paste0("SYMBOL_", sub(" \\([0-9]+\\)$", "", grp[k]))]] <- ""
  }
  # one cell suppressed, one to use with caution
  d$SYMBOL_INUK[d$GEO_NAME == "Lakeside" & d$CHARACTERISTIC_ID == "38"] <- "x"
  d$SYMBOL_FIRST_NATIONS[
    d$CHARACTERISTIC_ID == "560" &
      d$GEO_NAME == "Lakeside"
  ] <- "E"
  d
}

# The household table of the profile, for the housing component.
.ihh <- function() {
  ch <- data.frame(
    id = c(1, 37, 38, 52, 54),
    name = c(
      "Total - Private households by household size - 25% sample data",
      paste(
        "Total - Private households by number of persons per room - 25%",
        "sample data"
      ),
      "  One person or fewer per room",
      paste(
        "Total - Occupied private dwellings by dwelling condition - 25%",
        "sample data"
      ),
      "  Major repairs needed"
    ),
    f = c(1, 1, 0.9, 1, 0.15)
  )
  geo <- c("Canada", "Lakeside")
  dg <- c("2021A000011124", "2021A0005000001")
  size <- c(16000, 60)
  g <- expand.grid(c = seq_len(nrow(ch)), p = 1:2)
  d <- data.frame(
    DGUID = dg[g$p],
    GEO_NAME = geo[g$p],
    CHARACTERISTIC_ID = as.character(ch$id[g$c]),
    CHARACTERISTIC_NAME = ch$name[g$c],
    stringsAsFactors = FALSE
  )
  hf <- c(1, 0.06, 0.94)
  hn <- c("TOTAL_IND_HHLD (1)", "IND_HHLD (2)", "NON_IND_HHLD")
  for (k in 1:3) {
    d[[hn[k]]] <- as.character(
      5 * round(size[g$p] * hf[k] * ch$f[g$c] * (1 - 0.03 * k) / 5)
    )
    d[[paste0("SYMBOL_", sub(" \\([0-9]+\\)$", "", hn[k]))]] <- ""
  }
  d
}

# Health indicators in the long layout of scests(): two periods, an
# all-adults row and one age group, four groups.
.ih <- function() {
  ind <- c(
    "Perceived mental health, very good or excellent",
    "Mood disorder",
    "Anxiety disorder",
    "Diabetes",
    "Heavy drinking",
    "Has a regular healthcare provider",
    "Sense of belonging to local community, somewhat strong or very strong",
    "Perceived life stress, most days quite a bit or extremely stressful"
  )
  base <- c(55, 12, 13, 7, 20, 80, 70, 22)
  grp <- c("First Nations", "M\u00e9tis", "Inuit", "Non-Indigenous")
  gs <- c(-8, -4, -6, 0)
  g <- expand.grid(i = seq_along(ind), g = seq_along(grp), a = 1:2, p = 1:2)
  est <- base[g$i] +
    gs[g$g] * ifelse(base[g$i] > 50, 1, -0.3) +
    g$p -
    1 +
    (g$a - 1) * 2
  half <- 1 + 0.5 * (g$g != 4) + 0.3 * (g$a == 2)
  q <- ifelse(g$g == 3 & g$i == 4, "E", "ok")
  q[g$g == 3 & g$i == 5 & g$p == 2] <- "F"
  data.frame(
    period = c("2015/2018", "2019/2022")[g$p],
    group = grp[g$g],
    age = c("Total, 18 years and over", "18 to 34 years")[g$a],
    indicator = ind[g$i],
    estimate = ifelse(q == "F", NA, est),
    lower = ifelse(q == "F", NA, est - half),
    upper = ifelse(q == "F", NA, est + half),
    quality = q,
    stringsAsFactors = FALSE
  )
}

# Life expectancy, infant mortality and suicidal thoughts in the same layout.
.ile <- function() {
  data.frame(
    period = rep(c("2011/2016", "2016/2021"), each = 3),
    group = rep(c("First Nations", "Inuit", "Non-Indigenous"), 2),
    estimate = c(73, 70, 81.5, 74.2, 70.9, 81.8),
    lower = c(72.5, 69, 81.3, 73.6, 69.1, 81.6),
    upper = c(73.5, 71, 81.7, 74.8, 72.7, 82),
    quality = "ok",
    stringsAsFactors = FALSE
  )
}

.iim <- function() {
  data.frame(
    period = rep(c("1994/1996", "2014/2016"), each = 3),
    group = rep(c("First Nations", "Inuit", "Non-Indigenous"), 2),
    estimate = c(NA, NA, NA, 9.2, 12.3, 4.4),
    lower = c(NA, NA, NA, 8.1, 9.5, 4.2),
    upper = c(NA, NA, NA, 10.3, 15.1, 4.6),
    quality = c("..", "..", "..", "ok", "E", "ok"),
    stringsAsFactors = FALSE
  )
}

.isu <- function() {
  data.frame(
    period = "2017",
    group = c("First Nations", "M\u00e9tis", "Inuit"),
    age = "Total, 15 years and over",
    measure = "Seriously contemplated suicide in the past 12 months",
    estimate = c(4.1, 3.2, 4.5),
    lower = c(3.6, 2.7, 3.6),
    upper = c(4.6, 3.7, 5.4),
    quality = "ok",
    stringsAsFactors = FALSE
  )
}

# A result's printed lines as one string.
ls_flat <- function(x) gsub("\\s+", " ", paste(format(x), collapse = " "))

# The numeric and logical columns (or elements) of a result, without names,
# as the reference values of test-igv.R hold them.
.ignum <- function(x) {
  if (is.data.frame(x)) {
    x <- as.list(x)
  }
  if (is.list(x)) {
    x <- x[vapply(x, function(v) is.numeric(v) || is.logical(v), NA)]
    return(unname(lapply(x, function(v) unname(as.vector(v)))))
  }
  unname(as.vector(x))
}

# The tables of the source package's tests.
.hc <- function() {
  data.frame(
    g = rep(c("A", "B"), each = 4),
    ind = rep(c("well", "well", "ill", "ill"), 2),
    per = rep(c("t1", "t2"), 4),
    est = c(50, 55, 20, 18, 40, 41, 30, 35),
    lo = c(46, 51, 18, 16, 36, 37, 27, 32),
    hi = c(54, 59, 22, 20, 44, 45, 33, 38),
    q = c("", "", "", "", "", "F", "", "E"),
    stringsAsFactors = FALSE
  )
}


# shares and Wald limits that imply an effective sample size of exactly n
.wald <- function(p, n) {
  z <- stats::qnorm(0.975)
  s <- sqrt(p * (1 - p) / n)
  list(lo = 100 * (p - z * s), hi = 100 * (p + z * s))
}

.wb <- function() {
  data.frame(
    group = rep(c("First Nations", "Metis"), each = 4),
    indicator = rep(c("Hopeful", "Total", "Distress", "Satisfied"), 2),
    pct = c(60, 100, 20, 50, 65, 100, 15, 55),
    lo = c(55, 100, 17, 45, 60, 100, 12, 50),
    hi = c(65, 100, 23, 55, 70, 100, 18, 60),
    q = c("", "", "", "E", "", "", "", ""),
    stringsAsFactors = FALSE
  )
}
