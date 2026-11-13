# Parity with the functions crimlab's generic methods were ported from
# (rmorie, rmoriebricklayer). Each case is a function returning the source
# function's numbers on a small synthetic input defined here; attribute
# pkg names the package it needs. portref(key) runs the case when that
# package is installed and otherwise returns the value recorded in
# ref-port.rds by dev/ref-port.R, which sources this file.

pcase <- list()
pc <- function(key, pkg, f) {
  attr(f, "pkg") <- pkg
  pcase[[key]] <<- f
}

portref <- function(key) {
  f <- pcase[[key]]
  if (requireNamespace(attr(f, "pkg"), quietly = TRUE)) return(f())
  readRDS(test_path("ref-port.rds"))[[key]]
}

# rates: an age-by-group table, a two-period count, SIR inputs
pin <- list(
  dsr = data.frame(age = rep(c("y", "m", "o"), 2),
                   grp = rep(c("A", "B"), each = 3),
                   x = c(30, 22, 5, 41, 30, 12),
                   n = c(5000, 9000, 6000, 3000, 8000, 7000)),
  chg = data.frame(year = c(2020, 2021), count = c(30, 45), pop = c(1e4, 1.1e4)),
  obs = c(15, 39, 6, 30), expd = c(14, 20, 9, 31),
  tab = c(a = 15, b = 45, c = 40, d = 40),
  uof = data.frame(race = rep(c("W", "B", "E"), c(40, 25, 15)),
                   force = c(rep(c(1, 0), c(8, 32)), rep(c(1, 0), c(10, 15)),
                             rep(c(1, 0), c(3, 12)))),
  siu = data.frame(police_service = rep(c("North PS", "Lake PS"), c(7, 5)),
                   charges_recommended = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
                                           FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
)

blf <- function(f) getExportedValue("rmoriebricklayer", f)
rmf <- function(f) getExportedValue("rmorie", f)

pc("agestd", "rmoriebricklayer", function() {
  d <- pin$dsr
  a <- blf("agestd")(d$x, d$n, std = c(30, 45, 25), grp = d$grp, per = 1e4)
  g <- blf("rategap")(a, "B", "A")
  list(rate = a$rate, lower = a$lower, upper = a$upper, ratio = g$ratio,
       ratiolo = g$ratiolo, ratioup = g$ratioup, diff = g$difference,
       diflo = g$diflo, difup = g$difup)
})
pc("stdinc", "rmoriebricklayer", function() {
  s <- blf("stdinc")(pin$obs, pin$expd)
  e <- blf("ebrates")(pin$obs, pin$expd)
  f <- blf("fnllim")(pin$expd, lvls = 0.95)
  list(lower = s$lower, upper = s$upper, eb = e$eb, shrinkage = e$shrinkage,
       flo = f$lower, fhi = f$upper)
})
pc("expcnts", "rmoriebricklayer", function() {
  d <- pin$dsr
  blf("expcnts")(d$x, d$n, d$grp, d$age)$expected
})
pc("rtcalc", "rmoriebricklayer", function() {
  x <- pin$chg
  r <- blf("rtcalc")(x, count, pop, per = 1000)
  k <- blf("ratechg")(x, count, pop, year)
  list(rate = r$rate, lower = r$lower, upper = r$upper, rr = k$rateratio[2],
       pctlo = k$pctlo[2], pctup = k$pctup[2])
})
pc("riskratio", "rmorie", function() {
  t <- pin$tab
  a <- rmf("risk_ratio")(t[["a"]], t[["b"]], t[["c"]], t[["d"]])
  b <- rmf("morie_risk_ratio_ci")(matrix(t, 2, byrow = TRUE))
  d <- rmf("morie_risk_difference_ci")(matrix(t, 2, byrow = TRUE))
  w <- rmf("risk_difference")(t[["a"]], t[["b"]], t[["c"]], t[["d"]])
  z <- rmf("risk_ratio")(0, 10, 5, 5)
  list(rr = a$estimate, lo = a$ci_lower, hi = a$ci_upper, rr2 = b$rr,
       lo2 = b$ci_lower, hi2 = b$ci_upper, rd = d$rd, rdlo = d$ci_lower,
       rdhi = d$ci_upper, wrd = w$estimate, wlo = w$ci_lower, whi = w$ci_upper,
       zrr = z$estimate, zlo = z$ci_lower)
})
pc("uofdisp", "rmorie", function() {
  r <- rmf("mrm_uof_demographic_disparity")(pin$uof, "race", "force")
  p <- r$per_category
  list(category = vapply(p, function(e) e$category, ""),
       rate = vapply(p, function(e) e$rate, 0),
       lo = vapply(p, function(e) e$lo, 0), hi = vapply(p, function(e) e$hi, 0),
       rr = vapply(p, function(e) e$rr, 0), baseline = r$baseline)
})
pc("siusvc", "rmorie", function() {
  r <- rmf("morie_siu_by_police_service")(pin$siu)
  list(counts = unlist(r$payload$counts), charged = r$payload$charged,
       rows = r$tables[[1]]$rows)
})
pc("tpsuof", "rmorie", function() {
  r <- rmf("morie_tps_use_of_force")(c("a", "b", "a", NA, "c"), 250)
  list(rate = r$rate, n = r$n, counts = unlist(r$type_counts))
})

# contingency: OTIS-shaped aggregate tables and record-level layouts
local({
  # the draws leave the caller's random stream as it was
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, globalenv()))
  set.seed(11)
  reg <- c("Central", "East", "North", "Toronto", "West")
  g <- expand.grid(EndFiscalYear = 2021:2022, Race = c("Black", "Indigenous", "White", "Other"),
                   Region_MostRecentPlacement = reg, stringsAsFactors = FALSE)
  g$NumberIndividuals_RestrictiveConfinement <- rpois(nrow(g), 25)
  pin$c04 <<- g
  cs <- expand.grid(EndFiscalYear = 2021:2022, Race = c("Indigenous", "White"),
                    Gender = c("Female", "Male"), Religion = c("Christian", "Muslim", "None"),
                    Age_Category = c("18-24", "25-49", "50+"),
                    Alert_Type = c("Mental Health", "Suicide Risk", "No Alert"),
                    Region_MostRecentPlacement = reg[1:3], stringsAsFactors = FALSE)
  cs$NumberIndividuals_RestrictiveConfinement <- rpois(nrow(cs), 6)
  pin$cser <<- cs
  b4 <- expand.grid(EndFiscalYear = 2021:2022, Region_AtTimeOfPlacement = reg,
                    Gender = c("Female", "Male"), Measure = c("Max", "Median", "Mode"),
                    stringsAsFactors = FALSE)
  b4$NumberConsecutiveDays_Segregation <- rpois(nrow(b4), 9)
  pin$b04 <<- b4
  d2 <- expand.grid(Year = 2015:2020, Gender = c("Female", "Male"), stringsAsFactors = FALSE)
  d2$Number_CustodialDeaths <- rpois(nrow(d2), 4)
  pin$d02 <<- d2
  n <- 300
  at <- sample(reg, n, TRUE)
  now <- ifelse(runif(n) < 0.7, at, sample(reg, n, TRUE))
  pin$loc <<- data.frame(Region_AtTimeOfPlacement = at, Region_MostRecentPlacement = now)
  pin$loc2 <<- data.frame(at = sample(c("N", "S"), 60, TRUE), now = sample(c("N", "S"), 60, TRUE))
  pin$wpn <<- data.frame(weapon = sample(c("none", "knife", "gun", "other"), 200, TRUE,
                                         c(0.5, 0.2, 0.2, 0.1)),
                         force = sample(c("hands", "CEW", "firearm"), 200, TRUE))
  yn <- function() sample(c("Yes", "No"), 150, TRUE)
  pin$mort <<- data.frame(MentalHealth_Alert = yn(), SuicideRisk_Alert = yn(),
                          SuicideWatch_Alert = yn())
  pin$siuy <<- data.frame(
    date_of_incident_iso = sprintf("%d-03-01", sample(2018:2022, 120, TRUE)),
    charges_recommended = sample(c(TRUE, FALSE), 120, TRUE, c(0.3, 0.7)))
  pin$siuo <<- data.frame(
    police_service = sample(c("North PS", "Lake PS", "City PS"), 90, TRUE),
    director_decision_category = sample(c("charges", "no charges", "other"), 90, TRUE))
})

pc("otisxtab", "rmorie", function() {
  rows <- function(r) {
    t <- r$tables[[length(r$tables)]]
    list(head = t$headers, rows = lapply(t$rows, as.character))
  }
  list(c04 = rows(rmf("morie_otis_analyze_c04")(pin$c04)),
       b04 = rows(rmf("morie_otis_analyze_b04")(pin$b04)),
       d02 = lapply(rmf("morie_otis_analyze_d02")(pin$d02)$tables,
                    function(t) lapply(t$rows, as.character)))
})
pc("cchi2", "rmorie", function() {
  ds <- list(c03 = pin$cser, c04 = pin$cser, c05 = pin$cser, c06 = pin$cser,
             c07 = pin$cser, c08 = pin$cser, c09 = pin$cser)
  r <- rmf("morie_otis_analyze_c_chi2")(ds)
  lapply(r$payload$slices, function(s) unlist(s$stats[c("chi2", "dof", "pvalue", "cramer_v")]))
})
pc("locality", "rmorie", function() {
  a <- rmf("mrm_otis_region_locality")(pin$loc)
  b <- rmf("mrm_uof_region_locality")(pin$loc2, "at", "now")
  c2 <- rmf("mrm_uof_region_locality")(pin$loc, "Region_AtTimeOfPlacement",
                                      "Region_MostRecentPlacement")
  list(otis = c(a$chi2, a$df, a$p_value, a$cramers_v, a$diagonal_share),
       uof2 = c(b$chi2, b$pvalue, b$cramers_v, b$diagonal_share),
       uof5 = c(c2$chi2, c2$pvalue, c2$cramers_v, c2$diagonal_share))
})
pc("weapon", "rmorie", function() {
  w <- rmf("mrm_uof_weapon_diversity")(pin$wpn, "weapon", "force")
  list(stat = c(w$chi2, w$df, w$pvalue, w$cramers_v),
       top = vapply(w$top_residuals, function(e) e$std_residual, 0))
})
pc("mort", "rmorie", function() {
  as.list(rmf("mrm_otis_mortification_cooccurrence")(pin$mort))
})
pc("siuchi", "rmorie", function() {
  r <- rmf("morie_siu_charges_by_year_chi2")(pin$siuy)
  o <- rmf("mrm_siu_outcome_classifier")(pin$siuo)
  y <- rmf("morie_siu_by_year")(pin$siuy)
  list(stat = unlist(r$payload), out = o,
       years = unlist(y$payload$by_year))
})
