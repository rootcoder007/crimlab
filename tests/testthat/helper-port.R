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

bl <- function(f) getExportedValue("rmoriebricklayer", f)
rm <- function(f) getExportedValue("rmorie", f)

pc("agestd", "rmoriebricklayer", function() {
  d <- pin$dsr
  a <- bl("agestd")(d$x, d$n, std = c(30, 45, 25), grp = d$grp, per = 1e4)
  g <- bl("rategap")(a, "B", "A")
  list(rate = a$rate, lower = a$lower, upper = a$upper, ratio = g$ratio,
       ratiolo = g$ratiolo, ratioup = g$ratioup, diff = g$difference,
       diflo = g$diflo, difup = g$difup)
})
pc("stdinc", "rmoriebricklayer", function() {
  s <- bl("stdinc")(pin$obs, pin$expd)
  e <- bl("ebrates")(pin$obs, pin$expd)
  f <- bl("fnllim")(pin$expd, lvls = 0.95)
  list(lower = s$lower, upper = s$upper, eb = e$eb, shrinkage = e$shrinkage,
       flo = f$lower, fhi = f$upper)
})
pc("expcnts", "rmoriebricklayer", function() {
  d <- pin$dsr
  bl("expcnts")(d$x, d$n, d$grp, d$age)$expected
})
pc("rtcalc", "rmoriebricklayer", function() {
  x <- pin$chg
  r <- bl("rtcalc")(x, count, pop, per = 1000)
  k <- bl("ratechg")(x, count, pop, year)
  list(rate = r$rate, lower = r$lower, upper = r$upper, rr = k$rateratio[2],
       pctlo = k$pctlo[2], pctup = k$pctup[2])
})
pc("riskratio", "rmorie", function() {
  t <- pin$tab
  a <- rm("risk_ratio")(t[["a"]], t[["b"]], t[["c"]], t[["d"]])
  b <- rm("morie_risk_ratio_ci")(matrix(t, 2, byrow = TRUE))
  d <- rm("morie_risk_difference_ci")(matrix(t, 2, byrow = TRUE))
  w <- rm("risk_difference")(t[["a"]], t[["b"]], t[["c"]], t[["d"]])
  z <- rm("risk_ratio")(0, 10, 5, 5)
  list(rr = a$estimate, lo = a$ci_lower, hi = a$ci_upper, rr2 = b$rr,
       lo2 = b$ci_lower, hi2 = b$ci_upper, rd = d$rd, rdlo = d$ci_lower,
       rdhi = d$ci_upper, wrd = w$estimate, wlo = w$ci_lower, whi = w$ci_upper,
       zrr = z$estimate, zlo = z$ci_lower)
})
pc("uofdisp", "rmorie", function() {
  r <- rm("mrm_uof_demographic_disparity")(pin$uof, "race", "force")
  p <- r$per_category
  list(category = vapply(p, function(e) e$category, ""),
       rate = vapply(p, function(e) e$rate, 0),
       lo = vapply(p, function(e) e$lo, 0), hi = vapply(p, function(e) e$hi, 0),
       rr = vapply(p, function(e) e$rr, 0), baseline = r$baseline)
})
pc("siusvc", "rmorie", function() {
  r <- rm("morie_siu_by_police_service")(pin$siu)
  list(counts = unlist(r$payload$counts), charged = r$payload$charged,
       rows = r$tables[[1]]$rows)
})
pc("tpsuof", "rmorie", function() {
  r <- rm("morie_tps_use_of_force")(c("a", "b", "a", NA, "c"), 250)
  list(rate = r$rate, n = r$n, counts = unlist(r$type_counts))
})
