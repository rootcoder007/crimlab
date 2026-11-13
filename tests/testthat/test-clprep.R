ns <- asNamespace("crimlab")

# clprep(): the preset table, the recode language, filters, missing rows,
# recipes and refusals; then generic(clprep(data, preset)) against the
# rmorie analyses it generalises (helper-port.R: live when rmorie is
# installed, recorded values otherwise).

test_that("the preset list and its table", {
  p <- clprep()
  expect_named(p, c("preset", "generic", "args", "source"))
  expect_false(anyDuplicated(p$preset) > 0)
  expect_true(all(p$generic %in% c(
    "cntreg", "grprr", "ctabx", "fairdf", "tmpaud", "calaud", "scrdsp",
    "cbatt", "lodimp", "tpbhm", "stkflw", "bndxw"
  )))
  # every argument a preset fills is an argument of its generic
  tab <- ns$.presets()
  for (k in unique(tab$preset)) {
    r <- tab[tab$preset == k, ]
    a <- setdiff(r$arg, "recipe")
    expect_true(all(a %in% names(formals(r$generic[1]))), info = k)
  }
  # the two SIUs are kept apart
  expect_true(all(grepl(
    "Special Investigations Unit",
    p$source[startsWith(p$preset, "siu-")],
    fixed = TRUE
  )))
  expect_true(all(grepl(
    "Structured Intervention Units",
    p$source[startsWith(p$preset, "sdsiu-")],
    fixed = TRUE
  )))
})

test_that("recodes, filters and dropped rows", {
  rc <- ns$.recode
  x <- c(" Female", "male", "FEMALE", NA)
  expect_identical(rc(x, "female=1;*=0", 4), c(1, 0, 1, NA))
  expect_identical(rc(c("Age 50+", "25 to 49"), "has:50=1;*=0", 2), c(1, 0))
  expect_identical(
    rc(c("yes", "no", "maybe"), "yes=1;no=0;*=NA", 3),
    c(1, 0, NA)
  )
  expect_identical(rc(c("a", "b"), "a=x;b=y", 2), c("x", "y"))
  expect_error(rc(c("a", "c"), "a=1", 2), "has no recode: c")
  expect_identical(rc("2021-04-01", "sub:1,4", 1), "2021")
  expect_identical(rc(1:2, "fac", 2), factor(1:2))
  expect_identical(rc(c("1", "2"), "num", 2), c(1, 2))
  expect_identical(rc(NULL, "one", 3), c(1, 1, 1))
  expect_true(rc(NULL, "const:TRUE", 3))
  expect_identical(rc(NULL, "const:prop", 3), "prop")
  d <- data.frame(
    EndFiscalYear = c(2020, 2020, 2021, 2021, 2022),
    Gender = c("Female", "Male", "Female", NA, "Male"),
    Measure = "Median",
    Region_AtTimeOfPlacement = "Central",
    NumberConsecutiveDays_Segregation = c(4, 7, 5, 6, 8)
  )
  d$Measure[5] <- "Max"
  p <- clprep(d, "otis-b04-cnt")
  expect_s3_class(p, "clpre")
  expect_identical(attr(p, "drop"), 2L)
  expect_identical(p$y, c(4, 7, 5))
  expect_identical(p$grp, c(1, 0, 1))
  expect_named(p$adj, c("Region_AtTimeOfPlacement", "EndFiscalYear"))
  expect_identical(attr(p, "generic"), "cntreg")
  expect_output(print(p), "otis-b04-cnt for cntreg")
})

test_that("recipes: long, and the person-year reductions", {
  d <- data.frame(
    Alert_Type = c("a", "b"),
    Gender = "F",
    EndFiscalYear = 2020,
    Number_Segregation_Placements_With_Alert = c(3, 4),
    Number_Segregation_Placements_Without_Alert = c(9, 8)
  )
  p <- clprep(d, "otis-b07-cnt")
  expect_identical(p$y, c(3, 4, 9, 8))
  expect_identical(p$grp, c(1L, 1L, 0L, 0L))
  b <- ob01
  v <- ns$.acvm(b)
  i <- b$UniqueIndividual_ID == v$UniqueIndividual_ID[1] &
    b$EndFiscalYear == v$EndFiscalYear[1]
  combo <- paste(
    b$MentalHealth_Alert,
    b$SuicideRisk_Alert,
    b$SuicideWatch_Alert
  )[i]
  expect_identical(v$ac[1], length(unique(combo)))
  # alert flags given as logical or numbers read the same
  bl <- b
  bl$MentalHealth_Alert <- b$MentalHealth_Alert == "Yes"
  bl$SuicideRisk_Alert <- as.numeric(b$SuicideRisk_Alert == "Yes")
  expect_identical(ns$.acvm(bl), v)
  y <- ns$.pyr(b)
  expect_identical(sum(y$np), nrow(b))
  expect_error(ns$.recipe("nope", "", b), "unknown recipe nope")
  expect_error(ns$.pyr(b[1:3]), "has no column")
})

test_that("refusals", {
  errs(
    "clprep",
    list(
      list("`preset` must name a preset", data.frame(a = 1)),
      list("`preset` must be one name", data.frame(a = 1), c("a", "b")),
      list("no preset named zz", data.frame(a = 1), "zz"),
      list("`data` must be a data frame", 1:3, "otis-c01-cnt"),
      list("has no column", data.frame(a = 1), "otis-c01-cnt")
    )
  )
})

# ---- parity: generic(clprep(data, preset)) against rmorie

test_that("cntreg(clprep()) reproduces the OTIS aggregate count models", {
  b <- portref("otiscnt")
  for (k in setdiff(names(b), "vm")) {
    src <- sub("-RC", "", k, fixed = TRUE)
    pre <- sprintf(
      "otis-%s%s-cnt",
      src,
      if (grepl("RC", k, fixed = TRUE)) "rc" else ""
    )
    # year clusters with year effects leave the GEE singular (rmorie
    # reports its covariance as degenerate there too)
    if (k == "c01-RC") {
      expect_warning(r <- cntreg(clprep(ocnt[[src]], pre)), "GEE fit failed")
    } else {
      r <- cntreg(clprep(ocnt[[src]], pre))
    }
    co <- r$coef[r$coef$term == "1", ]
    got <- cbind(
      c(co$est[co$fit == "pois"], co$se[co$fit == "pois"], r$fit$aic[1]),
      c(co$est[co$fit == "nb"], co$se[co$fit == "nb"], r$fit$aic[2])
    )
    expect_equal(got, unname(b[[k]]), tolerance = 1e-6, info = k)
  }
  # the alert-complexity count: rmorie's Poisson row, as it rounds it;
  # its negative binomial fit fails on these counts (no finite theta)
  # where crimlab reports theta at MASS's iteration limit
  v <- cntreg(clprep(ob01, "otis-vm-cnt"))
  pr <- v$irr[v$irr$fit == "pois", ]
  expect_identical(as.character(round(pr$irr, 4)), b$vm[[1]][2])
  expect_identical(sprintf("[%.3f, %.3f]", pr$lo, pr$hi), b$vm[[1]][3])
  expect_identical(b$vm[[2]][2], "fit failed")
  expect_false(v$fit$conv[2])
})

test_that("ctabx(clprep()) reproduces the OTIS tables", {
  b <- portref("otisxtab")
  r <- ctabx(clprep(pin$c04, "otis-c04-xtab"))
  row <- b$c04$rows[[1]]
  expect_identical(as.character(r$tab[row[1], b$c04$head[2:6]]), row[2:6])
  r <- ctabx(clprep(pin$b04, "otis-b04-xtab"))
  row <- b$b04$rows[[1]]
  expect_identical(as.character(r$tab[row[1], b$b04$head[2:4]]), row[2:4])
  r <- ctabx(clprep(pin$d02, "otis-dyr-xtab"))
  expect_identical(
    as.character(r$tab[, 1]),
    vapply(b$d02[[1]], function(e) e[2], "")
  )
  ch <- portref("cchi2")
  s <- ctabx(clprep(pin$cser, "otis-c05-xtab"))$stat
  expect_equal(c(s$chisq, s$df, s$p, s$v), unname(ch$c05), tolerance = 1e-12)
  l <- portref("locality")
  s <- ctabx(clprep(pin$loc, "otis-loc-xtab"))$stat
  expect_equal(c(s$chisq, s$p, s$v, s$diag), l$uof5, tolerance = 1e-12)
  m <- portref("mort")
  s <- ctabx(clprep(pin$mort, "otis-mh-sr-xtab"))$stat
  expect_identical(
    c(round(s$chisq, 2), round(s$v, 4)),
    c(m$chi2[1], m$cramers_v[1])
  )
})

test_that("person-level OTIS descriptives through the recipes", {
  b <- portref("otisprim")
  ac <- ctabx(clprep(ob01, "otis-ac-xtab"))$tab
  expect_equal(
    as.numeric(ac[as.character(b$ac$ac), 1]),
    b$ac$n_persons,
    tolerance = 0
  )
  vm <- ctabx(clprep(ob01, "otis-nreg-xtab"))$tab
  expect_equal(as.numeric(vm[names(b$vm), 1]), as.numeric(b$vm), tolerance = 0)
  np <- clprep(ob01, "otis-np-xtab")
  expect_equal(
    c(length(np$row), mean(np$row), max(np$row), mean(np$row > 1)),
    unname(b$np),
    tolerance = 1e-14
  )
  d21 <- ob01[ob01$EndFiscalYear == 2021, ]
  rp <- ctabx(clprep(d21, "otis-rplace-xtab"))$tab
  expect_equal(
    rp[rownames(b$rplace), intersect(colnames(b$rplace), colnames(rp))],
    b$rplace[, intersect(colnames(b$rplace), colnames(rp))],
    tolerance = 0,
    ignore_attr = TRUE
  )
  tr <- ctabx(clprep(ob01, "otis-rctrnd-xtab"))$tab
  expect_equal(
    tr[cbind(as.character(b$rctrnd$year), b$rctrnd$region)],
    b$rctrnd$n_individuals,
    tolerance = 0
  )
  ny <- ctabx(clprep(ob01, "otis-nyear-xtab"))$tab
  expect_equal(
    as.numeric(ny[, 1]),
    vapply(b$nyear, function(e) as.numeric(e[[2]]), 0),
    tolerance = 0
  )
})

test_that("clprep() reproduces the SIU and ARSAU analyses", {
  s <- portref("siusvc")
  g <- grprr(clprep(pin$siu, "siu-svc-rr"))
  expect_equal(
    g$rates$n[match(names(s$counts), g$rates$group)],
    as.numeric(s$counts),
    tolerance = 0
  )
  q <- portref("siuchi")
  p <- clprep(pin$siuy, "siu-yrchg-xtab")
  t <- ctabx(p)$tab
  keep <- rownames(t)[apply(t > 0, 1, all)]
  i <- p$row %in% keep
  r <- ctabx(p$row[i], p$col[i], correct = TRUE)$stat
  expect_equal(c(r$chisq, r$df, r$p), unname(q$stat), tolerance = 1e-12)
  y <- ctabx(clprep(pin$siuy, "siu-year-xtab"))$tab
  expect_equal(as.numeric(y[, 1]), as.numeric(q$years), tolerance = 0)
  o <- ctabx(clprep(pin$siuo, "siu-outcome-xtab"))$tab
  expect_equal(
    o[cbind(as.character(q$out$service), as.character(q$out$outcome))],
    q$out$n_cases,
    tolerance = 0
  )
  u <- portref("uofdisp")
  a <- data.frame(
    Race = pin$uof$race,
    IndivInjuries_PhysicalInjuries = ifelse(pin$uof$force == 1, "Yes", "No")
  )
  a$IndivInjuries_PhysicalInjuries[1] <- "unknown"
  pr <- clprep(a, "arsau-indiv-race-rr")
  expect_identical(attr(pr, "drop"), 1L)
  w <- portref("weapon")
  aw <- data.frame(Weapon = pin$wpn$weapon, Location = pin$wpn$force)
  r <- ctabx(clprep(aw, "arsau-weapon-xtab"))$stat
  expect_equal(c(r$chisq, r$df, r$p, r$v), w$stat, tolerance = 1e-12)
  full <- grprr(
    clprep(
      data.frame(
        Race = pin$uof$race,
        IndivInjuries_PhysicalInjuries = pin$uof$force
      ),
      "arsau-indiv-race-rr"
    ),
    ref = u$baseline
  )
  expect_equal(
    full$rates$rate[match(u$category, full$rates$group)],
    u$rate,
    tolerance = 1e-14
  )
})

test_that("presets on the Structured Intervention Units and TPS layouts", {
  d <- data.frame(
    region = c("Atlantic", "Ontario", "Quebec"),
    siu_pct = c(10, 8.5, 37.3),
    pop_pct = c(9, 27.5, 19.8)
  )
  g <- grprr(clprep(d, "sdsiu-t12-rr"))
  expect_equal(
    g$rates$rep,
    d$siu_pct / d$pop_pct * sum(d$pop_pct) / sum(d$siu_pct),
    tolerance = 1e-14
  )
  r <- data.frame(region = c("Ontario", "Quebec"), short_stay_rate = c(4, 38))
  expect_equal(
    grprr(clprep(r, "sdsiu-t13-rr"))$cmp$ratio,
    9.5,
    tolerance = 1e-15
  )
  t <- data.frame(DIVISION = c("D11", "D14", "D11"))
  expect_identical(
    as.numeric(ctabx(clprep(t, "tps-division-xtab"))$tab[, 1]),
    c(2, 1)
  )
})

test_that("the risk-tool presets feed fairdf, tmpaud, calaud and scrdsp", {
  d <- fairsim(400, seed = 3)
  # a preset gives the same result as the columns given by hand
  expect_identical(
    calaud(clprep(d, "fairsim-cal"), nperm = 99, seed = 1),
    calaud(d$area, d$risk, d$det, d$group, nperm = 99, seed = 1)
  )
  expect_identical(
    fairdf(clprep(d, "fairsim-df")),
    fairdf(d$det, d$group, base = "A")
  )
  # an argument given in the call wins over the preset's
  expect_identical(
    fairdf(clprep(d, "fairsim-df"), base = "B"),
    fairdf(d$det, d$group, base = "B")
  )
  expect_identical(
    scrdsp(clprep(d, "predpol-generic-score")),
    scrdsp(d$risk, d$group)
  )
  d$pred <- d$det
  d$period <- rep(1:2, 200)
  d$city <- rep(c("north", "south"), each = 200)
  expect_identical(
    tmpaud(clprep(d, "predpol-generic-time")),
    tmpaud(d$pred, d$group, d$period, d$city)
  )
  d$outcome <- d$det
  d$population <- 1000
  expect_identical(
    calaud(clprep(d, "predpol-generic"), nperm = 99, seed = 2),
    calaud(d$area, d$risk, d$outcome, d$group, d$population,
           nperm = 99, seed = 2)
  )
  # the two cut-offs of a three-level score: low/medium against high, and
  # low against medium/high
  l <- data.frame(
    race = rep(c("a", "b"), each = 6),
    static_score = rep(c("low", "medium", "high"), 4)
  )
  mh <- clprep(l, "laniyonu-ard-static-mh")
  lm <- clprep(l, "laniyonu-ard-static-lm")
  expect_identical(as.numeric(mh$pred), rep(c(0, 0, 1), 4))
  expect_identical(as.numeric(lm$pred), rep(c(0, 1, 1), 4))
  expect_identical(fairdf(mh), fairdf(rep(c(0, 0, 1), 4), l$race))
  # every risk-tool preset names its source
  pr <- clprep()
  tool <- pr[pr$generic %in% c("fairdf", "tmpaud", "calaud", "scrdsp"), ]
  expect_gt(nrow(tool), 10)
  expect_true(all(nzchar(tool$source)))
})

test_that("a preset of one's own, and the flag and dum recodes", {
  mine <- data.frame(
    preset = "my-cnt", generic = "cntreg", arg = c("y", "grp", "adj"),
    column = c("count", "sex", "year"), recode = c("", "f=1;*=0", "fac")
  )
  e <- data.frame(count = c(4, 30, 6, 41, NA), sex = c("f", "m", "f", "m", "f"),
                  year = c(1, 1, 2, 2, 2))
  p <- clprep(e, mine)
  expect_identical(p$y, c(4, 30, 6, 41))
  expect_identical(p$grp, c(1, 0, 1, 0))
  expect_identical(p$adj, data.frame(year = factor(c(1, 1, 2, 2))))
  expect_identical(attr(p, "preset"), "my-cnt")
  expect_identical(attr(p, "drop"), 1L)
  expect_equal(cntreg(p), cntreg(p$y, p$grp, adj = p$adj))
  expect_error(clprep(e, mine[, 1:2]), "needs the columns preset, generic")
  expect_error(clprep(e, mine[0, ]), "needs the columns preset, generic")
  two <- mine
  two$preset[3] <- "other"
  expect_error(clprep(e, two), "must hold one preset")
  rc <- ns$.recode
  expect_identical(rc(c("Yes", "no", "2", "0"), "flag", 4), c(1L, 0L, 0L, 0L))
  expect_identical(rc(c(2, 0), "flag", 2), c(1L, 0L))
  m <- rc(c("b", "a", "c", "b"), "dum", 4)
  expect_identical(m, cbind(b = c(1, 0, 0, 1), c = c(0, 0, 1, 0)))
  expect_identical(dim(rc(c("a", "a"), "dum", 2)), c(2L, 0L))
})

test_that("the OTIS treatment-outcome pairs of cbatt", {
  b <- ob01
  set.seed(9)
  b$Number_Of_Placements <- rpois(nrow(b), 1.5)
  b$NumberConsecutiveDays_Segregation <- rexp(nrow(b), 0.1)
  b$Gender[3] <- NA
  dum <- function(x) {
    f <- factor(x)
    m <- vapply(levels(f)[-1L], function(l) as.numeric(f == l),
                numeric(length(f)))
    matrix(m, length(f))
  }
  flag <- function(s) {
    as.integer(tolower(s) %in% c("yes", "y", "true", "t", "1"))
  }
  # a01: Ruhela's reduction (acvm) as the cntreg preset uses it
  v <- ns$.acvm(b)
  p <- clprep(b, "otis-a01-cbatt")
  expect_identical(p$y, v$vm)
  expect_identical(p$d, v$hiac)
  expect_identical(names(p$x)[1], "GenderMale")
  x <- cbind(dum(v$Gender), dum(v$Age_Category), dum(v$EndFiscalYear))
  expect_identical(unname(as.matrix(p$x)), x)
  expect_equal(cbatt(p, reps = 0), cbatt(v$vm, v$hiac, x, reps = 0))
  # a01n: the naive reduction, from scratch
  cc <- b[stats::complete.cases(b[c(ns$.otcols, "Gender", "Age_Category")]), ]
  key <- paste(cc$UniqueIndividual_ID, cc$EndFiscalYear)
  nal <- flag(cc$MentalHealth_Alert) + flag(cc$SuicideRisk_Alert) +
    flag(cc$SuicideWatch_Alert)
  ch <- as.integer(cc$Region_AtTimeOfPlacement != cc$Region_MostRecentPlacement)
  n <- clprep(b, "otis-a01n-cbatt")
  k <- unique(key[order(cc$UniqueIndividual_ID, cc$EndFiscalYear)])
  expect_identical(n$d, as.integer(c(tapply(nal, key, max))[k] >= 2L))
  expect_identical(n$y, unname(c(tapply(ch, key, max))[k]))
  # pair a: row-level flags
  a <- clprep(b, "otis-pa-cbatt")
  ka <- !is.na(b$Gender)
  expect_identical(a$y, flag(b$SuicideRisk_Alert[ka]))
  expect_identical(a$d, flag(b$MentalHealth_Alert[ka]))
  nl <- vapply(b[ka, c("Gender", "Age_Category", "Region_AtTimeOfPlacement",
                       "Region_MostRecentPlacement")],
               function(z) length(unique(z)) - 1L, 1L)
  expect_identical(ncol(a$x), sum(nl))
  # pair b: two or more alerts, two or more placements
  pb <- clprep(b, "otis-pb-cbatt")
  bb <- b[ka, ]
  na3 <- flag(bb$MentalHealth_Alert) + flag(bb$SuicideRisk_Alert) +
    flag(bb$SuicideWatch_Alert)
  expect_identical(pb$d, as.integer(na3 >= 2L))
  expect_identical(pb$y, as.integer(bb$Number_Of_Placements >= 2))
  # pair c: region change, days capped at the 99th percentile
  pc <- clprep(b, "otis-pc-cbatt")
  y <- bb$NumberConsecutiveDays_Segregation
  expect_identical(pc$y, pmin(y, quantile(y, 0.99, names = FALSE)))
  expect_identical(
    pc$d,
    as.integer(bb$Region_AtTimeOfPlacement != bb$Region_MostRecentPlacement)
  )
  expect_identical(names(pc$x)[ncol(pc$x)], "MentalHealth_Alert")
  expect_error(clprep(b[1:5], "otis-pc-cbatt"), "has no column")
})

test_that("the gentrification recipe", {
  set.seed(10)
  d <- data.frame(inc = rnorm(50, 50), rent = rnorm(50, 900, 50),
                  colg = runif(50), rentg = runif(50))
  d$inc[4] <- NA
  g <- ns$.recipe("gent", "inc,rent,colg,rentg", d)$gent
  q <- function(x, p) quantile(x, p, na.rm = TRUE, names = FALSE)
  el <- d$inc <= q(d$inc, 0.5) & d$rent <= q(d$rent, 0.5)
  gr <- d$colg >= q(d$colg, 0.667) & d$rentg >= q(d$rentg, 0.667)
  ref <- ifelse(el, ifelse(gr, "gentrified", "eligible"), "ineligible")
  expect_identical(as.character(g), ref)
  expect_identical(levels(g), c("ineligible", "eligible", "gentrified"))
  g2 <- ns$.recipe("gent", "inc,rent,colg,rentg,0.3,0.9", d)$gent
  el <- d$inc <= q(d$inc, 0.3) & d$rent <= q(d$rent, 0.3)
  gr <- d$colg >= q(d$colg, 0.9) & d$rentg >= q(d$rentg, 0.9)
  expect_identical(
    as.character(g2),
    ifelse(el, ifelse(gr, "gentrified", "eligible"), "ineligible")
  )
  expect_error(ns$.recipe("gent", "inc,rent,colg,nope", d), "no column nope")
  # through a preset of one's own
  mine <- data.frame(preset = "gent-xtab", generic = "ctabx",
                     arg = c("recipe", "row", "col"),
                     column = c("gent", "gent", ""),
                     recode = c("inc,rent,colg,rentg", "", "one"))
  p <- clprep(d, mine)
  expect_identical(p$row, g[!is.na(g)])
})

test_that("the toxicology and taphonomy layouts", {
  d <- data.frame(case_id = c("a", "b", "c", "d"), conc = c(0.4, 0.05, NA, 1),
                  lod = c(0.1, 0.1, 0.2, 0.5))
  p <- clprep(d, "tox-lodimp")
  expect_identical(attr(p, "drop"), 1L)
  expect_identical(lodimp(p), lodimp(d$conc[-3], d$lod[-3]))
  set.seed(11)
  cols <- c("lime_treatment", "temp_c", "humidity_pct", "arid",
            "casket_sealed", "reinterment_count", "exposure_days",
            "decades_elapsed")
  b <- as.data.frame(matrix(runif(40 * 8), 40, 8, dimnames = list(NULL, cols)))
  b$preservation_score <- runif(40)
  t1 <- tpbhm(clprep(b, "taph-tpbhm"))
  expect_equal(t1, tpbhm(b$preservation_score, as.matrix(b[cols])))
  pc <- c("decomp_stage", "body_scoring_tbs", "accumulated_deg_days",
          "temp_c", "humidity_pct", "precipitation_mm", "burial_depth_cm",
          "soil_ph", "scavenger_activity", "insect_activity")
  e <- as.data.frame(matrix(rnorm(60 * 10), 60, 10, dimnames = list(NULL, pc)))
  e$pmi_days <- rexp(60, 0.05)
  expect_equal(tpbhm(clprep(e, "taph-pmi-tpbhm"), psd = 5),
               tpbhm(e$pmi_days, as.matrix(e[pc]), psd = 5))
})

test_that("the Chicago and New York count tables", {
  set.seed(12)
  n <- 300
  cr <- data.frame(
    primary_type = sample(c("THEFT", "BATTERY", "ASSAULT"), n, TRUE),
    year = sample(2019:2021, n, TRUE),
    community_area = sample(1:5, n, TRUE),
    arrest = sample(c(TRUE, FALSE), n, TRUE, prob = c(0.2, 0.8))
  )
  tb <- function(p, d) as.numeric(ctabx(clprep(d, p))$tab)
  expect_identical(tb("chicago-type-xtab", cr),
                   as.numeric(table(cr$primary_type)))
  expect_identical(tb("chicago-year-xtab", cr),
                   as.numeric(table(cr$year)))
  ar <- ifelse(cr$arrest, "arrest", "no arrest")
  expect_identical(tb("chicago-area-xtab", cr),
                   as.numeric(table(cr$community_area, ar)))
  ca <- data.frame(race = sample(c("W", "B", "H"), 80, TRUE))
  expect_identical(tb("chicago-race-xtab", ca),
                   as.numeric(table(ca$race)))
  ny <- data.frame(
    ofns_desc = sample(c("ASSAULT 3", "PETIT LARCENY"), n, TRUE),
    arrest_boro = sample(c("B", "K", "M", "Q", "S"), n, TRUE),
    perp_race = sample(c("BLACK", "WHITE", "ASIAN"), n, TRUE),
    law_cat_cd = sample(c("F", "M", "V"), n, TRUE)
  )
  expect_identical(tb("nypd-offense-xtab", ny),
                   as.numeric(table(ny$ofns_desc)))
  expect_identical(tb("nypd-boro-xtab", ny),
                   as.numeric(table(ny$arrest_boro)))
  fel <- ifelse(ny$law_cat_cd == "F", "felony", "other")
  expect_identical(tb("nypd-felony-xtab", ny),
                   as.numeric(table(ny$perp_race, fel)))
  cp <- data.frame(susp_race = sample(c("BLACK", "WHITE", NA), n, TRUE))
  p <- clprep(cp, "nypd-susprace-xtab")
  expect_identical(attr(p, "drop"), sum(is.na(cp$susp_race)))
  expect_identical(tb("nypd-susprace-xtab", cp),
                   as.numeric(table(cp$susp_race)))
})

test_that("the Doob affidavit tables", {
  # the published tables (Corrections and Conditional Release Statistical
  # Overview), as rmorie transcribed them
  t1 <- data.frame(
    type = c("Day Parole", "Full Parole", "Statutory Release"),
    revoke_violent = c(4.8, 5.0, 84.6),
    success = c(3085, 901.2, 3735.6),
    total = c(3392.2, 1023.8, 5829.0)
  )
  a <- grprr(clprep(t1, "doob-t1-rr"))
  expect_equal(a, grprr(t1$success, t1$total, t1$type, type = "prop",
                        per = 100))
  expect_equal(a$rates$rate, 100 * t1$success / t1$total, tolerance = 1e-14)
  v <- grprr(clprep(t1, "doob-t1v-rr"))
  expect_equal(v$rates$rate, 100 * t1$revoke_violent / t1$total,
               tolerance = 1e-14)
  t2 <- data.frame(
    year = c("2013-14", "2014-15", "2015-16", "2016-17", "2017-18"),
    avg_count = c(15342, 14886, 14712, 14159, 14092),
    admissions = c(5071, 4818, 4891, 4908, 4718),
    statutory_releases = c(5636, 5373, 5309, 4888, 4427)
  )
  f <- stkflw(clprep(t2, "doob-t2-flow"))
  expect_equal(f, stkflw(people = t2$admissions, adp = t2$avg_count,
                         prd = t2$year, rel = t2$statutory_releases,
                         early = 182.5))
  expect_equal(f$flow$cut, t2$statutory_releases * 182.5 / 365,
               tolerance = 1e-14)
  expect_equal(f$flow$alos, t2$avg_count * 365 / t2$admissions,
               tolerance = 1e-14)
  t3 <- data.frame(
    age_group = c("18-49", "50-59", "60+"),
    canada_adult_pop = c(15770626, 5305888, 8811720),
    csc_in_custody = c(10544, 2236, 1312),
    csc_admissions_2017_18 = c(3920, 548, 250)
  )
  r <- grprr(clprep(t3, "doob-t3-rr"))
  sh <- t3$csc_in_custody / sum(t3$csc_in_custody)
  sp <- t3$canada_adult_pop / sum(t3$canada_adult_pop)
  expect_equal(r$rates$rep, sh / sp, tolerance = 1e-14)
  expect_equal(r$rates$rate, 1e5 * t3$csc_in_custody / t3$canada_adult_pop,
               tolerance = 1e-14)
  # the affidavit's over-representation of the young, from rounded shares
  expect_equal(r$rates$rep[1], 74.8 / 52.8, tolerance = 2e-3)
  ad <- grprr(clprep(t3, "doob-t3a-rr"))
  expect_equal(ad$rates$x, t3$csc_admissions_2017_18)
})

test_that("the Victorian tables, the pop recipe and numeric constants", {
  set.seed(15)
  lga <- data.frame(
    Year = 2024,
    `Local Government Area` = c("Ararat", "Ballarat", "Casey"),
    `Incidents Recorded` = c(1200, 9800, 25000),
    `Rate per 100,000 population` = c(9800.5, 8650.2, 6900.1),
    check.names = FALSE
  )
  p <- clprep(lga, "vic-lga-rr")
  expect_identical(p$per, 1e5)
  expect_equal(p$n, lga[[3]] / lga[[4]] * 1e5, tolerance = 1e-15)
  g <- grprr(p)
  expect_equal(g$rates$rate, lga[[4]], tolerance = 1e-12)
  ind <- data.frame(
    Year = 2024,
    `Offence Division` = rep(c("A Crimes against the person",
                               "B Property and deception offences"),
                             each = 3),
    `Indigenous Status` = rep(c("Aboriginal and/or Torres Strait Islander",
                                "Non-Indigenous", "Total People"), 2),
    `Victim Reports` = c(120, 3000, 3120, 40, 2500, 2540),
    check.names = FALSE
  )
  x <- ctabx(clprep(ind, "vic-indig-xtab"))
  expect_identical(as.numeric(x$tab), c(120, 40, 3000, 2500))
  expect_identical(colnames(x$tab), c("aboriginal", "non-indigenous"))
  tr <- data.frame(Year = rep(2020:2022, each = 2),
                   `Offence Division` = rep(c("A", "B"), 3),
                   `Incidents Recorded` = c(5, 7, 6, 9, 8, 4),
                   check.names = FALSE)
  y <- ctabx(clprep(tr, "vic-trend-xtab"))
  expect_identical(as.numeric(y$tab), c(5, 6, 8, 7, 9, 4))
  rc <- ns$.recode
  expect_identical(rc(NULL, "const:2.5", 1), 2.5)
  expect_identical(rc(c("50", "100"), "pct", 2), c(0.5, 1))
  expect_error(ns$.recipe("pop", "a;b;100", data.frame(a = 1)),
               "no column b")
})

test_that("the Toronto neighbourhood crosswalks", {
  cw <- data.frame(
    hood_140 = c("001", "002", "002", "003"),
    hood_158 = c("001", "160", "161", "003"),
    pct_140_in_158 = c(100, 40, 60, 100),
    pct_158_in_140 = c(100, 100, 100, 100)
  )
  up <- bndxw(clprep(cw, "tps-hood158-140"), x = c(5, 2, 3, 7),
              unit = c("001", "160", "161", "003"))
  expect_identical(up$tab$value, c(5, 5, 7))
  dn <- bndxw(clprep(cw, "tps-hood140-158"), x = c(5, 10, 7),
              unit = c("001", "002", "003"))
  expect_identical(dn$tab$to, c("001", "003", "160", "161"))
  expect_equal(dn$tab$value, c(5, 7, 4, 6), tolerance = 1e-15)
  expect_identical(dn$map$to, c("001", "161", "003"))
})
