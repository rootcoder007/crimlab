ns <- asNamespace("crimlab")

# clprep(): the preset table, the recode language, filters, missing rows,
# recipes and refusals; then generic(clprep(data, preset)) against the
# rmorie analyses it generalises (helper-port.R: live when rmorie is
# installed, recorded values otherwise).

test_that("the preset list and its table", {
  p <- clprep()
  expect_named(p, c("preset", "generic", "args", "source"))
  expect_false(anyDuplicated(p$preset) > 0)
  expect_true(all(p$generic %in% c("cntreg", "grprr", "ctabx")))
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
