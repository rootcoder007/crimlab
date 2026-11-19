# The footprint of a computation: the values of the reference
# implementation (rmoriebricklayer 1.1.1, computed once on these inputs and
# printed with 17 significant digits), the energy of runs with fixed times
# through the kernel, the devices table, and the refusals.

src_ca <- paste(
  "Electricity Maps (2025). Canada 2024 Yearly Carbon Intensity Data",
  "(Version January 27, 2025). Electricity Maps (via Green Algorithms data",
  "v3.1)"
)
src_ember <- paste(
  "Ember (2026) Yearly Electricity Data via Our World in Data, lifecycle",
  "gCO2e/kWh, CC BY"
)
ci_names <- c("g_per_kwh", "location", "name", "year", "basis", "source",
              "note")

test_that("co2int() gives the reference rows", {
  r <- co2int("CA-ON")
  expect_s3_class(r, "co2int")
  expect_named(r, ci_names)
  # the intensity is a weighted sum: its last digit differs between systems
  expect_equal(
    unclass(r),
    list(g_per_kwh = 90.969999999999999, location = "CA-ON",
         name = "Canada, Ontario", year = 2024L, basis = "lifecycle",
         source = src_ca, note = character(0)),
    tolerance = 1e-14
  )
  for (k in c("Canada", "can")) {
    expect_identical(co2int(k)[1:6], co2int("CA")[1:6])
  }
  expect_identical(co2int("Canada")$g_per_kwh, 190.72)
  expect_identical(co2int("Canada")$source, src_ember)
  expect_identical(co2int("CA-XX")$note,
                   "no row for zone CA-XX; its country CA used")
  expect_identical(co2int("CA-XX")$location, "CA")
  u <- co2int(" us-cal-ciso ")
  expect_identical(u$location, "US-CAL-CISO")
  expect_identical(u$g_per_kwh, 227.13999999999999)
  expect_identical(co2int()$g_per_kwh, 458.29000000000002)
  expect_identical(
    unclass(co2int(120)),
    list(g_per_kwh = 120, location = "user", name = "user-supplied",
         year = NA_integer_, basis = "user",
         source = "user-supplied carbon intensity", note = character(0))
  )
})

test_that("co2inttab() is the reference table", {
  tab <- co2inttab()
  expect_identical(dim(tab), c(378L, 7L))
  expect_named(
    tab,
    c("location", "iso3", "name", "g_per_kwh", "year", "basis", "source")
  )
  expect_identical(
    unname(vapply(tab, class, "")),
    c("character", "character", "character", "numeric", "integer",
      "character", "character")
  )
  expect_identical(sum(tab$g_per_kwh), 152916.19907999999)
  expect_identical(tab$location[c(1, 50, 378)], c("WORLD", "HR", "AU-NT"))
  expect_identical(tab$g_per_kwh[c(1, 50, 378)],
                   c(458.29000000000002, 158.5, 494.56999999999999))
  expect_identical(tab$year[c(1, 50, 378)], c(2025L, 2025L, 2024L))
  expect_true(all(tab$basis == "lifecycle"))
  expect_false(anyDuplicated(tab$location) > 0)
})

test_that("dtctloc() finds the reference locations", {
  ref <- list(
    list(list(envloc = "", codloc = "", tz = "Europe/Paris"),
         "FR", "system time zone Europe/Paris"),
    list(list(envloc = "", codloc = "", tz = "", locale = "en_IN.UTF-8"),
         "IN", "locale en_IN.UTF-8"),
    list(list(envloc = "CA-QC", tz = "Europe/Paris"),
         "CA-QC", "CRIMLAB_LOCATION environment variable"),
    list(list(envloc = "", codloc = "", tz = "UTC", locale = "C"),
         "WORLD",
         "fallback, no location found in environment, time zone or locale"),
    list(list(envloc = "", codloc = "DEU", tz = "UTC"),
         "DE", "CODECARBON_COUNTRY_ISO_CODE environment variable"),
    list(list(envloc = "", codloc = "", tz = "Europe/Monaco"),
         "FR", paste("system time zone Europe/Monaco (territory MC has no",
                     "row of its own; grid of FR)")),
    list(list(envloc = "", codloc = "", tz = "Antarctica/Troll"),
         "WORLD", paste("system time zone Antarctica/Troll (territory AQ",
                        "has no row of its own; world average)")),
    list(list(envloc = "", codloc = "", tz = "posix/Asia/Calcutta"),
         "IN", "system time zone posix/Asia/Calcutta"),
    list(list(envloc = "", codloc = "", tz = "Asia/Tokyo",
              locale = "fr_FR.UTF-8"),
         "JP", "system time zone Asia/Tokyo")
  )
  for (r in ref) {
    expect_identical(do.call(dtctloc, r[[1L]]),
                     list(location = r[[2L]], method = r[[3L]]))
  }
  expect_warning(
    d <- dtctloc(envloc = "Atlantis", codloc = "", tz = "Europe/Paris"),
    "CRIMLAB_LOCATION='Atlantis' is not a location"
  )
  expect_identical(d$location, "FR")
  expect_error(dtctloc(tz = c("a", "b")),
               "^dtctloc: `tz` must be one string \\(or NULL\\)\\.$")
})

test_that("fpeqv() divides as the reference does", {
  e <- fpeqv(1200)
  expect_s3_class(e, "fpeqv")
  expect_identical(e$carkm, 6.8571428571428568)
  expect_identical(e$treemon, 1.3086150490730644)
  expect_named(e$sources, c("car", "tree"))
  expect_identical(fpeqv(0)$carkm, 0)
  expect_error(fpeqv(c(10, 20)),
               "^fpeqv: `co2eg` must be a single non-negative number\\.$")
  expect_error(fpeqv(-1), "non-negative")
  expect_output(print(e), "Car travel")
})

test_that("raplok() reads a powercap tree", {
  skip_on_os("windows")
  d <- file.path(tempfile(), "intel-rapl:0")
  dir.create(d, recursive = TRUE)
  writeLines("package-0", file.path(d, "name"))
  writeLines("123456789", file.path(d, "energy_uj"))
  expect_true(raplok(dirname(d)))
  writeLines("not a number", file.path(d, "energy_uj"))
  expect_false(raplok(dirname(d)))
  expect_false(raplok(file.path(tempdir(), "no-such-dir")))
  expect_error(raplok(NA_character_),
               "^raplok: `rapldir` must be a single directory path\\.$")
})

# The source's energy for fixed times (its .fp_energy() on constructed
# proc.time() values): wall, CPU seconds, utilisation, then the CPU, memory
# and total kWh, PUE included.
fpe <- function(wall, cpu, cores, mode, u, pw, memgb, pue, curve = 0) {
  inp <- c(wall, cpu, cores, mode, u, NA, NA, NA, NA, pw, curve, memgb,
           0.3725, pue, 100, 0)
  .Call(C_clab_r_fpe, as.numeric(inp), numeric(0), numeric(0), numeric(0),
        numeric(0), numeric(0), integer(0))[[1L]]
}

test_that("the energy of a run equals the reference", {
  ref <- list(
    list(fpe(3600, 7200, 4, 0, NA, 100, 10, 1.5),
         c(3600, 0.5, 0.075000000000000011, 0.0055875000000000005,
           0.080587500000000006)),
    list(fpe(12.5, 3, 8, 2, 0.25, NA, 16, 1),
         c(12.5, 0.25, 0.00014756944444444445, 2.0694444444444442e-05,
           0.00016826388888888888)),
    list(fpe(0.002, 0.9, 2, 0, NA, 65, 4, 1.67, 1),
         c(0.45000000000000001, 1, 1.3568750000000001e-05, 3.110375e-07,
           1.3879787500000001e-05)),
    list(fpe(60, 30, 4, 2, 0.5, 45, 8, 1, 1),
         c(60, 0.5, 0.00037500000000000001, 4.9666666666666663e-05,
           0.00042466666666666667))
  )
  for (r in ref) {
    expect_equal(r[[1L]][c(1L, 2L, 6L, 7L, 9L)], r[[2L]], tolerance = 1e-12)
  }
  # the flags behind the reference's assumptions: the unknown CPU, then
  # codecarbon's curve
  expect_identical(ref[[2L]][[1L]][5L], 8)
  expect_identical(ref[[3L]][[1L]][5L], 16)
})

test_that("co2fp() has the reference's shape and assumptions", {
  fp <- co2fp(1 + 1, loc = "CA-ON", cpupwrw = 45, memgb = 8, usage = 0.5,
              mthd = "model")
  expect_s3_class(fp, "co2fp")
  expect_named(
    fp,
    c("value", "durs", "cputms", "usage", "usagemode", "cores", "memgb",
      "enekwh", "co2eg", "co2int", "method", "assumptions", "citations")
  )
  expect_named(fp$enekwh, c("cpu", "memory", "total"))
  expect_identical(fp$value, 2)
  expect_identical(fp$assumptions, character(0))
  expect_identical(fp$citations[3L], paste0(
    "carbon intensity: ", src_ca, " (Canada, Ontario, 2024)"
  ))
  h <- fp$durs / 3600
  expect_equal(fp$enekwh$memory, h * 8 * 0.3725 / 1000, tolerance = 1e-12)
  expect_equal(fp$enekwh$cpu, h * 45 * 0.5 / 1000, tolerance = 1e-12)
  # the product is formed in C++, where a fused multiply-add may round
  # the last bit otherwise than R does
  expect_equal(fp$co2eg, fp$enekwh$total * 90.969999999999999,
               tolerance = 1e-14)
  fp2 <- co2fp(1 + 1, loc = 30, cpupwrw = 65, memgb = 16, pue = 1.67,
               mthd = "model", loadcur = "codecarbon")
  expect_identical(
    fp2$assumptions,
    c("codecarbon process load curve applied to the CPU power",
      "PUE 1.67 applied")
  )
  expect_length(fp2$citations, 2L)
  expect_identical(fp2$method, "model")
  expect_output(print(fp2), "Computation footprint \\(modelled\\)")
})

test_that("co2fp() adds the energy of other devices", {
  dev <- data.frame(
    kind = c("gpu", "cpu", "gpu", "npu"),
    vendor = "example",
    name = c("card 0", "package 0", "card 1", "unit 0"),
    watts = c(120, 35, NA, NA),
    energy = c(NA, NA, 0.002, NA),
    how = "simulated reading",
    stringsAsFactors = FALSE
  )
  fp <- co2fp(1 + 1, loc = 30, cpupwrw = 45, memgb = 8, usage = 0.5,
              mthd = "model", devices = dev)
  h <- fp$durs / 3600
  expect_identical(fp$devices$counted, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(fp$devices$kwh, c(h * 120 / 1000, NA, 0.002, NA),
               tolerance = 1e-12)
  expect_identical(
    fp$devices$note,
    c("average power over the run", "already counted as CPU or memory energy",
      "measured energy", "no power or energy reading")
  )
  expect_equal(fp$enekwh$devices, h * 120 / 1000 + 0.002, tolerance = 1e-12)
  expect_equal(fp$enekwh$total,
               fp$enekwh$cpu + fp$enekwh$memory + fp$enekwh$devices,
               tolerance = 1e-12)
  expect_identical(fp$co2eg, fp$enekwh$total * 30)
  expect_false("devices" %in% names(co2fp(1, loc = 30, mthd = "model")))
  expect_error(co2fp(1, devices = dev[, c("kind", "name")]),
               "^co2fp: `devices` is missing columns: vendor, how\\.$")
  expect_error(co2fp(1, devices = dev[, c("kind", "vendor", "name", "how")]),
               "needs a `watts` or an `energy` column")
  bad <- dev
  bad$watts[1L] <- -1
  expect_error(co2fp(1, devices = bad), "`devices\\$watts` must be")
})

test_that("co2fp() refuses bad arguments before the expression runs", {
  ran <- FALSE
  expect_error(co2fp(ran <- TRUE, pue = 0.9),
               "^co2fp: `pue` must be a single number >= 1\\.$")
  expect_false(ran)
  expect_error(co2fp(1, loc = "Atlantis"), "^co2fp: Unknown location")
  expect_error(co2fp(1, cores = 1.5), "`cores` must be a whole number")
  expect_error(co2fp(1, cpupwrw = -1),
               "`cpupwrw` must be a single positive number")
  expect_error(co2fp(1, usage = 2), "`usage` must be")
  expect_error(co2fp(1, usage = "both"), "`usage` must be one of")
  expect_error(co2fp(1, mthd = "rapl", rapldir = tempfile()),
               "RAPL energy counters are not readable here")
  expect_error(co2int("Atlantis"), "^co2int: Unknown location 'Atlantis'")
  expect_error(co2int(c(1, 2)), "numeric `loc` must be a single")
})

test_that("co2int() prints its row and refuses a loc it cannot use", {
  expect_identical(
    capture.output(print(co2int("CA-XX"))),
    c("Carbon intensity of electricity: 190.72 gCO2e/kWh",
      "  Location: CA (Canada), 2025", "  Basis: lifecycle",
      paste0("  Source: ", src_ember),
      "  Note: no row for zone CA-XX; its country CA used")
  )
  # a user value has no year
  expect_identical(
    capture.output(print(co2int(120)))[2L],
    "  Location: user (user-supplied)"
  )
  for (bad in list(TRUE, c("a", "b"))) {
    expect_error(
      co2int(bad),
      paste("co2int: `loc` must be a single location code or name, a",
            "number in gCO2e/kWh, or NULL."),
      fixed = TRUE
    )
  }
})

test_that("a locale without a known country falls back to the world", {
  fb <- list(
    location = "WORLD",
    method = "fallback, no location found in environment, time zone or locale"
  )
  for (lc in c(NA_character_, "en_ZZ.UTF-8")) {
    expect_identical(dtctloc(envloc = "", codloc = "", tz = "", locale = lc),
                     fb)
  }
  # the binding's proxy for a territory: its neighbour's grid, and the
  # world average ("") for one the proxy list does not name
  expect_identical(.Call(C_clab_r_proxy, "MC"), "FR")
  expect_identical(.Call(C_clab_r_proxy, "ZZ"), "")
})

test_that("machine utilisation is read from the CPU accounting file", {
  ps <- tempfile()
  on.exit(unlink(ps), add = TRUE)
  mfp <- function(expr, ...) {
    co2fp(expr, loc = 30, cpupwrw = 45, memgb = 8, usage = "machine",
          procstat = ps, mthd = "model", ...)
  }
  # idle is the fourth and fifth numbers: 850 of 1000 ticks, then 1600 of
  # 1900, so 1 - 750 / 900 of the machine was busy
  writeLines("cpu  100 0 50 800 50 0 0 0 0 0", ps)
  fp <- mfp(writeLines("cpu  200 0 100 1500 100 0 0 0 0 0", ps),
            loadcur = "codecarbon")
  expect_identical(fp$usagemode, "machine")
  expect_equal(fp$usage, 1 / 6, tolerance = 1e-12)
  expect_identical(fp$assumptions,
                   "codecarbon machine load curve applied to the CPU power")
  # no change in the file: no time passed on the machine, utilisation 0
  fp <- mfp(1)
  expect_identical(fp$usage, 0)
  expect_identical(fp$assumptions, character(0))
  unread <- paste0("machine CPU accounting unreadable (no /proc/stat); ",
                   "process utilisation used")
  for (l in c("intr 5", "cpu 1 2 3", "cpu 1 2 x 4")) {
    writeLines(l, ps)
    fp <- mfp(1)
    expect_identical(fp$usagemode, "process")
    expect_identical(fp$assumptions, unread)
  }
  # the file is gone by the end of the run
  writeLines("cpu  100 0 50 800", ps)
  fp <- mfp(unlink(ps))
  expect_identical(fp$usagemode, "process")
  expect_identical(
    fp$assumptions,
    "machine CPU accounting became unreadable; process utilisation used"
  )
})

test_that("the kernel applies the machine curve and its fallback", {
  # an hour on 45 W at machine utilisation 1/6: codecarbon's curve gives
  # 45 (0.1 + 0.9 u^3) W; with the second reading missing, the process
  # utilisation 1800 / (3600 x 2) and the linear 45 u W instead
  run <- function(i1, i2) {
    inp <- c(3600, 1800, 2, 1, NA, 850, 1000, i1, i2, 45, 1, 0, 0.3725, 1,
             100, 0)
    .Call(C_clab_r_fpe, as.numeric(inp), numeric(0), numeric(0),
          numeric(0), numeric(0), numeric(0), integer(0))[[1L]]
  }
  o <- run(1600, 1900)
  expect_equal(o[2L], 1 / 6, tolerance = 1e-12)
  expect_equal(o[6L], 45 * (0.1 + 0.9 / 216) / 1000, tolerance = 1e-12)
  expect_identical(o[c(3L, 5L)], c(1, 16))
  o <- run(NA, NA)
  expect_identical(o[c(3L, 5L)], c(0, 17))
  expect_equal(o[2L], 0.25, tolerance = 1e-12)
  # the process curve: power floored at 10% of the CPU's
  expect_equal(o[6L], 45 * 0.25 / 1000, tolerance = 1e-12)
})

test_that("each argument of co2fp() is checked by its own message", {
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(co2fp(1, rapldir = NA), "co2fp: `rapldir` must be a single path.")
  ef(co2fp(1, procstat = 1), "co2fp: `procstat` must be a single path.")
  ef(co2fp(1, usage = c(0.1, 0.2)),
     paste0("co2fp: `usage` takes one value: \"process\", \"machine\" or a ",
            "single number in [0, 1]."))
  expect_identical(
    co2fp(1, loc = 30, cpupwrw = 45, memgb = 8, usage = NULL,
          mthd = "model")$usagemode,
    "process"
  )
  ef(co2fp(1, devices = 1),
     "co2fp: `devices` must be a data frame, one row per device.")
  dev <- data.frame(kind = "gpu", vendor = "v", name = "n", how = "h",
                    watts = "5", stringsAsFactors = FALSE)
  ef(co2fp(1, devices = dev), "co2fp: `devices$watts` must be numeric.")
  # a table with watts only: no energy column, so no measured energy
  dev$watts <- 5
  fp <- co2fp(1, loc = 30, cpupwrw = 45, memgb = 8, usage = 0.5,
              mthd = "model", devices = dev)
  expect_identical(fp$devices$note, "average power over the run")
  expect_equal(fp$devices$kwh, fp$durs / 3600 * 5 / 1000, tolerance = 1e-12)
})

test_that("unknown cores and memory are assumed and said so", {
  local_mocked_bindings(.memgb = function() NA_real_)
  fp <- co2fp(1, loc = 30, cpupwrw = 45, cores = 2, usage = 0.5,
              mthd = "model")
  expect_identical(fp$memgb, 8)
  expect_identical(fp$assumptions, "memory size unknown: 8 GB assumed")
  local_mocked_bindings(detectCores = function(...) NA_integer_,
                        .package = "parallel")
  fp <- co2fp(1, loc = 30, cpupwrw = 45, memgb = 4, usage = 0.5,
              mthd = "model")
  expect_identical(fp$cores, 1L)
  expect_identical(fp$assumptions, "core count unknown; 1 assumed")
})

test_that("RAPL counters are measured, wrapped, or given up", {
  skip_on_os("windows")
  d <- file.path(tempfile(), "intel-rapl:0")
  on.exit(unlink(dirname(d), recursive = TRUE), add = TRUE)
  dir.create(d, recursive = TRUE)
  writeLines("package-0", file.path(d, "name"))
  uj <- file.path(d, "energy_uj")
  rfp <- function(expr) {
    co2fp(expr, loc = 30, cpupwrw = 45, memgb = 8, usage = 0.5,
          rapldir = dirname(d))
  }
  # 3.6 J between the readings is 1e-6 kWh
  writeLines("1000000", uj)
  fp <- rfp(writeLines("4600000", uj))
  expect_identical(fp$method, "rapl")
  expect_equal(fp$enekwh$cpu, 1e-6, tolerance = 1e-12)
  writeLines("1000000", uj)
  fp <- rfp(writeLines("x", uj))
  expect_identical(fp$method, "model")
  expect_identical(
    fp$assumptions,
    "RAPL counters became unreadable during the run; modelled instead"
  )
  # a counter that went down, with no range to wrap it through
  writeLines("5000000", uj)
  fp <- rfp(writeLines("1000000", uj))
  expect_identical(fp$method, "model")
  expect_identical(
    fp$assumptions,
    "a RAPL counter wrapped and its range was unknown; modelled instead"
  )
  # with its range: 1e6 - 5e6 + 9e6 microjoules, 5 J
  writeLines("9000000", file.path(d, "max_energy_range_uj"))
  writeLines("5000000", uj)
  fp <- rfp(writeLines("1000000", uj))
  expect_identical(fp$method, "rapl")
  expect_equal(fp$enekwh$cpu, 5 / 3.6e6, tolerance = 1e-12)
})
