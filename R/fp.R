#' Carbon intensity of electricity by location, offline, anywhere
#'
#' The grams of CO2-equivalent emitted per kilowatt-hour of electricity at a
#' location, with the year and source of the figure. The table is compiled
#' into the package, so it works offline in any country: every country at the
#' latest year published by Our World in Data from Ember's yearly electricity
#' data (the latest year each grid reports: 2025 or 2024 for nearly all, older
#' for a few grids, Ukraine (2022) among them, named in \code{year}),
#' sub-national zones from Electricity Maps' 2024 yearly data as redistributed
#' in the Green Algorithms data release v3.1 (Canadian provinces and
#' territories, US balancing authorities, Australian states, Indian and
#' Japanese regions, Brazilian and Chilean systems, ...), and the world
#' average. All figures are lifecycle intensities (generation plus upstream),
#' the basis both sources publish.
#'
#' A location is matched by its code (\code{"CA"}, \code{"CA-ON"},
#' \code{"US-CAL-CISO"}), its ISO 3166-1 alpha-3 code (\code{"CAN"}) or its
#' name (\code{"Canada"}), case-insensitively. An unknown sub-national code
#' falls back to its country and says so in \code{note}. \code{NULL} detects
#' the location with \code{\link{dtctloc}}.
#'
#' Computing has a footprint of its own, and research that reports one
#' (as journals and funders increasingly ask) needs the grid's intensity:
#' \code{\link{co2fp}} multiplies a computation's energy by it.
#'
#' @section Provenance of the table:
#' The country rows are Ember's yearly electricity data as published by Our
#' World in Data; the sub-national zones and the Canadian provinces come from
#' the Electricity Maps yearly files redistributed in the Green Algorithms
#' data repository (v3.1). Every row names its source and year. Values are
#' reproduced unchanged: Ember reports three small diesel grids (Montserrat,
#' Saint Helena, Uzbekistan) at exactly 1000 gCO2e/kWh and Turkmenistan's gas
#' fleet at 1306, and the table does not round or cap them. Kosovo appears
#' once, under the code \code{XK}.
#'
#' Errors: a location the table does not know; a number that is not a single
#' non-negative value; text that is not valid in its encoding.
#'
#' @param loc A location code, alpha-3 code or name; a number in
#'   gCO2e/kWh, returned as a user-supplied value; or \code{NULL} to detect
#'   the location from the environment, time zone or locale.
#' @param x A result of \code{co2int()}.
#' @inheritParams freqtbl
#' @return \code{co2int()}: a list of class \code{co2int} with
#'   \code{g_per_kwh}, \code{location}, \code{name}, \code{year},
#'   \code{basis}, \code{source} and \code{note} (a character vector, empty
#'   unless a fallback or detection was involved). \code{co2inttab()}: the
#'   full table as a data frame with the columns \code{location},
#'   \code{iso3}, \code{name}, \code{g_per_kwh}, \code{year}, \code{basis}
#'   and \code{source}.
#' @references Ember (2026). Yearly Electricity Data, via Our World in Data.
#'
#'   Electricity Maps (2025). 2024 Yearly Carbon Intensity Data, via Green
#'   Algorithms data v3.1.
#'
#'   Lannelongue, L., Grealey, J. and Inouye, M. (2021). Green Algorithms:
#'   quantifying the carbon footprint of computation. Advanced Science 8(12),
#'   2100707. \doi{10.1002/advs.202100707}
#' @seealso \code{\link{co2fp}}, \code{\link{dtctloc}}.
#' @examples
#' co2int("CA-ON")
#' co2int("Canada")$year
#' co2int("CA-XX")$note
#' co2int(120)$source
#' nrow(co2inttab())
#' # the Canadian provinces in the table, lowest intensity first
#' tab <- co2inttab()
#' ca <- tab[
#'   startsWith(tab$location, "CA-"),
#'   c("location", "name", "g_per_kwh", "year")
#' ]
#' head(ca[order(ca$g_per_kwh), ])
#' # a location the table does not know is an error, not a silent world
#' # average
#' try(co2int("Atlantis"))
#' @export
co2int <- function(loc = "WORLD") .ciget(loc, "co2int")

#' @rdname co2int
#' @export
co2inttab <- function() .citab()

#' @rdname co2int
#' @export
print.co2int <- function(x, ...) {
  cat(
    "Carbon intensity of electricity: ", format(x$g_per_kwh), " gCO2e/kWh\n",
    "  Location: ", x$location, " (", x$name, ")",
    if (length(x$year) && !is.na(x$year)) paste0(", ", x$year), "\n",
    "  Basis: ", x$basis, "\n",
    "  Source: ", x$source, "\n",
    sep = ""
  )
  for (n in x$note) cat("  Note: ", n, "\n", sep = "")
  invisible(x)
}

.fpcache <- new.env(parent = emptyenv())

#' Where is this computation running?
#'
#' Finds a location code for \code{\link{co2int}} without any network access,
#' in this order: \enumerate{ \item the \code{CRIMLAB_LOCATION} environment
#' variable (any code or name the table knows); \item
#' \code{CODECARBON_COUNTRY_ISO_CODE} (alpha-3), the codecarbon variable, so a
#' machine already configured for codecarbon needs nothing more; \item the
#' system time zone, mapped to its country through the IANA zone tables
#' (\code{America/Toronto} is Canada, \code{Europe/Berlin} Germany); \item
#' the territory in the locale (\code{en_CA.UTF-8} is Canada). } When none of
#' these yields a location in the table the answer is \code{"WORLD"}, and the
#' method says so, so a report never silently applies one country's grid to
#' another's computation.
#'
#' The zone table, built from the IANA database, holds every zone of
#' \code{zone.tab}, which names exactly one country per zone
#' (\code{Europe/Stockholm} is Sweden even though the database now keeps its
#' rules under \code{Europe/Berlin}), plus the aliases and old names of the
#' \code{backward} file under their own country: the country of the zone the
#' alias stood for before the database merged it (\code{Iceland} is
#' \code{Atlantic/Reykjavik}, Iceland, not Abidjan's Ivory Coast), the link
#' target's country when that target serves a single country
#' (\code{Asia/Calcutta}, \code{US/Eastern}), and a reviewed list for the
#' rest (\code{Atlantic/Jan_Mayen}, \code{America/Coral_Harbour},
#' \code{Canada/Eastern}). Every geographic zone \code{OlsonNames()} lists
#' resolves to a country or territory; a territory the carbon table has no
#' row of its own for (Monaco, Jersey, Vaduz, the Antarctic stations) uses
#' the grid of the country it draws on or the world average, and the method
#' says which. The \code{Etc/}, \code{UTC}, fixed-offset and generic zones
#' (\code{CET}, \code{EST}, \code{WET}) name no country and fall through to
#' the locale. A variable set to a value the table does not know is ignored
#' with a warning.
#'
#' Errors: an input that is not one string or \code{NULL}.
#'
#' @param envloc,codloc,tz,locale The inputs, exposed so a caller or a test
#'   can supply them; the defaults read the running session.
#' @return A list with \code{location} (a code present in
#'   \code{co2inttab()}) and \code{method} (how it was found).
#' @seealso \code{\link{co2int}}, \code{\link{co2fp}}.
#' @examples
#' dtctloc()
#' dtctloc(envloc = "", codloc = "", tz = "Europe/Paris")
#' dtctloc(envloc = "", codloc = "", tz = "", locale = "en_IN.UTF-8")
#' # an explicit setting wins over everything else
#' dtctloc(envloc = "CA-QC", tz = "Europe/Paris")
#' # nothing usable: the world average, and the method says so
#' dtctloc(envloc = "", codloc = "", tz = "UTC", locale = "C")
#' @export
dtctloc <- function(
  envloc = Sys.getenv("CRIMLAB_LOCATION"),
  codloc = Sys.getenv("CODECARBON_COUNTRY_ISO_CODE"),
  tz = Sys.timezone(),
  locale = Sys.getlocale("LC_TIME")
) {
  f <- "dtctloc"
  args <- list(envloc = envloc, codloc = codloc, tz = tz, locale = locale)
  for (nm in names(args)) {
    v <- args[[nm]]
    if (!is.null(v) && !(is.character(v) && length(v) <= 1L)) {
      .ckstop(f, "`%s` must be one string (or NULL)", nm)
    }
    .cktxt(v, nm, f)
  }
  r <- .locenv(envloc, codloc, f)
  if (is.null(r)) r <- .loctz(tz)
  if (is.null(r)) r <- .locloc(locale)
  if (is.null(r)) {
    r <- list(
      location = "WORLD",
      method = "fallback, no location found in environment, time zone or locale"
    )
  }
  r
}

#' Is CPU energy measurable here?
#'
#' \code{TRUE} when the operating system exposes readable RAPL (Running
#' Average Power Limit) energy counters for at least one CPU package: Linux
#' with the \code{intel_rapl} powercap driver (Intel and, since kernel 5.8,
#' AMD), with read permission on \code{energy_uj}. On most systems since
#' Linux 5.10 that file is readable by root only, so \code{FALSE} is the
#' common answer for an ordinary user, and \code{\link{co2fp}} then falls
#' back to its model.
#'
#' @details
#' A domain counts when it is a directory named \code{intel-rapl:N} (the
#' \code{intel-rapl:N:M} subdomains are parts of a package and would count
#' the same energy twice) whose \code{name} file starts with
#' \code{"package"} and whose \code{energy_uj} file holds a number. Every such
#' domain must be readable: if any is not, the answer is \code{FALSE}. Only
#' files under \code{rapldir} are read; nothing is written. On macOS and
#' Windows there is no such directory and the answer is \code{FALSE}. Errors
#' when \code{rapldir} is not a single non-missing string.
#'
#' @param rapldir Character scalar, the powercap directory; default
#'   \code{"/sys/class/powercap"}, the Linux location. Another path is useful
#'   for testing against a copy of the directory tree.
#' @return A logical scalar: \code{TRUE} when every CPU-package energy counter
#'   under \code{rapldir} is readable, \code{FALSE} otherwise (including when
#'   the directory does not exist).
#' @seealso \code{\link{co2fp}}
#' @examples
#' raplok()
#'
#' # a fake powercap tree with one readable package domain (Linux names the
#' # domain "intel-rapl:0"; a colon is not allowed in a Windows file name)
#' if (.Platform$OS.type != "windows") {
#'   d <- file.path(tempfile(), "intel-rapl:0")
#'   dir.create(d, recursive = TRUE)
#'   writeLines("package-0", file.path(d, "name"))
#'   writeLines("123456789", file.path(d, "energy_uj"))
#'   raplok(dirname(d))
#' }
#'
#' # a folder that does not exist is simply FALSE; a non-path is an error
#' raplok(file.path(tempdir(), "no-such-dir"))
#' try(raplok(NA_character_))
#' @export
raplok <- function(rapldir = "/sys/class/powercap") {
  if (!is.character(rapldir) || length(rapldir) != 1L || is.na(rapldir)) {
    .ckstop("raplok", "`rapldir` must be a single directory path")
  }
  .ckpath(rapldir, "raplok", "rapldir")
  !is.null(.raplrd(.rapldom(rapldir)))
}

#' Energy and carbon footprint of a computation
#'
#' Evaluates \code{expr} and reports the electricity it used and the
#' CO2-equivalent that electricity implies, by one of two methods, for the
#' grid it ran on, offline. Use it to report the energy and emissions of an
#' analysis next to its results (as journals and funders increasingly ask),
#' or to compare two ways of running the same computation. The result keeps
#' the value of \code{expr}, and lists every figure that was assumed rather
#' than measured.
#'
#' @details
#' \strong{Measured} (\code{meth = "rapl"}): the CPU package energy
#' counters (RAPL, microjoules, read from \code{rapldir}) before and after,
#' with wrap-around handled through \code{max_energy_range_uj}. This is what
#' codecarbon does on Linux; the counters cover the CPU packages only, so
#' memory is still modelled. Available when \code{\link{raplok}()} is
#' \code{TRUE}.
#'
#' \strong{Modelled} (\code{meth = "model"}): the Green Algorithms formula
#' (Lannelongue, Grealey and Inouye 2021),
#' \deqn{E = t \,(P_{cpu}\, u + M\, P_{mem})\, PUE}
#' with \eqn{t} the running time, \eqn{P_{cpu}} the CPU's thermal design
#' power in watts, \eqn{u} its utilisation, \eqn{M} the memory in GB,
#' \eqn{P_{mem} = 0.3725} W/GB and \eqn{PUE} the data-centre overhead. When
#' the CPU's power is unknown (\code{cpupwrw = NULL}) the fallback is
#' codecarbon's: 50\% of a constant 85 W, that is 42.5 W charged for the whole
#' running time whatever the measured utilisation (the utilisation is still
#' reported, but not applied).
#'
#' \strong{Utilisation} is measured, not assumed, in one of two ways.
#' \code{usage = "process"} (the default) is this R process's own CPU time,
#' children included, divided by the wall time and the core count:
#' codecarbon's process tracking mode, which charges only what the
#' computation itself used. \code{usage = "machine"} is the whole machine's
#' busy share over the run, from \code{/proc/stat} on Linux: codecarbon's
#' machine tracking mode, right when the computation is the only thing the
#' machine is doing or when the question is what the machine drew. Where
#' \code{/proc/stat} is absent the process figure is used and the assumptions
#' say so. A number between 0 and 1 fixes the utilisation, as the Green
#' Algorithms calculator asks for it. \code{loadcur} chooses how power scales
#' with that utilisation: linearly (Green Algorithms) or along codecarbon's
#' curves, a floor at 10\% of \eqn{P_{cpu}} for a process and
#' \eqn{0.1 + 0.9 u^3} for a machine, so a figure can be made comparable with
#' either tool.
#'
#' \strong{Other devices} (\code{devices}): the CPU counters and the model
#' cover the CPU and memory only. The energy of a graphics card or another
#' accelerator, measured by its own tools, is added through a table with one
#' row per device and the columns \code{kind} (such as \code{"gpu"}),
#' \code{vendor}, \code{name}, \code{watts} (the device's average power over
#' the run, in W) and \code{energy} (its energy over the run, in kWh) (one of
#' the two may be left out, and a missing value means not measured) and
#' \code{how} (how the figure was obtained). A device's energy is its
#' \code{energy}, or else \code{watts} times the running time; that energy is
#' added to the CPU and memory energy before the PUE and the carbon
#' intensity are applied. A row of kind \code{"cpu"} or \code{"memory"} is
#' not counted again, nor is a row with neither figure; the result's
#' \code{devices} table says which rows were counted and why the others were
#' not. Without \code{devices} the result is as described below, with no
#' \code{devices} element.
#'
#' Emissions are \eqn{E \times} the carbon intensity of the location
#' (\code{\link{co2int}}: a table compiled into the package of every
#' country at its latest published year and of sub-national zones, so no
#' network access is needed anywhere). With \code{loc = NULL} the location is
#' detected from the environment, the time zone or the locale
#' (\code{\link{dtctloc}}) and the detection is recorded. Every default that
#' stood in for a measurement is listed in \code{assumptions}, so a report can
#' say what was measured and what was assumed.
#'
#' The energy figures include the PUE factor (CPU, memory and total alike).
#' The running time has a floor of the CPU time divided by the cores, so a
#' computation faster than the clock's resolution does not come out at zero.
#' Short computations (well under a second) give figures dominated by timer
#' resolution; measure a representative workload.
#'
#' Every argument is checked, and the carbon intensity looked up, before
#' \code{expr} is evaluated, so a bad argument never costs the computation.
#' Errors: an unknown \code{loc}; \code{pue} below 1 or not a single finite
#' number; \code{cpupwrw} or \code{cores} not a single positive number
#' (\code{cores} also whole); \code{memgb} or \code{mempwrwpergb} negative;
#' \code{usage} other than \code{"process"}, \code{"machine"} or one number in
#' \eqn{[0, 1]}; \code{meth = "rapl"} where the counters are not readable;
#' a \code{devices} table without its columns or with a negative or infinite
#' figure. An error inside \code{expr} propagates unchanged. Nothing is
#' written and nothing is fetched.
#'
#' @param expr The expression to evaluate. It is evaluated once, in the
#'   calling environment, when the function reaches it (lazily, as a
#'   promise).
#' @param loc The location, passed to \code{\link{co2int}}. It is a location
#'   code, alpha-3 code or name, a number in gCO2e/kWh, or \code{NULL} (the
#'   default) to detect it.
#' @param meth Character scalar: \code{"auto"} (the default; RAPL when
#'   readable, else the model), \code{"rapl"} (measured; an error when the
#'   counters are not readable) or \code{"model"} (always the Green
#'   Algorithms model).
#' @param cpupwrw The CPU's total thermal design power in watts (the model's
#'   \eqn{P_{cpu}}); \code{NULL} uses codecarbon's unknown-CPU fallback.
#' @param usage \code{"process"} (measured from this process's CPU time),
#'   \code{"machine"} (measured from the whole machine's CPU accounting), or a
#'   number in between 0 and 1. \code{NULL} means \code{"process"}.
#' @param loadcur Character scalar: \code{"linear"} (the default; power
#'   proportional to utilisation, as Green Algorithms) or \code{"codecarbon"}
#'   (codecarbon's curves, see Details). Applies only when \code{cpupwrw} is
#'   given.
#' @param cores Number of logical cores the process utilisation is measured
#'   against; \code{NULL} detects it.
#' @param memgb Memory charged to the computation in GB; \code{NULL} detects
#'   the machine's total memory where it can, else assumes 8 GB.
#' @param mempwrwpergb Memory power per GB (default 0.3725 W/GB, Green
#'   Algorithms; codecarbon v2 used 0.375).
#' @param pue Power usage effectiveness: 1 for a laptop or desktop (default);
#'   1.67 is the Green Algorithms default for a data centre, 1.56 the Uptime
#'   Institute 2024 global average. Must be a single number of at least 1.
#' @param rapldir Character scalar, the powercap directory holding the RAPL
#'   domains; default \code{"/sys/class/powercap"} (Linux).
#' @param procstat Character scalar, the kernel CPU accounting file read for
#'   \code{usage = "machine"}; default \code{"/proc/stat"} (Linux).
#' @param devices \code{NULL} (the default), or a data frame of the other
#'   devices' measured energy, one row per device: \code{kind},
#'   \code{vendor}, \code{name}, \code{watts} and/or \code{energy} (kWh), and
#'   \code{how} (see Details).
#' @param x A result of \code{co2fp()} or \code{fpeqv()}.
#' @inheritParams freqtbl
#' @return A list of class \code{co2fp} with elements \code{value} (the
#'   result of \code{expr}), \code{durs} (numeric, wall seconds),
#'   \code{cputms} (numeric, CPU seconds of this process and its children),
#'   \code{usage} (numeric in \eqn{[0, 1]}, the utilisation),
#'   \code{usagemode} (\code{"process"}, \code{"machine"} or
#'   \code{"fixed"}), \code{cores} (integer), \code{memgb} (numeric),
#'   \code{enekwh} (a list of numerics \code{cpu}, \code{memory} and
#'   \code{total}, in kWh, PUE included, and \code{devices} when
#'   \code{devices} is given), \code{co2eg} (numeric, grams of
#'   CO2-equivalent), \code{co2int} (the list returned by
#'   \code{\link{co2int}}, with its year and source), \code{method}
#'   (\code{"rapl"} or \code{"model"}, the method actually used),
#'   \code{assumptions} (character vector of every default or fallback that
#'   stood in for a measurement) and \code{citations} (character vector).
#'   With \code{devices}, also \code{devices}: a data frame with \code{kind},
#'   \code{vendor}, \code{name}, \code{how}, \code{kwh} (the device's energy
#'   before PUE, \code{NA} when not counted), \code{counted} (logical) and
#'   \code{note} (how the energy was taken, or why it was not counted). The
#'   print method shows a short summary.
#' @references Lannelongue, L., Grealey, J. and Inouye, M. (2021). Green
#'   Algorithms: quantifying the carbon footprint of computation. Advanced
#'   Science 8(12), 2100707. \doi{10.1002/advs.202100707}
#'
#'   Courty, B. et al. (2024). CodeCarbon: estimate and track carbon
#'   emissions from machine learning computing. \doi{10.5281/zenodo.4658424}
#'
#'   Lacoste, A., Luccioni, A., Schmidt, V. and Dandres, T. (2019).
#'   Quantifying the carbon emissions of machine learning. arXiv:1910.09700.
#' @seealso \code{\link{co2int}}, \code{\link{raplok}}, \code{\link{fpeqv}}.
#' @examples
#' fp <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON")
#' fp$enekwh$total
#' fp$co2eg
#' fp$co2int$year
#' fp$assumptions
#' print(fp)
#' # wherever you are, offline: the location is detected and recorded
#' co2fp(1 + 1, cpupwrw = 45, memgb = 8)$assumptions
#' # the whole machine's utilisation, on codecarbon's machine curve
#' co2fp(
#'   1 + 1,
#'   loc = "FR",
#'   usage = "machine",
#'   loadcur = "codecarbon",
#'   cpupwrw = 45,
#'   memgb = 8
#' )$usagemode
#' # a fixed utilisation and a data-centre PUE, with a carbon intensity given
#' # as a number
#' fp2 <- co2fp(
#'   sum(sqrt(seq_len(2e5))),
#'   loc = 30,
#'   usage = 0.5,
#'   cpupwrw = 65,
#'   memgb = 16,
#'   pue = 1.67,
#'   meth = "model"
#' )
#' fp2$enekwh
#' fpeqv(fp2$co2eg)
#' # a graphics card measured at 120 W on average (a simulated reading) and
#' # a CPU row, which is not counted twice
#' dev <- data.frame(
#'   kind = c("gpu", "cpu"), vendor = c("example", "example"),
#'   name = c("card 0", "package 0"), watts = c(120, 35),
#'   how = c("simulated reading", "simulated reading")
#' )
#' fp3 <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON", cpupwrw = 65,
#'              memgb = 16, meth = "model", devices = dev)
#' fp3$devices
#' # a bad argument is refused before the expression runs
#' try(co2fp(Sys.sleep(10), pue = 0.9))
#' @export
co2fp <- function(
  expr,
  loc = NULL,
  meth = c("auto", "rapl", "model"),
  cpupwrw = NULL,
  usage = "process",
  loadcur = c("linear", "codecarbon"),
  cores = NULL,
  memgb = NULL,
  mempwrwpergb = 0.3725,
  pue = 1,
  rapldir = "/sys/class/powercap",
  procstat = "/proc/stat",
  devices = NULL
) {
  f <- "co2fp"
  meth <- .ckarg(meth, c("auto", "rapl", "model"), "meth", f)
  loadcur <- .ckarg(loadcur, c("linear", "codecarbon"), "loadcur", f)
  ci <- .ciget(loc, f)
  # every argument is checked before `expr` runs, so a bad argument never
  # costs the computation
  a <- .fpargs(pue, cpupwrw, cores, memgb, mempwrwpergb, rapldir, procstat,
               usage, f)
  dv <- .fpdev(devices, f)
  rp <- .fprapl(meth, rapldir, f)
  md <- .fpmode(a$mode, procstat, ci$note)
  mode <- md$mode
  t0 <- proc.time()
  value <- expr
  t1 <- proc.time()
  after <- if (!is.null(rp$before)) .raplrd(rp$domains) else NULL
  st1 <- if (mode == "machine") .pstat(procstat) else NULL
  parts <- c("user.self", "sys.self", "user.child", "sys.child")
  d <- unname(t1[parts]) - unname(t0[parts])
  cpu <- max(sum(d[is.finite(d)]), 0)
  cr <- .fpcores(a$cores, md$asm)
  mem <- .fpmem(a$memgb)
  inp <- c(
    unname(t1[["elapsed"]] - t0[["elapsed"]]), cpu, cr$cores,
    match(mode, c("process", "machine", "fixed")) - 1L,
    if (mode == "fixed") a$usage else NA_real_,
    if (is.null(md$st0)) c(NA_real_, NA_real_) else md$st0[c("idle", "total")],
    if (is.null(st1)) c(NA_real_, NA_real_) else st1[c("idle", "total")],
    if (is.null(a$cpupwrw)) NA_real_ else a$cpupwrw,
    as.numeric(loadcur == "codecarbon"), mem$memgb, a$mempwrwpergb, pue,
    ci$g_per_kwh, as.numeric(!is.null(rp$before))
  )
  r <- .fpcall(inp, rp$before, after, dv)
  o <- r[[1L]]
  umode <- c("process", "machine", "fixed")[o[3L] + 1L]
  asm <- c(.fpflags(as.integer(o[5L]), umode, cr$asm), mem$note)
  if (pue > 1) asm <- c(asm, sprintf("PUE %.2f applied", pue))
  out <- list(
    value = value, durs = o[1L], cputms = cpu, usage = o[2L],
    usagemode = umode, cores = as.integer(cr$cores), memgb = mem$memgb,
    enekwh = .fpen(o, devices), co2eg = o[10L], co2int = ci,
    method = c("rapl", "model")[o[4L] + 1L], assumptions = asm,
    citations = .fpcites(ci)
  )
  if (!is.null(devices)) out$devices <- .fpdevout(dv, r[[2L]])
  structure(out, class = c("co2fp", "list"))
}

#' @rdname co2fp
#' @export
print.co2fp <- function(x, ...) {
  ci <- x$co2int
  cat(
    sprintf(
      "Computation footprint (%s)\n",
      if (identical(x$method, "rapl")) "measured, RAPL" else "modelled"
    ),
    sprintf(
      "  Time: %.3f s wall, %.3f s CPU; utilisation %.2f (%s) on %d cores\n",
      x$durs, x$cputms, x$usage, x$usagemode, x$cores
    ),
    sprintf(
      "  Energy: %.3g kWh (cpu %.3g, memory %.3g%s)\n",
      x$enekwh$total, x$enekwh$cpu, x$enekwh$memory,
      if (is.null(x$enekwh$devices)) "" else
        sprintf(", devices %.3g", x$enekwh$devices)
    ),
    sprintf(
      "  CO2e: %.3g g at %.4g gCO2e/kWh (%s%s)\n",
      x$co2eg, ci$g_per_kwh, ci$name,
      if (is.na(ci$year)) "" else sprintf(", %d", ci$year)
    ),
    sep = ""
  )
  if (length(x$assumptions)) {
    cat("  Assumed: ", paste(x$assumptions, collapse = "; "), "\n", sep = "")
  }
  invisible(x)
}

#' Everyday equivalents of a quantity of CO2-equivalent
#'
#' Converts grams of CO2-equivalent into two everyday comparisons: kilometres
#' driven in an average European passenger car (175 gCO2/km) and tree-months
#' of carbon sequestration (917 gCO2 per tree per month), the two comparisons
#' the Green Algorithms calculator prints, with their sources. Use it to make
#' the \code{co2eg} of \code{\link{co2fp}} easier to read in a report.
#'
#' @details
#' The conversions are plain divisions: \eqn{\mathrm{km} = g / 175} and
#' \eqn{\mathrm{tree\ months} = g / 917}. Both factors are averages and the
#' comparisons are illustrations, not offsets: a tree-month is the carbon one
#' mature tree takes up in a month, not a promise that planting one cancels
#' the emission. Errors when \code{co2eg} is not a single finite, non-negative
#' number.
#'
#' @param co2eg Numeric scalar, grams of CO2-equivalent; finite and
#'   non-negative (for example the \code{co2eg} element of a
#'   \code{\link{co2fp}} result).
#' @param x A result of \code{fpeqv()}.
#' @inheritParams freqtbl
#' @return A list of class \code{fpeqv} with \code{carkm} (numeric,
#'   kilometres), \code{treemon} (numeric) and \code{sources} (named character
#'   vector, \code{car} and \code{tree}, naming each factor's source).
#' @references Green Algorithms data v3.1, context table: car factor from
#'   Helmers et al. (2019) Atmospheric Environment and UK BEIS 2019 conversion
#'   factors; tree sequestration from Bernal et al. (2018) via the same table.
#' @seealso \code{\link{co2fp}}
#' @examples
#' fpeqv(1200)
#'
#' # the footprint of a small computation, in everyday terms
#' fp <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON", cpupwrw = 45,
#'             memgb = 8)
#' unlist(fpeqv(fp$co2eg)[c("carkm", "treemon")])
#'
#' # a vector is refused: convert one total at a time
#' try(fpeqv(c(10, 20)))
#' @export
fpeqv <- function(co2eg) {
  if (
    !is.numeric(co2eg) || length(co2eg) != 1L || !is.finite(co2eg) ||
      co2eg < 0
  ) {
    .ckstop("fpeqv", "`co2eg` must be a single non-negative number")
  }
  e <- .Call(C_clab_r_fpeqv, as.numeric(co2eg))
  structure(
    list(
      carkm = e[1L],
      treemon = e[2L],
      sources = c(
        car = paste0(
          "175 gCO2/km, European passenger car ",
          "(Green Algorithms context.csv)"
        ),
        tree = "917 gCO2 per tree-month (Green Algorithms context.csv)"
      )
    ),
    class = c("fpeqv", "list")
  )
}

#' @rdname fpeqv
#' @export
print.fpeqv <- function(x, ...) {
  cat(
    "Footprint in everyday terms\n",
    "  Car travel: ", format(x$carkm, digits = 4), " km\n",
    "  Tree uptake: ", format(x$treemon, digits = 4), " tree-months\n",
    sep = ""
  )
  invisible(x)
}

.citab <- function() {
  if (is.null(.fpcache$tab)) .fpcache$tab <- .Call(C_clab_r_citab)
  .fpcache$tab
}

.cimatch <- function(key) .Call(C_clab_r_cimatch, enc2utf8(key))

.ciget <- function(loc, fn) {
  .cktxt(loc, "loc", fn)
  if (is.numeric(loc)) {
    return(.ciuser(loc, fn))
  }
  note <- character(0)
  if (is.null(loc)) {
    det <- dtctloc()
    loc <- det$location
    note <- sprintf("location %s via %s", det$location, det$method)
  }
  if (!is.character(loc) || length(loc) != 1L || is.na(loc)) {
    .ckstop(
      fn,
      paste(
        "`loc` must be a single location code or name, a number in",
        "gCO2e/kWh, or NULL"
      )
    )
  }
  tab <- .citab()
  i <- .cimatch(loc)
  if (is.na(i) && grepl("-", loc, fixed = TRUE)) {
    i <- .cimatch(sub("-.*$", "", loc))
    if (!is.na(i)) {
      note <- c(
        note,
        sprintf(
          "no row for zone %s; its country %s used",
          toupper(loc),
          tab$location[i]
        )
      )
    }
  }
  if (is.na(i)) {
    .ckstop(
      fn,
      paste(
        "Unknown location '%s'. Use a code, alpha-3 code or name from",
        "co2inttab(), or pass a number in gCO2e/kWh"
      ),
      loc
    )
  }
  .cirow(tab, i, note)
}

# .ciget(): a carbon intensity the caller gave as a number
.ciuser <- function(loc, fn) {
  if (length(loc) != 1L || !is.finite(loc) || loc < 0) {
    .ckstop(
      fn,
      "a numeric `loc` must be a single non-negative gCO2e/kWh value"
    )
  }
  structure(
    list(
      g_per_kwh = as.numeric(loc),
      location = "user",
      name = "user-supplied",
      year = NA_integer_,
      basis = "user",
      source = "user-supplied carbon intensity",
      note = character(0)
    ),
    class = c("co2int", "list")
  )
}

# .ciget(): row i of the table of intensities
.cirow <- function(tab, i, note) {
  structure(
    list(
      g_per_kwh = tab$g_per_kwh[i],
      location = tab$location[i],
      name = tab$name[i],
      year = tab$year[i],
      basis = tab$basis[i],
      source = tab$source[i],
      note = note
    ),
    class = c("co2int", "list")
  )
}

.fpknown <- function(key) {
  if (is.null(key) || length(key) != 1L || is.na(key) || !nzchar(key)) {
    return(NA_character_)
  }
  i <- .cimatch(key)
  if (is.na(i)) NA_character_ else .citab()$location[i]
}

.locign <- function(val, var, fn) {
  if (length(val) == 1L && !is.na(val) && nzchar(val)) {
    warning(
      sprintf(
        paste(
          "%s: %s='%s' is not a location co2inttab() knows",
          "(code, alpha-3 or name); it is ignored."
        ),
        fn,
        var,
        val
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.locenv <- function(envloc, codloc, fn) {
  loc <- .fpknown(envloc)
  if (!is.na(loc)) {
    return(list(
      location = loc,
      method = "CRIMLAB_LOCATION environment variable"
    ))
  }
  .locign(envloc, "CRIMLAB_LOCATION", fn)
  loc <- .fpknown(codloc)
  if (!is.na(loc)) {
    return(list(
      location = loc,
      method = "CODECARBON_COUNTRY_ISO_CODE environment variable"
    ))
  }
  .locign(codloc, "CODECARBON_COUNTRY_ISO_CODE", fn)
  NULL
}

.loctz <- function(tz) {
  if (!(length(tz) == 1L && !is.na(tz) && nzchar(tz))) {
    return(NULL)
  }
  terr <- .Call(C_clab_r_tzterr, enc2utf8(tz))
  loc <- .fpknown(terr)
  if (!is.na(loc)) {
    return(list(location = loc, method = sprintf("system time zone %s", tz)))
  }
  if (is.na(terr)) {
    return(NULL)
  }
  # a territory the carbon table has no row for (Monaco, Jersey, the
  # Antarctic stations): its administering or neighbouring country's grid,
  # or the world average, never the locale's country
  proxy <- .Call(C_clab_r_proxy, terr)
  loc <- if (nzchar(proxy)) .fpknown(proxy) else NA_character_
  if (is.na(loc)) loc <- "WORLD"
  list(
    location = loc,
    method = sprintf(
      "system time zone %s (territory %s has no row of its own; %s)",
      tz,
      terr,
      if (identical(loc, "WORLD")) "world average" else paste("grid of", loc)
    )
  )
}

.locloc <- function(locale) {
  if (!(length(locale) == 1L && !is.na(locale))) {
    return(NULL)
  }
  m <- regmatches(locale, regexpr("^[A-Za-z]{2,3}_[A-Z]{2}", locale))
  if (length(m) != 1L) {
    return(NULL)
  }
  loc <- .fpknown(sub("^.*_", "", m))
  if (is.na(loc)) {
    return(NULL)
  }
  list(location = loc, method = sprintf("locale %s", locale))
}

# The CPU-package domains under a powercap directory.
.rapldom <- function(dir) {
  if (!dir.exists(dir)) {
    return(character(0))
  }
  d <- list.files(dir, pattern = "^intel-rapl:[0-9]+$", full.names = TRUE)
  keep <- vapply(
    d,
    function(p) {
      nm <- tryCatch(
        suppressWarnings(
          readLines(file.path(p, "name"), n = 1L, warn = FALSE)
        ),
        error = function(e) ""
      )
      file.exists(file.path(p, "energy_uj")) &&
        length(nm) == 1L &&
        startsWith(nm, "package")
    },
    logical(1)
  )
  d[keep]
}

# One number from the first line of a file, NA when there is none.
.rd1 <- function(path) {
  v <- tryCatch(
    suppressWarnings(as.numeric(readLines(path, n = 1L, warn = FALSE))),
    error = function(e) NA_real_
  )
  if (length(v) != 1L) NA_real_ else v
}

# The counters (microjoules) and their ranges; NULL when one is unreadable.
.raplrd <- function(domains) {
  if (!length(domains)) {
    return(NULL)
  }
  vals <- vapply(file.path(domains, "energy_uj"), .rd1, numeric(1))
  if (anyNA(vals)) {
    return(NULL)
  }
  maxes <- vapply(file.path(domains, "max_energy_range_uj"), .rd1, numeric(1))
  list(uj = unname(vals), max = unname(maxes))
}

# The machine's idle and total ticks from the kernel's CPU accounting.
.pstat <- function(path) {
  l <- tryCatch(
    suppressWarnings(readLines(path, n = 1L, warn = FALSE)),
    error = function(e) character(0)
  )
  if (length(l) != 1L || !startsWith(l, "cpu ")) {
    return(NULL)
  }
  v <- suppressWarnings(as.numeric(strsplit(trimws(l), "\\s+")[[1L]][-1L]))
  if (length(v) < 4L || anyNA(v)) {
    return(NULL)
  }
  c(idle = v[4L] + if (length(v) >= 5L) v[5L] else 0, total = sum(v))
}

# The machine's total memory in GB, NA when it cannot be read.
.memgb <- function() {
  if (file.exists("/proc/meminfo")) {
    l <- tryCatch(
      readLines("/proc/meminfo", warn = FALSE),
      error = function(e) character(0)
    )
    m <- grep("^MemTotal:", l, value = TRUE)
    kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", m)))
    if (length(kb) == 1L && is.finite(kb) && kb > 0) return(kb / 1024^2)
  }
  sysctl <- Sys.which("sysctl")
  if (identical(Sys.info()[["sysname"]], "Darwin") && nzchar(sysctl)) {
    b <- tryCatch(
      suppressWarnings(as.numeric(system2(
        sysctl,
        c("-n", "hw.memsize"),
        stdout = TRUE,
        timeout = 10
      ))),
      error = function(e) NA_real_
    )
    if (length(b) == 1L && is.finite(b) && b > 0) return(b / 1024^3)
  }
  NA_real_
}

# co2fp(): the energy by part, the devices' when there are devices
.fpen <- function(o, devices) {
  en <- list(cpu = o[6L], memory = o[7L], total = o[9L])
  if (!is.null(devices)) en$devices <- o[8L]
  en
}

# co2fp(): the RAPL domains and their counters before the run (NULL when
# unreadable, an error when RAPL was asked for)
.fprapl <- function(meth, rapldir, f) {
  domains <- if (meth != "model") .rapldom(rapldir) else character(0)
  before <- if (meth != "model") .raplrd(domains) else NULL
  if (meth == "rapl" && is.null(before)) {
    .ckstop(
      f,
      paste0(
        "RAPL energy counters are not readable here (see raplok()); use ",
        "meth = \"model\""
      )
    )
  }
  list(domains = domains, before = before)
}

# co2fp(): the utilisation mode and the machine's CPU accounting before the
# run, process utilisation when /proc/stat cannot be read
.fpmode <- function(mode, procstat, asm) {
  st0 <- if (mode == "machine") .pstat(procstat) else NULL
  if (mode == "machine" && is.null(st0)) {
    mode <- "process"
    asm <- c(
      asm,
      paste0(
        "machine CPU accounting unreadable (no /proc/stat); process ",
        "utilisation used"
      )
    )
  }
  list(mode = mode, st0 = st0, asm = asm)
}

# co2fp(): the cores charged, 1 when the count is unknown
.fpcores <- function(cores, asm) {
  if (is.null(cores)) {
    cores <- tryCatch(
      parallel::detectCores(logical = TRUE),
      error = function(e) NA_integer_
    )
    if (is.na(cores)) {
      cores <- 1L
      asm <- c(asm, "core count unknown; 1 assumed")
    }
  }
  list(cores = cores, asm = asm)
}

# co2fp(): the memory charged, the machine's total by default, 8 GB when
# that is unknown
.fpmem <- function(memgb) {
  note <- NULL
  if (is.null(memgb)) {
    memgb <- .memgb()
    if (is.na(memgb)) {
      memgb <- 8
      note <- "memory size unknown: 8 GB assumed"
    } else {
      note <- sprintf("memory charged: the machine's total %.1f GB", memgb)
    }
  }
  list(memgb = memgb, note = note)
}

# co2fp(): the energy and emissions from the core
.fpcall <- function(inp, before, after, dv) {
  nr <- if (is.null(before)) 0L else length(before$uj)
  r1 <- if (is.null(after)) rep(NA_real_, nr) else after$uj
  .Call(
    C_clab_r_fpe, as.numeric(unname(inp)),
    if (nr) as.numeric(before$uj) else numeric(0),
    as.numeric(r1),
    if (nr) as.numeric(before$max) else numeric(0),
    dv$watts, dv$energy, dv$kind
  )
}

# co2fp(): the assumptions the core's flags record
.fpflags <- function(flags, umode, asm) {
  if (bitwAnd(flags, 1L)) {
    asm <- c(asm, paste(
      "machine CPU accounting became unreadable; process utilisation used"
    ))
  }
  if (bitwAnd(flags, 4L)) {
    asm <- c(
      asm,
      "RAPL counters became unreadable during the run; modelled instead"
    )
  }
  if (bitwAnd(flags, 2L)) {
    asm <- c(
      asm,
      "a RAPL counter wrapped and its range was unknown; modelled instead"
    )
  }
  if (bitwAnd(flags, 8L)) {
    asm <- c(asm, "CPU power unknown: 85 W x 50% (codecarbon fallback)")
  }
  if (bitwAnd(flags, 16L)) {
    asm <- c(
      asm,
      sprintf(
        "codecarbon %s load curve applied to the CPU power",
        if (umode == "machine") "machine" else "process"
      )
    )
  }
  asm
}

.fpnum <- function(x, what, fn, positive = FALSE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x < 0 ||
      x > 1e9 ||
      (positive && x <= 0)
  ) {
    .ckstop(
      fn,
      "`%s` must be a single %s number",
      what,
      if (positive) "positive" else "non-negative"
    )
  }
  as.numeric(x)
}

.fpargs <- function(pue, cpupwrw, cores, memgb, mempwrwpergb, rapldir,
                    procstat, usage, fn) {
  if (!is.numeric(pue) || length(pue) != 1L || !is.finite(pue) || pue < 1) {
    .ckstop(fn, "`pue` must be a single number >= 1")
  }
  cpupwrw <- .fpnum(cpupwrw, "cpupwrw", fn, positive = TRUE)
  cores <- .fpnum(cores, "cores", fn, positive = TRUE)
  if (!is.null(cores) && cores != round(cores)) {
    .ckstop(fn, "`cores` must be a whole number")
  }
  memgb <- .fpnum(memgb, "memgb", fn)
  mempwrwpergb <- .fpnum(mempwrwpergb, "mempwrwpergb", fn)
  if (!is.character(rapldir) || length(rapldir) != 1L || is.na(rapldir)) {
    .ckstop(fn, "`rapldir` must be a single path")
  }
  if (!is.character(procstat) || length(procstat) != 1L || is.na(procstat)) {
    .ckstop(fn, "`procstat` must be a single path")
  }
  .ckpath(rapldir, fn, "rapldir")
  .ckpath(procstat, fn, "procstat")
  if (is.null(usage)) usage <- "process"
  if (length(usage) != 1L) {
    .ckstop(
      fn,
      paste0(
        "`usage` takes one value: \"process\", \"machine\" or a single ",
        "number in [0, 1]"
      )
    )
  }
  mode <- if (is.character(usage)) {
    .ckarg(usage, c("process", "machine"), "usage", fn)
  } else if (is.numeric(usage) && is.finite(usage) && usage >= 0 &&
               usage <= 1) {
    "fixed"
  } else {
    .ckstop(
      fn,
      "`usage` must be \"process\", \"machine\" or a single number in [0, 1]"
    )
  }
  list(
    cpupwrw = cpupwrw,
    cores = cores,
    memgb = memgb,
    mempwrwpergb = mempwrwpergb,
    usage = if (mode == "fixed") as.numeric(usage) else usage,
    mode = mode
  )
}

# The devices table, checked: kinds as codes (1 for a CPU or memory row,
# counted already), figures as numbers (NA: not measured).
.fpdev <- function(d, fn) {
  if (is.null(d)) {
    return(list(kind = integer(0), watts = numeric(0), energy = numeric(0)))
  }
  if (!is.data.frame(d)) {
    .ckstop(fn, "`devices` must be a data frame, one row per device")
  }
  need <- c("kind", "vendor", "name", "how")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    .ckstop(fn, "`devices` is missing columns: %s", toString(miss))
  }
  if (!any(c("watts", "energy") %in% names(d))) {
    .ckstop(fn, "`devices` needs a `watts` or an `energy` column")
  }
  num <- function(cl) {
    if (!cl %in% names(d)) {
      return(rep(NA_real_, nrow(d)))
    }
    v <- d[[cl]]
    if (!(is.numeric(v) || (is.logical(v) && all(is.na(v))))) {
      .ckstop(fn, "`devices$%s` must be numeric", cl)
    }
    v <- as.numeric(v)
    if (any(is.infinite(v) | (!is.na(v) & v < 0))) {
      .ckstop(
        fn,
        "`devices$%s` must be non-negative and finite (NA: not measured)",
        cl
      )
    }
    v
  }
  w <- num("watts")
  e <- num("energy")
  kind <- as.character(d$kind)
  list(
    kind = as.integer(tolower(trimws(kind)) %in% c("cpu", "memory")),
    watts = w,
    energy = e,
    tab = d
  )
}

.fpdevout <- function(dv, kwh) {
  d <- dv$tab
  # the first that holds: counted already, measured energy, average power
  note <- rep("no power or energy reading", length(dv$kind))
  note[!is.na(dv$watts)] <- "average power over the run"
  note[!is.na(dv$energy)] <- "measured energy"
  note[dv$kind == 1L] <- "already counted as CPU or memory energy"
  data.frame(
    kind = as.character(d$kind),
    vendor = as.character(d$vendor),
    name = as.character(d$name),
    how = as.character(d$how),
    kwh = kwh,
    counted = !is.na(kwh),
    note = note,
    stringsAsFactors = FALSE
  )
}

.fpcites <- function(ci) {
  cites <- c(
    paste0(
      "Lannelongue, Grealey & Inouye (2021) Adv Sci 8:2100707 (energy ",
      "model, 0.3725 W/GB)"
    ),
    paste0(
      "codecarbon (RAPL method; process/machine tracking; 85 W x 50% ",
      "unknown-CPU fallback)"
    ),
    sprintf("carbon intensity: %s (%s, %d)", ci$source, ci$name, ci$year)
  )
  if (is.na(ci$year)) cites <- cites[-3L]
  cites
}
