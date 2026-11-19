# Air pollution and the footprint of a computation

Dispersion of a pollutant, the health burden attributed to it, and the
energy and carbon cost of running an analysis.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Air pollution

### Atmospheric dispersion: sigmas, plume rise, plume, puff, grid and particles

[`pgsig()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
[`briggs()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
[`gplume()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md),
[`gspuff()`](https://rootcoder007.github.io/crimlab/reference/pgsig.md)
and others: Six standard models of how a pollutant released into the air
spreads downwind, for exposure estimates around a stack, a fire or an
industrial site (for example when relating air quality to health or
justice outcomes by neighbourhood):

``` r

# plume widths 1 km downwind: neutral, then very stable, urban
pgsig(1000, "D")
#> $sigmay
#> [1] 76.27701
#> 
#> $sigmaz
#> [1] 37.94733
pgsig(c(500, 1000, 2000), "F", set = "urban")
#> $sigmay
#> [1]  50.20790  92.96697 163.97832
#> 
#> $sigmaz
#> [1] 30.23716 50.59644 80.00000

# rise of a hot plume, and the effective height it gives a 50 m stack
rise <- briggs(c(100, 1000), u = 4, diam = 1, exitvel = 12,
               stktmp = 420, ambtmp = 285)
rise$rise
#> [1] 18.22335 28.88281
heff <- 50 + rise$finalrise

# ground-level concentration (g/m^3) from 100 g/s along the centreline
rec <- cbind(x = c(500, 1000, 2000, 8000), y = 0, z = 0)
gplume(100, u = 4, h = heff, rcpt = rec)
#> [1] 2.120279e-05 3.168716e-04 3.826280e-04 1.051326e-04

# under a 120 m mixing lid the far receptor sees a trapped plume
gplume(100, u = 4, h = heff, rcpt = rec, mixh = 120)
#> [1] 2.120279e-05 3.172063e-04 4.073068e-04 1.738501e-04

# a 1 kg puff, 10 minutes after release, at the puff centre
gspuff(1000, u = 3, h = 10, rcpt = rbind(c(1800, 0, 0)), t = 600)
#> [1] 0.0001266794

# a small grid: a blob moved by wind and spread by diffusion
c0 <- matrix(0, 15, 15)
c0[4, 8] <- 1
f <- advdif2d(c0, u = 1, v = 0, kx = 0.5, ky = 0.5,
              dx = 1, dy = 1, dt = 0.2, nsteps = 20)
f$stable
#> [1] TRUE
round(f$mass, 3)
#> [1] 0.962

# 2,000 simulated particles, counted on a grid
p <- lagrprt(2000, x0 = 0, y0 = 0, u = 1, v = 0, kx = 1, ky = 1,
             dt = 1, nsteps = 10, seed = 1,
             grid = c(-10, 30, -20, 20, 8, 8))
c(p$mean_x, p$meany)
#> [1] 10.033977334 -0.003237509

# inputs that cannot be right are refused
try(pgsig(100, "G"))
#> Error : pgsig: `stab` must be one of the Pasquill classes "A" to "F".
try(briggs(100, u = 4, diam = 1, exitvel = 12, stktmp = 155, ambtmp = 20))
#> Error : briggs: temperatures are in kelvin: stktmp 155 and ambtmp 20 look like degrees Celsius (ambient air is about 290 K).
```

### A Gaussian-plume concentration field on a grid

[`plmfld()`](https://rootcoder007.github.io/crimlab/reference/plmfld.md):
Evaluates `gplume` at every point of an `x` by `y` ground-level grid (or
at height `z`) and returns the field as a matrix with the axes. It maps
where the plume from a single stack touches the ground, how far downwind
the highest concentration lies, and how wide the affected area is.

``` r

f <- plmfld(100, 5, 50, x = seq(100, 3000, by = 100),
            y = seq(-500, 500, by = 50))
dim(f$field)
#> [1] 30 21
# where the ground-level maximum lies, downwind
f$x[which.max(apply(f$field, 1, max))]
#> [1] 800
# the coordinates must increase
try(plmfld(100, 5, 50, x = c(300, 100), y = c(-50, 0, 50)))
#> Error : plmfld: `x` and `y` must be strictly increasing receptor coordinates.
```

### Concentration-response functions for PM2.5 and NO2

[`crfpm25()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md),
[`crfno2()`](https://rootcoder007.github.io/crimlab/reference/crfpm25.md):
A concentration-response function (CRF) turns an air-pollution
concentration into the relative risk (RR) of a health outcome compared
with a counterfactual “clean air” level: an RR of 1.05 means 5% more
deaths (or cases) than at the reference concentration. These are the
first step of any burden estimate; `polbrdn` calls them.

``` r

# All-cause mortality at an annual PM2.5 of 12 micrograms per cubic metre
r <- crfpm25(12)
c(rr = r$rr, lower = r$rrlo, upper = r$rrup)
#>       rr    lower    upper 
#> 1.048873 1.036787 1.054883
# the cause-specific IER curve for ischaemic heart disease
crfpm25(12, y = "ihd")$rr
#> [1] 1.553634
# NO2 and respiratory deaths; at or below the reference no excess risk
crfno2(25, y = "respiratory")$rr
#> [1] 1.045336
crfno2(8)$rr
#> [1] 1
# a vector of exposures (one per area): the mean RR and each area's RR
crfpm25(c(5, 12, 20))$extra$rrperunit
#> [1] 1.000000 1.048873 1.115480
# a coefficient of the reader's own: RR 1.04 per 10 micrograms
crfno2(30, betaper10 = log(1.04))$rr
#> [1] 1.0816
# the environment as a driver of offending: a simulated daily NO2 series
# and simulated daily counts of assaults for three areas, joined by area
# and day
set.seed(1)
d <- data.frame(area = rep(c("a", "b", "c"), each = 30), day = rep(1:30, 3))
d$no2 <- stats::rlnorm(90, log(20), 0.3)
d$assaults <- stats::rpois(90, 3 + 0.05 * d$no2)
# the relative risk each area-day's exposure carries, and the mean count
# on days above and below its median
d$rr <- crfno2(d$no2)$extra$rrperunit
tapply(d$assaults, d$rr > stats::median(d$rr), mean)
#>    FALSE     TRUE 
#> 3.822222 4.222222
# an outcome without a published function is refused
try(crfpm25(12, y = "asthma"))
#> Error : crfpm25: Unknown outcome 'asthma'. Available: allcause, ihd, stroke.
```

### Population attributable fraction

[`attrfrac()`](https://rootcoder007.github.io/crimlab/reference/attrfrac.md):
The share of all cases in a population that is due to an exposure: the
fraction that would not occur if the exposure were removed and nothing
else changed. It combines how strongly the exposure raises risk (the
relative risk) with how common it is (the prevalence); multiplied by the
number of cases it gives the attributable cases, as `polbrdn` does.

``` r

attrfrac(1.5, 0.4)
#> Attributable fraction 0.1667 (relative risk 1.5, 0.4 exposed)
attrfrac(2, 1) # everyone exposed: 1 - 1/RR
#> Attributable fraction 0.5 (relative risk 2, 1 exposed)
paf <- attrfrac(1.08, 1, rrlo = 1.06, rrup = 1.09)
attr(paf, "ci")
#>      lower      upper 
#> 0.05660377 0.08256881
attrfrac(0.8, 0.5) # protective
#> Attributable fraction -0.1111 (relative risk 0.8, 0.5 exposed)
#> Negative: the exposure protects.
try(attrfrac(1.5, 40))
#> Error : attrfrac: expoprev must be in [0,1]; got 40.
```

### Deaths avoided under a counterfactual exposure reduction

[`mortdisp()`](https://rootcoder007.github.io/crimlab/reference/mortdisp.md):
The expected number of deaths (or cases) a year that a given drop in
exposure would avoid in a population, from the baseline rate and a
log-linear concentration-response coefficient. This is the health-impact
function used in regulatory benefit assessments.

``` r

# a 10-unit cut, a million people, 800 deaths per 100,000 a year
mortdisp(10, 1e6, 0.008, 0.0039)
#> [1] 305.9943
# the PM2.5 all-cause coefficient is published per 10 micrograms
mortdisp(5, 2.8e6, 0.0075, log(1.08) / 10)
#> [1] 792.7406
# a rare outcome must be declared
mortdisp(5, 1e6, 5e-6, 0.01, rare = TRUE)
#> [1] 0.2438529
# a rate per 100,000 passed as a rate per person is refused
try(mortdisp(10, 1e6, 800, 0.0039))
#> Error : mortdisp: `brate` is a rate per person-year, between 0 and 1 (0.008 is 800 per 100,000); 800 looks like a rate per 100,000.
```

### End-to-end pollution burden

[`polbrdn()`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md):
Chains the concentration-response function, the attributable fraction
and the baseline case count into annual attributable cases, as the GBD
risk factor estimates do (GBD 2019 Risk Factors Collaborators 2020). It
serves a single population (a city, a region) whose mean exposure,
baseline rate and size are known; for many areas at once see
`polbrdarea`.

``` r

# ambient NO2 of 25 micrograms per cubic metre over a million people
b <- polbrdn(25, 1, 0.008, 1e6, poll = "NO2")
b$attrcases
#> [1] 234.1369
# PM2.5 and ischaemic heart disease deaths, against the WHO guideline
polbrdn(12, 1, 0.0012, 2.8e6, poll = "PM2.5", y = "ihd",
        refconc = 5)$attrcases
#> [1] 1235.817
# childhood asthma counts incident cases, not deaths
polbrdn(30, 1, 0.012, 5e5, y = "chiast")$extra$unit
#> [1] "incident cases"
try(polbrdn(40, 1, 0.008, 1e6, poll = "O3"))
#> Error : polbrdn: Unknown pollutant 'O3'. Use 'PM2.5' or 'NO2'.
```

### Per-area pollution burden

[`polbrdarea()`](https://rootcoder007.github.io/crimlab/reference/polbrdarea.md):
Applies `polbrdn` to every row of a table with one row per area (forward
sortation area, census tract, neighbourhood, …) and returns the rows
sorted by attributable cases, worst first: where in a city or region the
pollution burden falls.

``` r

# simulated areas
t <- data.frame(fsa = c("A1", "B2"), exposure = c(28, 18),
                population = c(40000, 60000), brate = 0.008)
polbrdarea(t)
#> Pollution burden by area: 2 areas, 18.75 attributable of 800 cases (2.34%)
#>  fsa exposure population brate    rr     paf attrcases basecases
#>   A1       28      40000 0.008 1.036 0.03502    11.205       320
#>   B2       18      60000 0.008 1.016 0.01572     7.544       480
try(polbrdarea(t[, c("fsa", "exposure")]))
#> Error : polbrdarea: areatab missing columns: population, brate.
```

### Concentration index of exposure by income

[`expconcidx()`](https://rootcoder007.github.io/crimlab/reference/expconcidx.md):
Measures whether exposure to pollution is concentrated among
lower-income or higher-income people. The index runs from -1 (all the
exposure borne by the poorest) through 0 (exposure shared evenly across
income) to +1 (all of it borne by the richest); environmental-justice
analyses report it next to the burden.

``` r

d <- data.frame(exposure = c(30, 25, 20, 15, 10), income = 1:5)
ci <- expconcidx(d, "exposure", "income")
ci$concidx
#> [1] -0.2
ci$interpretation
#> [1] "CI = -0.2000. Pro-poor exposure burden: lower-income individuals bear disproportionately higher pollution."
try(expconcidx(data.frame(exposure = 1:3, income = c("a", "b", "c")),
               "exposure", "income"))
#> Error : expconcidx: `income` must be numeric, got character.
```

### Exposure-response estimate with a bootstrap sensitivity interval

[`exprspplr()`](https://rootcoder007.github.io/crimlab/reference/exprspplr.md):
Fits the partially linear exposure-response model (the effect of an
exposure on an outcome, the confounders partialled out by cross-fitted
regressions), or any estimator with the same interface, then re-fits it
on `nboot` non-parametric resamples of the rows for a percentile
interval (Efron and Tibshirani 1993, chapter 13). The bootstrap interval
does not lean on the estimator’s own standard error, so it is a check on
it: when the two disagree, the analytic one should not be trusted alone.

``` r

# simulated areas: exposure driven by x1, the outcome by the exposure and
# x2
set.seed(1)
d <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
d$exposure <- 20 + 2 * d$x1 + rnorm(200)
d$asthma <- 5 + 0.3 * d$exposure + d$x2 + rnorm(200)
r <- exprspplr(d, y = "asthma", expo = "exposure", confs = c("x1", "x2"),
               nboot = 20)
r$ate
#> [1] 0.2536464
c(r$cilobs, r$ciupbs)
#> [1] 0.1365177 0.3476568

# Any estimator with the same interface: here an ordinary regression.
ols <- function(data, outcome, treatment, covariates) {
  f <- stats::reformulate(c(treatment, covariates), outcome)
  cf <- summary(stats::lm(f, data = data))$coefficients
  list(ate = cf[treatment, 1], se = cf[treatment, 2], method = "OLS")
}
exprspplr(d, "asthma", "exposure", c("x1", "x2"), nboot = 20,
          estr = ols)$method
#> [1] "OLS + percentile bootstrap"

# Fewer than 10 resamples cannot give an interval.
try(exprspplr(d, "asthma", "exposure", c("x1", "x2"), nboot = 5))
#> Error : exprspplr: Only 5 successful bootstrap fits; need >=10. Inspect the model's convergence on your resamples.
```

### Run the pollution-to-health pipeline and report it

[`vrfpol()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md),
[`polrpttxt()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md):
Resolves the exposure (simulated demo data, a CSV with an `exposure`
column and an optional `income` column, a NAPS-style pull with a `value`
column and a `unit`, or scalar arguments), logs the assumptions, and
when they all hold runs the concentration-response,
attributable-fraction, displaced-mortality, burden and equity stages.
[`polrpttxt()`](https://rootcoder007.github.io/crimlab/reference/vrfpol.md)
renders the result as plain text. It is the one-call version of the air
pollution functions, made for scripts: for bad data or a value it cannot
use it does not stop with an R error, but returns a report whose status
says what went wrong. Only `poll`, which names the report, is checked
first and refused with an R error (as the command line refuses it) when
it is not one pollutant name: missing, not text, `NULL`, several names,
or text that is not valid UTF-8.

``` r

r <- vrfpol("no2", demo = TRUE)
r$status
#> [1] "ok"
r$pipeline$paf
#> [1] 0.02733115
cat(polrpttxt(r))
#> ==================================================================
#>   verify-pollution -- NO2 -> allcause
#>   data:   demo (synthetic)
#> ==================================================================
#> 
#> Inputs
#>   exposure mean:      24.261
#>   exposure prevalence:0.981
#>   baseline rate/100k: 500.00
#>   population:         1,000,000
#>   reference conc:     10
#> 
#> Assumption log
#>   [PASS] exposure finite and non-negative -- mean 24.26113743 vs ref 10
#>   [PASS] exposure within the concentration-response range (at most 10,000 ug/m3) -- exposure_mean=24.26114
#>   [PASS] prevalence in [0,1] -- exposure_prevalence=0.981
#>   [PASS] brate finite: zero, or between 1 and 10,000 per 100k -- brate=500 per 100k per year
#>   [PASS] population finite, positive and at most 10 billion -- population=1,000,000
#>   [PASS] reference finite and non-negative -- reference=10
#>   [PASS] pollutant supported by the CRFs -- Current CRFs: NO2 (log-linear), PM2.5 (log-linear all-cause; Burnett IER for IHD and stroke). Other pollutants reject.
#> 
#> Concentration-response
#>   RR:       1.0286
#>   source:   Huangfu & Atkinson (2020) Environ Int 144:105998; WHO (2021) Global AQ Guidelines
#> 
#> Attributable fraction (PAF): 0.0273  (95% CI 0.0138 to 0.0534, from the risk's limits)
#> 
#> Mortality displaced
#>   expected avoided deaths: 136.6
#> 
#> Burden of pollution
#>   attributable deaths:   136.7
#> 
#> Equity analysis
#>   concentration index: -0.0647
#> 
#> STATUS: ok

# Scalar inputs: a city of 2.8 million at a PM2.5 mean of 9.
s <- vrfpol("pm25", expomean = 9, expoprev = 1, brate = 750, pop = 2.8e6,
            reg = "Toronto")
s$pipeline$burden$attrcases
#> [1] 510.8617
attr(s, "exitstat")
#> [1] 0

# A CSV of simulated measurements with an income column.
f <- tempfile(fileext = ".csv")
utils::write.csv(data.frame(exposure = c(12, 18, 25, 31, 9), income = 1:5),
                 f, row.names = FALSE)
vrfpol("no2", expcsv = f)$pipeline$equity$concidx
#> [1] 0.02947368

# A rate per person-year passed as per 100,000 fails an assumption ...
bad <- vrfpol("no2", expomean = 20, expoprev = 1, brate = 0.008)
bad$status
#> [1] "assumption_failure"
# ... and an argument it cannot use is an error report, not an R error.
vrfpol("no2")$error
#> [1] "no exposure data: pass expcsv, or expomean and expoprev, or demo = TRUE"
```

## The footprint of a computation

### Energy and carbon footprint of a computation

[`co2fp()`](https://rootcoder007.github.io/crimlab/reference/co2fp.md):
Evaluates `expr` and reports the electricity it used and the
CO2-equivalent that electricity implies, by one of two methods, for the
grid it ran on, offline. Use it to report the energy and emissions of an
analysis next to its results (as journals and funders increasingly ask),
or to compare two ways of running the same computation. The result keeps
the value of `expr`, and lists every figure that was assumed rather than
measured.

``` r

fp <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON")
fp$enekwh$total
#> [1] 4.026175e-08
fp$co2eg
#> [1] 3.662611e-06
fp$co2int$year
#> [1] 2024
fp$assumptions
#> [1] "CPU power unknown: 85 W x 50% (codecarbon fallback)"
#> [2] "memory charged: the machine's total 15.6 GB"
print(fp)
#> Computation footprint (modelled)
#>   Time: 0.003 s wall, 0.003 s CPU; utilisation 0.25 (process) on 4 cores
#>   Energy: 4.03e-08 kWh (cpu 3.54e-08, memory 4.85e-09)
#>   CO2e: 3.66e-06 g at 90.97 gCO2e/kWh (Canada, Ontario, 2024)
#>   Assumed: CPU power unknown: 85 W x 50% (codecarbon fallback); memory charged: the machine's total 15.6 GB
# wherever you are, offline: the location is detected and recorded
co2fp(1 + 1, cpupwrw = 45, memgb = 8)$assumptions
#> [1] "location WORLD via fallback, no location found in environment, time zone or locale"
# the whole machine's utilisation, on codecarbon's machine curve
co2fp(
  1 + 1,
  loc = "FR",
  usage = "machine",
  loadcur = "codecarbon",
  cpupwrw = 45,
  memgb = 8
)$usagemode
#> [1] "machine"
# a fixed utilisation and a data-centre PUE, with a carbon intensity given
# as a number
fp2 <- co2fp(
  sum(sqrt(seq_len(2e5))),
  loc = 30,
  usage = 0.5,
  cpupwrw = 65,
  memgb = 16,
  pue = 1.67,
  mthd = "model"
)
fp2$enekwh
#> $cpu
#> [1] 4.522917e-08
#> 
#> $memory
#> [1] 8.294333e-09
#> 
#> $total
#> [1] 5.35235e-08
fpeqv(fp2$co2eg)
#> Footprint in everyday terms
#>   Car travel: 9.175e-09 km
#>   Tree uptake: 1.751e-09 tree-months
# a graphics card measured at 120 W on average (a simulated reading) and
# a CPU row, which is not counted twice
dev <- data.frame(
  kind = c("gpu", "cpu"), vendor = c("example", "example"),
  name = c("card 0", "package 0"), watts = c(120, 35),
  how = c("simulated reading", "simulated reading")
)
fp3 <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON", cpupwrw = 65,
             memgb = 16, mthd = "model", devices = dev)
fp3$devices
#>   kind  vendor      name               how   kwh counted
#> 1  gpu example    card 0 simulated reading 1e-07    TRUE
#> 2  cpu example package 0 simulated reading   NaN   FALSE
#>                                      note
#> 1              average power over the run
#> 2 already counted as CPU or memory energy
# a bad argument is refused before the expression runs
try(co2fp(Sys.sleep(10), pue = 0.9))
#> Error : co2fp: `pue` must be a single number >= 1.
```

### Is CPU energy measurable here?

[`raplok()`](https://rootcoder007.github.io/crimlab/reference/raplok.md):
`TRUE` when the operating system exposes readable RAPL (Running Average
Power Limit) energy counters for at least one CPU package: Linux with
the `intel_rapl` powercap driver (Intel and, since kernel 5.8, AMD),
with read permission on `energy_uj`. On most systems since Linux 5.10
that file is readable by root only, so `FALSE` is the common answer for
an ordinary user, and `co2fp` then falls back to its model.

``` r

raplok()
#> [1] FALSE

# a fake powercap tree with one readable package domain (Linux names the
# domain "intel-rapl:0"; a colon is not allowed in a Windows file name)
if (.Platform$OS.type != "windows") {
  d <- file.path(tempfile(), "intel-rapl:0")
  dir.create(d, recursive = TRUE)
  writeLines("package-0", file.path(d, "name"))
  writeLines("123456789", file.path(d, "energy_uj"))
  raplok(dirname(d))
}
#> [1] TRUE

# a folder that does not exist is simply FALSE; a non-path is an error
raplok(file.path(tempdir(), "no-such-dir"))
#> [1] FALSE
try(raplok(NA_character_))
#> Error : raplok: `rapldir` must be a single directory path.
```

### Carbon intensity of electricity by location, offline, anywhere

[`co2int()`](https://rootcoder007.github.io/crimlab/reference/co2int.md),
[`co2inttab()`](https://rootcoder007.github.io/crimlab/reference/co2int.md):
The grams of CO2-equivalent emitted per kilowatt-hour of electricity at
a location, with the year and source of the figure. The table is
compiled into the package, so it works offline in any country: every
country at the latest year published by Our World in Data from Ember’s
yearly electricity data (the latest year each grid reports: 2025 or 2024
for nearly all, older for a few grids, Ukraine (2022) among them, named
in `year`), sub-national zones from Electricity Maps’ 2024 yearly data
as redistributed in the Green Algorithms data release v3.1 (Canadian
provinces and territories, US balancing authorities, Australian states,
Indian and Japanese regions, Brazilian and Chilean systems, …), and the
world average. All figures are lifecycle intensities (generation plus
upstream), the basis both sources publish.

``` r

co2int("CA-ON")
#> Carbon intensity of electricity: 90.97 gCO2e/kWh
#>   Location: CA-ON (Canada, Ontario), 2024
#>   Basis: lifecycle
#>   Source: Electricity Maps (2025). Canada 2024 Yearly Carbon Intensity Data (Version January 27, 2025). Electricity Maps (via Green Algorithms data v3.1)
co2int("Canada")$year
#> [1] 2025
co2int("CA-XX")$note
#> [1] "no row for zone CA-XX; its country CA used"
co2int(120)$source
#> [1] "user-supplied carbon intensity"
nrow(co2inttab())
#> [1] 378
# the Canadian provinces in the table, lowest intensity first
tab <- co2inttab()
ca <- tab[
  startsWith(tab$location, "CA-"),
  c("location", "name", "g_per_kwh", "year")
]
head(ca[order(ca$g_per_kwh), ])
#>     location                              name g_per_kwh year
#> 236    CA-NL Canada, Newfoundland and Labrador     27.60 2024
#> 227    CA-QC                    Canada, Quebec     35.04 2024
#> 322    CA-MB                  Canada, Manitoba     52.89 2024
#> 367    CA-BC          Canada, British Columbia     70.34 2024
#> 241    CA-PE      Canada, Prince Edward Island     70.94 2024
#> 313    CA-ON                   Canada, Ontario     90.97 2024
# a location the table does not know is an error, not a silent world
# average
try(co2int("Atlantis"))
#> Error : co2int: Unknown location 'Atlantis'. Use a code, alpha-3 code or name from co2inttab(), or pass a number in gCO2e/kWh.
```

### Where is this computation running?

[`dtctloc()`](https://rootcoder007.github.io/crimlab/reference/dtctloc.md):
Finds a location code for `co2int` without any network access, in this
order:

``` r

dtctloc()
#> $location
#> [1] "WORLD"
#> 
#> $method
#> [1] "fallback, no location found in environment, time zone or locale"
dtctloc(envloc = "", codloc = "", tz = "Europe/Paris")
#> $location
#> [1] "FR"
#> 
#> $method
#> [1] "system time zone Europe/Paris"
dtctloc(envloc = "", codloc = "", tz = "", locale = "en_IN.UTF-8")
#> $location
#> [1] "IN"
#> 
#> $method
#> [1] "locale en_IN.UTF-8"
# an explicit setting wins over everything else
dtctloc(envloc = "CA-QC", tz = "Europe/Paris")
#> $location
#> [1] "CA-QC"
#> 
#> $method
#> [1] "CRIMLAB_LOCATION environment variable"
# nothing usable: the world average, and the method says so
dtctloc(envloc = "", codloc = "", tz = "UTC", locale = "C")
#> $location
#> [1] "WORLD"
#> 
#> $method
#> [1] "fallback, no location found in environment, time zone or locale"
```

### Everyday equivalents of a quantity of CO2-equivalent

[`fpeqv()`](https://rootcoder007.github.io/crimlab/reference/fpeqv.md):
Converts grams of CO2-equivalent into two everyday comparisons:
kilometres driven in an average European passenger car (175 gCO2/km) and
tree-months of carbon sequestration (917 gCO2 per tree per month), the
two comparisons the Green Algorithms calculator prints, with their
sources. Use it to make the `co2eg` of `co2fp` easier to read in a
report.

``` r

fpeqv(1200)
#> Footprint in everyday terms
#>   Car travel: 6.857 km
#>   Tree uptake: 1.309 tree-months

# the footprint of a small computation, in everyday terms
fp <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON", cpupwrw = 45,
            memgb = 8)
unlist(fpeqv(fp$co2eg)[c("carkm", "treemon")])
#>        carkm      treemon 
#> 6.594603e-09 1.258512e-09

# a vector is refused: convert one total at a time
try(fpeqv(c(10, 20)))
#> Error : fpeqv: `co2eg` must be a single non-negative number.
```
