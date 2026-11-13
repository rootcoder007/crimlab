# crimlab <img src="man/figures/logo.png" align="right" height="139" alt="crimlab hex logo" />

[![check](https://github.com/rootcoder007/crimlab/actions/workflows/check.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/check.yml)
[![coverage](https://github.com/rootcoder007/crimlab/actions/workflows/coverage.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/coverage.yml)
[![srr](https://github.com/rootcoder007/crimlab/actions/workflows/srr.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/srr.yml)
[![pkgdown](https://github.com/rootcoder007/crimlab/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/rootcoder007/crimlab/actions/workflows/pkgdown.yml)

Statistical methods for courts, corrections, policing, forensics and
sociolegal research. The computations run in a compiled C++17 core with a
plain C interface (`src/crimlab.h`); the R functions check their arguments
and call it. Where an identity or a bound has a machine-checked proof in
Lean, `proof()` names the theorem.

![crimlab in one minute: a veil-of-darkness table, a reported odds ratio checked against every labelling of its table, a clearance paradox, a categorical audit, a posterior from a likelihood ratio, hot spots, spatial autocorrelation, the command-line program, and the theorems behind a function](man/figures/demo.svg)

It is written for researchers and analysts in criminology, sociolegal
studies, policing, corrections and forensic science who need these
quantities computed correctly and want to know what each one does and does
not show. `vignette("crimlab")` walks through one example from each family.

## What it takes

The functions work on summary data, not on raw records:

* counts by group: stops, searches, uses of force, offences recorded and
  cleared, residents, offenders and victims;
* individual outcomes coded 0 or 1, with a group or a decision for each
  person: release, sentence, reconviction;
* allele frequencies and genotypes at one or more loci;
* likelihood ratios from validation comparisons whose source is known;
* networks as edge lists of labels (who offended with whom, who passed
  what to whom), or incident tables that list the offenders of each
  incident.
* point locations of incidents, planar or longitude and latitude (vectors,
  tables or sf points), with times, weights or counts by area;
* risk scores with outcomes and groups, trace-evidence measurements,
  event time windows, reported figures, matched-pair differences and
  counts over time.

## What it asks

* Disparities: does visibility change who is stopped, what an odds ratio
  says about risk, whether hit rates can show bias, how a disparity per
  resident splits into contact and force.
* Forensics: how rare a DNA profile is, what a likelihood ratio is worth
  and what it is not, how well a system's ratios are calibrated, kinship,
  paternity, mixtures and database searches.
* Policing: why an aggregate clearance rate can fall while every offence
  type clears better, what a proxy for exposure allows, what
  misrecorded group membership does, regression to the mean at hot spots.
* Corrections: what incapacitation prevents, what a release rule's
  failure rate can be when only the released are seen, mean time in
  custody from current stays, and what observational data alone say about
  a sentence's effect.
* Networks: who co-offends with whom, who sits between groups, which
  offenders form the core and the communities of a network, and how
  cleanly it splits in two.
* Places: where incidents centre and how far they spread, whether they
  cluster more than chance, where the hot spots are, whether places close
  together have similar rates (Geary's C, the Getis-Ord G, local Moran),
  whether incidents close in space are close in time, and where to search
  first for the base of a series.
* Methods from the forensic and sociolegal literature: the fairness of a
  risk tool by group, the weight of glass and other trace evidence, when
  an event known only within a window happened, whether reported figures
  follow Benford's law, how strong hidden bias must be to explain a
  matched result away, and group-based trajectories of offending.
* Counts and rates: rate ratios between groups from Poisson and negative
  binomial regression with an exposure, adjustment and cluster-robust or
  GEE standard errors (`cntreg()`); rates and proportions by group with
  their ratios and differences, exact and approximate intervals and
  over-representation (`grprr()`); directly standardised rates
  (`dastd()`); observed over expected counts with funnel limits and
  empirical Bayes shrinkage (`sirex()`).
* Tables: cross-tabulation of records or of an aggregate table's counts,
  with chi-square and likelihood-ratio tests, Fisher's exact test,
  Cramer's V with its interval, and the share on the diagonal (`ctabx()`).
* Presets: how the columns of a published dataset (Ontario segregation and
  deaths-in-custody tables, use-of-force files, Special Investigations
  Unit cases, Structured Intervention Unit tables, Toronto incident
  layers) feed these functions, so that `cntreg(clprep(data,
  "otis-c01-cnt"))` gives the published analysis (`clprep()`). Presets
  describe layouts and ship no data.
* Checks: which labelling of a table reproduces a reported odds ratio,
  risk ratio or risk difference (`relchk()`), and what the functions here
  would make of categorical inputs before an analysis (`catchk()`).

## Algorithms

Every method is a published one, cited on its help page; none is new, and
the package does not claim to be the first implementation of any of them in
R. What it adds is a single checked implementation: arguments are checked
in R and again in the compiled kernels, quantities that overflow or cancel
in the obvious formula are computed in a stable form, and most functions
rest on identities or bounds proved in Lean, which `proof()` names.

## Speed

`dev/bench.R` times each function against every package found that
computes the same quantity, on fixed seeded inputs and at sizes met in
practice, after checking that the two give the same values; where no
package does, against a base-R implementation of the formula. The table it
writes, with the machine and the date, is
[inst/bench/bench.md](inst/bench/bench.md) (and `bench.csv`); the bench
workflow runs it again every week. The lists below are written by the same
script from that table.

<!-- bench start: written by dev/bench.R -->
Faster, by at least a tenth, than every package that computes the same quantity, at every size timed: `vodtst()`, `orrisk()`, `outtst()`, `dspbm()`, `dyadrt()`, `gtmp()`, `kinidx()`, `patidx()`, `mixlr()`, `cllr()`, `cllrmin()`, `lrcal()`, `eceplt()`, `expbnd()`, `misrd()`, `rtmhs()`, `incapr()`, `sellab()`, `losbnd()`, `sntbnd()`, `hzdsel()`.

Not faster than such a package at some size timed: `incap()`.

No package found that computes the same quantity, so timed against a base-R implementation of the formula: `lrpost()`, `dbsrch()`, `tipplt()`, `clrsmp()`, `proof()`; of these, not faster than base R at some size: `lrpost()`, `dbsrch()`.
<!-- bench end -->

## Life cycle

Experimental. This is a development version: function names, arguments and
results may change before a first release.

## Installation

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("rootcoder007/crimlab")
```

The package needs a C++17 compiler and no R package beyond those that ship
with R.

## Examples

### Disparities in policing

```r
library(crimlab)

# veil of darkness: does visibility change who is stopped?
vodtst(day = c(white = 600, black = 400), night = c(white = 700, black = 300))

# the same counts from a data frame: groups are matched by name, never by
# position, and a table with day and night columns can be given whole
stops <- data.frame(race = c("white", "black", "white", "black"),
                    light = c("day", "day", "night", "night"),
                    n = c(600, 400, 700, 300))
vodtst(xtabs(n ~ race + light, stops), ref = "white")

# the risks, and the relative risk, behind an odds ratio
orrisk(or = 3, rate = 0.4, share = 0.5)

# outcome test: stop and hit rates at a search threshold
outtst(q = c(0.1, 0.5, 0.9), t = 0.3)

# force per resident = contact per resident x force per contact
dspbm(pop = c(a = 9000, b = 1000), contact = c(a = 900, b = 300),
      force = c(a = 9, b = 6))

# offending rates per offender and per offender-victim pair
pop <- c(a = 800, b = 200)
dyadrt(count = outer(pop, pop) / 1000, pop = pop)
```

### Clearance rates

```r
# both offence types clear better in B, yet B's aggregate rate is lower
clrsmp(clra = c(minor = 90, serious = 20), reca = c(minor = 100, serious = 50),
       clrb = c(minor = 19, serious = 60), recb = c(minor = 20, serious = 140))
```

### Corrections

```r
# custody against a community sentence, with no assumption about who got which
y <- c(1, 0, 1, 1, 0, 0, 1, 0)
z <- c("community", "community", "custody", "custody", "custody",
       "community", "custody", "community")
sntbnd(y, z, treat = "custody")

# a release rule scored from the outcomes of those the judge released
sellab(y = c(0, 1, 0, NA, NA), rel = c(TRUE, TRUE, TRUE, FALSE, FALSE),
       rule = c(TRUE, FALSE, TRUE, TRUE, FALSE))
```

### DNA match statistics

```r
# match probability with the Balding-Nichols correction
gtmp(p = c(D3 = 0.1, vWA = 0.2), q = c(0.3, NA), hom = c(FALSE, TRUE),
     theta = 0.03)

# a likelihood ratio needs a prior before it is a probability
lrpost(lr = 1e6, prior = 1 / 5e6)

# a match found by searching a database
dbsrch(mp = 1e-9, n = 5e6, pop = 3e7)
```

### Kinship, paternity and mixtures

```r
kinidx(g1 = c(12, 14), g2 = c(12, 14), freq = c(0.1, 0.2), rel = "sibling")
patidx(child = c(12, 14), father = c(14, 15), freq = c(0.1, 0.2),
       mother = c(12, 12))

p <- c(a = 0.1, b = 0.2, c = 0.05, d = 0.15)
mixlr(p, np = character(0), xp = 0, nd = c("c", "d"), xd = 1)
```

### Performance of likelihood ratios

```r
same <- c(50, 200, 8, 0.5)
diff <- c(0.01, 0.2, 1.5, 0.001)
cllr(same, diff)
cllrmin(same, diff)
lrcal(same, diff)
tipplt(same, diff)
eceplt(same, diff)
```

### Co-offending networks

```r
# offenders recorded together in an incident are tied, weighted by the
# number of incidents they share
inc <- c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5)
off <- c("ana", "ben", "cy", "ana", "ben", "cy", "dee", "dee", "eli", "fay")
g <- coofg(inc, off)
e <- g$edges
ntdeg(e$from, e$to, w = e$w, nodes = g$nodes)
ntbtw(e$from, e$to, nodes = g$nodes)
ntlvn(e$from, e$to, w = e$w, nodes = g$nodes)
ntbrk(e$from, e$to, grp = c(ana = "x", ben = "x", cy = "x", dee = "y",
                            eli = "y", fay = "z"), nodes = g$nodes)
```

### Spatial input

```r
# planar coordinates (metres, say): two vectors, a matrix or a data frame
# with columns x and y; longitude and latitude only with lonlat = TRUE
set.seed(1)
p <- data.frame(x = rnorm(100, 500, 50), y = rnorm(100, 800, 30))
cgstat(p)
sdell(p, c = 2)

ll <- data.frame(lon = c(-79.40, -79.38, -79.35), lat = c(43.64, 43.66, 43.65))
webpt(ll, lonlat = TRUE)          # the centre comes back in degrees
hvdist(ll, c(-79.38, 43.65))      # great-circle metres
eqproj(ll)                        # equal-area metres for the planar functions

# sf points go in directly; their CRS must agree with lonlat
s <- sf::st_as_sf(ll, coords = c("lon", "lat"), crs = 4326)
cgstat(s, lonlat = TRUE)
# and plain results go back to sf
sf::st_as_sf(eqproj(ll), coords = c("x", "y"))
```

The planar functions assume Cartesian coordinates in one length unit and
return lengths in that unit; functions with a `lonlat` argument project
longitude and latitude to metres themselves (an equal-area projection about
the points' centre) and return centres in degrees; functions that take a
study rectangle (`win`) take planar coordinates only. All are
two-dimensional.

### Hot spots, space-time and journey to crime

```r
set.seed(2)
x <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
y <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
fzmode(x, y, r = 0.5)                               # points within 0.5 of each
nnhclu(x, y, area = 225, p = 0.05, minpts = 5)$clusters
kmpart(x, y, k = 3, seed = 1)$cen
stacel(x, y, win = c(-5, 10, -5, 10), r = 0.75, top = 10, minpts = 5)$clusters
mantst(x, y, t = seq_along(x), nsim = 199, seed = 1)
mantst(x, y, t = seq_along(x), nsim = 199, seed = 1, metric = "manhattan")
hdbang(x, y, v = fzmode(x, y, r = 0.5))      # head-bang smoothing
crwalk(cbind(x = x[1:10], y = y[1:10]), h = 2, seed = 1)$msd
# clusters relative to a population at risk; zones as counts at centroids
rnnhcl(x, y, base = cbind(x = runif(500, -5, 10), y = runif(500, -5, 10)),
       h = 1, minpts = 5)$clusters
fzmode(x[1:20], y[1:20], r = 1, count = rpois(20, 3))

# a search priority over a grid, not evidence about any person
s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30, ny = 30,
            type = "negexp", par = c(1, 1), mode = "lik")
jtacc(s$prob, win = c(0, 6, 0, 6), base = c(3.2, 2.9))
```


### Spatial weights and autocorrelation

```r
set.seed(4)
p <- cbind(x = runif(60), y = runif(60))
v <- ifelse(p[, 1] > 0.7 & p[, 2] > 0.7, 5, 1) + rnorm(60, sd = 0.3)
w <- spwt(p, k = 5)                    # or type = "band", or "contig"
gearc(v, w)                            # Geary's C
gorg(v, spwt(p, type = "band", d = 0.2, style = "raw"))   # global G
head(lmorn(v, w, nsim = 199, seed = 1))                    # local Moran
locgi(v, spwt(p, type = "band", d = 0.2, style = "raw"))  # local G*
dirmn(p, origin = c(0.5, 0.5))         # directional mean from an origin
```
### Risk tools, trace evidence, time and data screens

```r
set.seed(3)
g <- rep(c("a", "b"), each = 100)
y <- rbinom(200, 1, ifelse(g == "a", 0.3, 0.5))
s <- plogis(rnorm(200, ifelse(y == 1, 0.8, -0.8)))
fair(s, y, g)$met                     # error rates and predictive values by group

rswp(rnorm(50, 0.5), gam = c(1, 1.5, 2))   # how much hidden bias explains it away
bnft(bdig(exp(runif(500, 0, 10)))$dig)$stat
gstr(matrix(rpois(100, 2), 10))[1:3, ]
sptst(rpois(30, 10), rpois(30, 10), nboot = 199, seed = 1)$glob
### Checking reports and inputs

```r
# which reference group and event give a reported odds ratio of 0.44?
tab <- matrix(c(20, 10, 80, 90), 2,
              dimnames = list(c("A", "B"), c("yes", "no")))
relchk(tab, 0.44, ref = "B", event = "yes")

# missing values, unused and rare levels, labels that differ by case, and
# the reference level a default would take
d <- data.frame(sex = c("Male", "male", "Female", "Female", NA),
                custody = c(TRUE, FALSE, TRUE, TRUE, FALSE))
catchk(d, nmin = 2)
```

### Proofs

```r
proof("vodtst")
```

## Checking the proofs

The Lean 4 files are in `inst/lean` (installed under
`system.file("lean", package = "crimlab")`). With
[elan](https://github.com/leanprover/elan) installed:

```sh
cd inst/lean
lake exe cache get   # prebuilt Mathlib, about 5 GB
lake build
```

A build with no error and no `sorry` is the check.

## Citation

`citation("crimlab")` gives the reference for the installed version:

> Ruhela, V. S. (2026). *crimlab: Criminology and Sociolegal Statistics.*
> R package. <https://github.com/rootcoder007/crimlab>

## Licence

Apache License 2.0 with LLVM Exceptions; see `LICENSE`. The SPDX identifier is
`Apache-2.0 WITH LLVM-exception`. The exception lets compiled code that embeds
the kernels be redistributed without the attribution conditions, and settles
linking with GPL-2 software.

## Python

`py/` holds a Python package, also named crimlab, that compiles the same
kernels from `src/` and gives the same functions, argument names and
results; its tests check every number against the R package. It has no
runtime dependencies. Install it from a clone with `pip install ./py`; see
`py/README.md` for how inputs are matched by name and the few differences
from R.
