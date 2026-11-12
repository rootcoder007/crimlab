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
* likelihood ratios from validation comparisons whose source is known.

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

## Algorithms

Every method is a published one, cited on its help page; none is new, and
the package does not claim to be the first implementation of any of them in
R. What it adds is a single checked implementation: arguments are checked
in R and again in the compiled kernels, quantities that overflow or cancel
in the obvious formula are computed in a stable form, and most functions
rest on identities or bounds proved in Lean, which `proof()` names.

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
