# crimlab

Criminology and Sociolegal Statistics: statistical methods for courts,
corrections, policing, forensics and sociolegal research. The computations
run in a compiled C++17 core with a plain C interface (`src/crimlab.h`); the R
functions check their arguments and call it. Where an identity or a bound
has a machine-checked proof in Lean, `proof()` names the theorem.

This is a development version.

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
sntbnd(y, z)

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

## Licence

AGPL-3. The C++ core (`src/`) carries the SPDX identifier
AGPL-3.0-or-later.
