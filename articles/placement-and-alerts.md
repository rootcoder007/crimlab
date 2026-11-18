# Placement volatility and alert complexity

How often a person in custody is moved within a year, how many alerts a
record carries, and what follows from both.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Placement volatility of each person-year

[`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md):
Counts how often a person’s placement region changes within a period (a
fiscal year), from records that each carry the region at the time of
placement `ra` and the region of the most recent placement `rb`. For a
person-year with records r = 1, …, R in a given order the count is

``` r

ra <- c("East", "East", "West", "North", "North")
rb <- c("East", "West", "West", "North", "East")
plcvol(ra, rb, id = c(1, 1, 1, 2, 2))
#> plcvol: 2 person-years, records in region order; mean volatility 1.5 (within 1, across 0.5)
#>  id n within across vm
#>   1 3      1      1  2
#>   2 2      1      0  1
plcvol(ra, rb, id = c(1, 1, 1, 2, 2), ord = "given")
#> plcvol: 2 person-years, records in given order; mean volatility 1.5 (within 1, across 0.5)
#>  id n within across vm
#>   1 3      1      1  2
#>   2 2      1      0  1
```

## Alert complexity of each person-year

[`alrtcx()`](https://rootcoder007.github.io/crimlab/reference/alrtcx.md):
Counts the distinct alert configurations among a person-year’s records.
Each record carries flags (in Ontario’s records the mental health,
suicide risk and suicide watch alerts); its configuration is the tuple
of the flags, (1, 0, 0) for a mental health alert alone, coded as the
text `"100"`. Alert complexity is the number of distinct configurations
the person-year’s records span, and `cx`, complexity of two or more, is
the treatment of Ruhela (2026, Appendix H).

``` r

f <- data.frame(mh = c("Yes", "No", "No", "Yes"),
                sr = c("No", "No", "Yes", "Yes"),
                sw = c("No", "No", "No", "Yes"))
alrtcx(f, id = c(1, 1, 2, 3))
#> alrtcx: 3 person-years; 1 (33.3%) span two or more alert configurations
#>  id n ac cx anyal
#>   1 2  2  1     1
#>   2 1  1  0     1
#>   3 1  1  0     1
```

## Person-year table of a restrictive confinement release

[`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md):
Reduces a placement-level release of Ontario’s restrictive confinement
records (table a01 of the open data, one row per person, fiscal year,
pair of regions and alert configuration, with a count of placements) to
one row per person and fiscal year, carrying placement volatility
([`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md)),
alert complexity
([`alrtcx()`](https://rootcoder007.github.io/crimlab/reference/alrtcx.md))
and the covariates of Ruhela (2026): gender, age group, fiscal year and
regional cluster (the regions the person-year’s records name, sorted and
joined by `" + "`).

``` r

d <- data.frame(
  EndFiscalYear = c(2025, 2025, 2025, 2026, 2026),
  UniqueIndividual_ID = c("a", "a", "b", "a", "c"),
  Region_AtTimeOfPlacement = c("East", "West", "East", "North", "East"),
  Region_MostRecentPlacement = c("East", "West", "East", "East", "East"),
  Gender = "Male", Age_Category = "25 to 49",
  MentalHealth_Alert = c("Yes", "No", "No", "No", "Yes"),
  SuicideRisk_Alert = "No", SuicideWatch_Alert = "No",
  Number_Of_Placements = c(3, 1, 2, 4, 1)
)
otpyr(d)
#> otpyr: 4 person-years, fiscal years 2025-2026 (2)
#>  year py  cx  vm
#>  2025  2 0.5 0.5
#>  2026  2 0.0 0.5
```

## Count regression with a random intercept per cluster

[`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md):
Fits the log-linear model of a count with a normal random intercept for
each cluster, log E\[y_ij \| u_j\] = x_ij’ beta + u_j, u_j ~ N(0,
sigma^2), with a Poisson or a negative binomial (variance
`mu + mu^2 / theta`) response, by maximum likelihood with the Laplace
approximation of each cluster’s integral (Breslow and Clayton 1993), as
‘glmmTMB’ and ‘lme4’ fit it. It is the main model of Ruhela (2026):
placement volatility on alert complexity with a random intercept per
regional cluster.

``` r

set.seed(2)
g <- rep(1:12, each = 40)
d <- rbinom(480, 1, 0.3)
y <- rpois(480, exp(-1 + 0.4 * d + rnorm(12, 0, 0.5)[g]))
cntmix(y, data.frame(d = d), g, fam = "poisson")
#> cntmix: poisson counts with a random intercept per cluster
#> 480 counts in 12 clusters; sigma 0.5101; log-likelihood -459.3308
#>         term     est     se      z         p    irr     lo    hi
#>  (Intercept) -0.8465 0.1695 -4.994 5.919e-07 0.4289 0.3077 0.598
#>            d  0.4339 0.1215  3.573 3.535e-04 1.5433 1.2164 1.958
```

## Alert complexity and placement volatility, pooled or by fiscal year

[`cxvol()`](https://rootcoder007.github.io/crimlab/reference/cxvol.md):
The analysis of Ruhela (2026) on any release of Ontario’s restrictive
confinement records: from the person-year table
([`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md))
it gives, for the years pooled or for each fiscal year, the share of
person-years with complex alert profiles (`cx = 1`), mean placement
volatility by group, the crude Poisson rate ratio (the ratio of the two
means, infinite when every count without complex alerts is 0), the
one-to-one nearest-neighbour propensity matching without replacement of
each complex person-year to a control on age group, gender and fiscal
year
([`cbatt()`](https://rootcoder007.github.io/crimlab/reference/cbatt.md),
in ‘MatchIt’’s default order), the count model on the matched frame with
a random intercept per regional cluster
([`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md)),
and the bounds on the average effect on the matched frame with no
assumption and under monotone treatment response and selection
([`mtrbnd()`](https://rootcoder007.github.io/crimlab/reference/mtrbnd.md)).

``` r

set.seed(4)
n <- 600
p <- data.frame(year = sample(2025:2026, n, TRUE),
                gender = sample(c("Male", "Female"), n, TRUE, c(9, 1)),
                age = sample(c("18 to 24", "25 to 49", "50+"), n, TRUE),
                rc = sample(c("A", "B", "A + B", "C"), n, TRUE),
                cx = rbinom(n, 1, 0.2))
p$vm <- rpois(n, exp(-2 + 0.3 * p$cx + (p$rc == "A + B")))
cxvol(p, fam = "poisson")
#> cxvol: alert complexity (cx) and placement volatility (vm), years pooled
#>     sample   n     cx    vm1    vm0 crude matched    irr     lo    hi     p    bound
#>  2025-2026 600 0.1917 0.2174 0.2062 1.054     230 0.8477 0.5012 1.434 0.538 -0.06087
#> irr: rate ratio of the matched count model (poisson); bound: upper end of [0, bound] under monotone response and
#> selection (below 0, the data contradict the two assumptions)
```

## Double machine learning of the effect of alert complexity

[`cxdml()`](https://rootcoder007.github.io/crimlab/reference/cxdml.md):
The interactive regression model of Chernozhukov et al. (2018) for the
effect of complex alert profiles (`cx`) on placement volatility (`vm`)
on the person-year table, by
[`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md):
random forests (or another learner of
[`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md))
for the outcome regressions and the propensity, `nfold`-fold
cross-fitting, propensities kept within `trim` of 0 and 1, covariates
gender, age group, regional cluster and (pooled, with `yterm`) fiscal
year. Both the average effect (ATE) and the effect on the treated (ATTE)
are given, from one set of cross-fitted nuisances. The same refusals as
[`cxvol()`](https://rootcoder007.github.io/crimlab/reference/cxvol.md)
apply.

``` r

set.seed(5)
n <- 400
p <- data.frame(year = sample(2025:2026, n, TRUE),
                gender = sample(c("Male", "Female"), n, TRUE),
                age = sample(c("18 to 24", "25 to 49"), n, TRUE),
                rc = sample(c("A", "B"), n, TRUE), cx = rbinom(n, 1, 0.3))
p$vm <- rpois(n, exp(-1.5 + 0.4 * p$cx))
cxdml(p, lrn = "glm")
#> cxdml: interactive model, effect of cx on vm (additive scale)
#>     sample   n estimand    est      se      lo     hi
#>  2025-2026 400      ATE 0.1571 0.06122 0.03713 0.2771
#>  2025-2026 400     ATTE 0.1336 0.05927 0.01742 0.2498
```

## Placement volatility events by fiscal year

[`cxyoy()`](https://rootcoder007.github.io/crimlab/reference/cxyoy.md):
The year-over-year table of placement volatility (Ruhela 2026, Appendix
E, Table E15) from the person-year table: for each fiscal year the
regional changes (`vm` summed), placements, person-years, the rates of
changes per 1,000 placements and per 1,000 person-years with exact
Poisson intervals, each year’s share of all changes, and the percent
changes from the previous year of the count and of both rates. The rate
intervals are `qgamma(a/2, n)` to `qgamma(1 - a/2, n + 1)` over the
exposure, so a count of 0 has lower limit 0.

``` r

p <- data.frame(year = rep(2024:2026, c(50, 60, 70)), vm = rpois(180, 0.2),
                np = rpois(180, 20) + 1)
cxyoy(p)
#> cxyoy: regional changes per fiscal year
#>  year vm   np py   rnp rnplo rnphi   rpy rpylo rpyhi  share    dvm   drnp   drpy
#>  2024 12 1045 50 11.48 5.934 20.06 240.0 124.0 419.2 0.2667     NA     NA     NA
#>  2025 17 1245 60 13.65 7.954 21.86 283.3 165.1 453.6 0.3778 41.667  18.91  18.06
#>  2026 16 1467 70 10.91 6.234 17.71 228.6 130.6 371.2 0.3556 -5.882 -20.13 -19.33
```

## Simulation of the alert complexity estimators

[`cxsim()`](https://rootcoder007.github.io/crimlab/reference/cxsim.md):
Draws person-year tables from the model of Ruhela (2026) and reports the
bias, spread and interval coverage of the estimators of
[`cxvol()`](https://rootcoder007.github.io/crimlab/reference/cxvol.md).
Each table has `n` person-years in `ncl` regional clusters with normal
intercepts (standard deviation `sdcl`), an age group (three levels) and
a gender, a treatment drawn from a logistic model in which the
covariates and the cluster raise the chance of a complex profile, and a
negative binomial count (dispersion `theta`) whose log mean is
`b0 + beta * cx` plus the covariates and the cluster. With `conf` the
person-year also has an unobserved trait that raises both the chance of
treatment and the count (by `conf` on the logit and on the log mean):
ignorability given the observed covariates then fails.

``` r

cxsim(reps = 1, n = 100, ncl = 5, b0 = 0, seed = 1)
#> cxsim: 1 tables of 100 person-years, log rate ratio 0.2624
#>      est   bias sd   rmse cover ok
#>    crude 0.4090 NA 0.4090    NA  1
#>    mixed 0.1381 NA 0.1381     1  1
#>  matched 0.0113 NA 0.0113     1  1
```
