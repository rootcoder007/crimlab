# Simulation of the alert complexity estimators

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

## Usage

``` r
cxsim(
  reps = 200,
  n = 2000,
  beta = log(1.3),
  ncl = 24,
  sdcl = 1,
  theta = 2,
  b0 = -2.5,
  conf = 0,
  seed = 1
)
```

## Arguments

- reps:

  The replicates: a whole number from 1 to 10000000.

- n:

  The person-years of each table: a whole number from 50 to 10000000.

- beta:

  The log rate ratio of `cx`, from -20 to 20.

- ncl, sdcl:

  The clusters and the standard deviation of their intercepts.

- theta:

  The negative binomial dispersion.

- b0:

  The intercept of the log mean, from -20 to 20.

- conf:

  0, or the effect of the unobserved trait, from -20 to 20.

- seed:

  The seed of the whole run, set once before the first table is drawn:
  one whole number. The session's own random stream is left as it was.

## Value

An object of class `cxsim`, a list: `res`, one row per estimator (`est`,
`bias`, `sd`, `rmse`, `cover`, `ok` the replicates fitted); `draws`,
every replicate's estimates; and the settings.

## Details

The estimators: `crude`, the Poisson rate ratio of `cx` alone; `mixed`,
[`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md)
on the full table with the covariates (log rate ratio, Wald interval);
`matched`,
[`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md)
on the frame matched as in
[`cxvol()`](https://rootcoder007.github.io/crimlab/reference/cxvol.md).
The target is `beta`, the log rate ratio given the covariates and the
cluster, which the model holds constant
([`proof()`](https://rootcoder007.github.io/crimlab/reference/proof.md)
names the log-link identity). Results are `bias` (mean estimate less
`beta`), `sd`, `rmse`, and `cover`, the share of 95% intervals holding
`beta` (`NA` for `crude`, which reports no interval). Each replicate is
a fresh table from `seed`; a fit that fails is left out and counted.

## Examples

``` r
cxsim(reps = 1, n = 100, ncl = 5, b0 = 0, seed = 1)
#> cxsim: 1 tables of 100 person-years, log rate ratio 0.2624
#>      est   bias sd   rmse cover ok
#>    crude 0.4090 NA 0.4090    NA  1
#>    mixed 0.1381 NA 0.1381     1  1
#>  matched 0.0113 NA 0.0113     1  1
```
