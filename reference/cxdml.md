# Double machine learning of the effect of alert complexity

The interactive regression model of Chernozhukov et al. (2018) for the
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

## Usage

``` r
cxdml(
  x,
  by = c("pool", "year"),
  yterm = TRUE,
  lrn = "forest",
  ntree = 500,
  nfold = 3,
  trim = 0.01,
  seed = 1,
  ...
)
```

## Arguments

- x:

  A person-year table from
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md),
  or what
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  takes (a data frame of the release or a table key).

- by:

  `"pool"` for the years together, `"year"` for each year.

- yterm:

  `TRUE` to adjust for the fiscal year in the pooled propensity and
  count models.

- lrn:

  The learner of
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md),
  for both nuisances.

- ntree:

  The trees of a forest.

- nfold:

  The folds of the cross-fitting: a whole number from 2.

- trim:

  Propensities are kept within `trim` of 0 and 1.

- seed:

  The seed of the folds and of the learners.

- ...:

  Passed to
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  when `x` is not a person-year table.

## Value

An object of class `cxdml`: a data frame with one row per sample and
estimand (`sample`, `n`, `estimand`, `est`, `se`, `lo`, `hi`).

## References

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
Newey, W. and Robins, J. (2018). Double/debiased machine learning for
treatment and structural parameters. *The Econometrics Journal*, 21(1),
C1-C68. [doi:10.1111/ectj.12097](https://doi.org/10.1111/ectj.12097)

## Examples

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
