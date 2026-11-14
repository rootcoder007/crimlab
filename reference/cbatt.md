# Design-based causal effects of a 0/1 treatment

The battery of estimators of the effect of a 0/1 treatment `d` on an
outcome `y` given covariates `x` that do not learn the nuisances by
cross-fitting: Hajek-normalised inverse probability weighting,
g-computation, subclassification on the propensity, nearest-neighbour
matching on the propensity with a caliper, and entropy balancing,
covariate-balancing weights for the effect on the treated.
[`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md)
gives the cross-fitted ones (augmented inverse probability weighting,
targeted maximum likelihood, the interactive model).

## Usage

``` r
cbatt(
  y,
  d,
  x,
  est = "ipw",
  tgt = "ate",
  ps = NULL,
  k = 1,
  calip = NULL,
  nsub = 5,
  reps = 499,
  seed = NULL,
  level = 0.95
)
```

## Arguments

- y:

  The outcome: finite numbers, one per unit.

- d:

  The treatment: 0 or 1 (or `FALSE` or `TRUE`), two units in each arm at
  least.

- x:

  The covariates: a numeric matrix, a data frame or list of numeric
  columns, or one numeric vector, one row per unit.

- est:

  `"ipw"`, `"gcomp"`, `"sub"`, `"psm"` or `"ebal"`.

- tgt:

  `"ate"`, `"att"` or `"atc"`.

- ps:

  `NULL` or the propensity of each unit, strictly between 0 and 1.

- k:

  The controls matched to each unit: a whole number from 1.

- calip:

  `NULL` or the caliper in standard deviations of the propensity.

- nsub:

  The subclasses: a whole number from 2.

- reps:

  The bootstrap samples: a whole number from 0 (none).

- seed:

  The seed of the bootstrap: `NULL` to draw it from R's generator, or a
  whole number from 0.

- level:

  The confidence level of the normal interval.

## Value

An object of class `cbatt`, a list: `res`, one row (`est`, `se`, `lo`,
`hi`, `med`, `n1`, `n0`); `ps`, the propensities (`NA` for `"gcomp"` and
`"ebal"`); `w`, the weights (`"ipw"`, `"ebal"`; else `NULL`); `match`
(`"psm"`: the rows of each `focal` unit and its `match`; else `NULL`);
`est`; `tgt`; and `level`.

## Details

The propensity is the logistic regression of `d` on an intercept and `x`
(the learner of
[`lgtreg()`](https://rootcoder007.github.io/crimlab/reference/lgtreg.md))
unless `ps` gives it. `tgt` picks the average effect (`"ate"`), the
effect on the treated (`"att"`) or on the controls (`"atc"`).

- `"ipw"`: the weighted mean of the treated less that of the controls,
  each normalised to sum to one (Hajek), with weights `1 / e` and
  `1 / (1 - e)` (ATE), `1` and `e / (1 - e)` (ATT), `(1 - e) / e` and
  `1` (ATC).

- `"gcomp"`: least squares of `y` on `x` in each arm; the mean
  difference of the two predictions over the target units.

- `"sub"`: `nsub` subclasses cut at the quantiles of the propensity
  among the target units; the subclasses' differences of means weighted
  by their target units.

- `"psm"`: greedy `k`:1 matching without replacement on the propensity,
  the treated (ATT) or the controls (ATC) in decreasing propensity
  ('MatchIt”s default order), within `calip` standard deviations of the
  propensity. The estimate is the mean over the matched units of the
  treated outcome less its controls' mean; `med` is their median, the
  median causal effect of the matched pairs; the standard error is that
  of a mean of independent matched sets.

- `"ebal"`: entropy balancing (Hainmueller 2012): control weights
  proportional to `exp(l'x)`, the nearest to uniform in entropy whose
  weighted means equal the treated means exactly; the ATT is the treated
  mean less the weighted control mean.

Standard errors (but `"psm"`) are the standard deviations over `reps`
bootstrap samples of the units, with the propensity fitted again in
each; intervals are normal.

A `clpre` list from
[`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
may be given as `y`: its elements fill the arguments by name (the
`otis-*-cbatt` presets).

## References

Hainmueller, J. (2012). Entropy balancing for causal effects: a
multivariate reweighting method to produce balanced samples in
observational studies. *Political Analysis*, 20(1), 25-46.
[doi:10.1093/pan/mpr025](https://doi.org/10.1093/pan/mpr025)

Imbens, G. W. and Rubin, D. B. (2015). *Causal Inference for Statistics,
Social, and Biomedical Sciences*. Cambridge University Press.
[doi:10.1017/CBO9781139025751](https://doi.org/10.1017/CBO9781139025751)

Austin, P. C. (2011). An introduction to propensity score methods for
reducing the effects of confounding in observational studies.
*Multivariate Behavioral Research*, 46(3), 399-424.
[doi:10.1080/00273171.2011.568786](https://doi.org/10.1080/00273171.2011.568786)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(400), 200)
d <- rbinom(200, 1, plogis(x[, 1]))
y <- 2 * d + x[, 1] - x[, 2] + rnorm(200)
cbatt(y, d, x, reps = 99, seed = 1)$res
#>        est        se       lo      hi med n1  n0
#> 1 2.077253 0.1680322 1.747917 2.40659  NA 99 101
cbatt(y, d, x, est = "ebal", tgt = "att", reps = 99, seed = 1)$res
#>        est        se       lo       hi med n1  n0
#> 1 2.067763 0.1830696 1.708953 2.426573  NA 99 101
cbatt(y, d, x, est = "psm", tgt = "att", calip = 0.2)$res
#>        est        se       lo       hi      med n1  n0
#> 1 2.011912 0.2508233 1.520307 2.503517 2.228522 99 101
```
