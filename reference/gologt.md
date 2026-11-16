# The generalised ordered logit

An ordered outcome (a risk level, a charge's seriousness, a sentence
band) regressed on covariates without forcing each covariate to act the
same way at every cut of the scale. For an outcome with levels
`1, ..., M`, each threshold `j < M` has its own intercept and slopes,
`P(y > j) = F(a_j + x b_j)` with `F` the logistic distribution function
(Peterson and Harrell 1990; Williams 2006, the parameterisation of
'gologit2'); a covariate named in `par` keeps one slope for all the
thresholds, as in the proportional-odds model, so that `par = x` is that
model. O'Connell and Laniyonu (2025) find racial disparity in a risk
tool concentrated at the low-to-medium cut and absent at the
medium-to-high one, which a proportional-odds fit averages away.

## Usage

``` r
gologt(data, y, x, par = NULL, lev = NULL, level = 0.95, maxit = 500)

# S3 method for class 'gologt'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit.

- y:

  The name of the ordered outcome column: an ordered factor or a factor
  (its levels in order), numbers (in increasing order) or text (with
  `lev`).

- x:

  The names of the covariate columns: numeric or logical, at least one.

- par:

  `NULL`, or the names of covariates in `x` held parallel (one slope for
  every threshold).

- lev:

  `NULL`, or the levels of the outcome from lowest to highest, each
  present in the data; needed for text.

- level:

  The confidence level of the Wald intervals: one number strictly
  between 0 and 1.

- maxit:

  The quasi-Newton steps: one whole number from 1.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `gologt`, a list: `coef`, one row per parameter
(`term`, `thr` the threshold, as `"low|med"`, or `"all"` for a parallel
slope; `est`, `se`, `z`, `p`, `lo`, `hi`); `wald`, one row per covariate
(`term`, `stat`, `df`, `p`; `NA` for a parallel one); `lr`, the test
against proportional odds (`stat`, `df`, `p`; `NA` when every slope is
parallel); `fit` (`loglik`, `aic`, `iter`, `conv`, 1 when converged);
`vcov`; `lev`; and `n`.

## Details

The parameters are estimated by maximum likelihood (quasi-Newton steps
with the analytic gradient, in the compiled core), with standard errors
from the inverse of the Hessian. A free covariate's slopes may cross, so
that a category's fitted probability is negative for some `x`; the fit
refuses such points within the data. `wald` tests, for each free
covariate, that its slopes are equal at every threshold (the test
'gologit2' uses to choose the parallel ones), and `lr` is the
likelihood-ratio test of the fit against the proportional-odds model.

## References

Peterson, B. and Harrell, F. E. (1990). Partial proportional odds models
for ordinal response variables. *Journal of the Royal Statistical
Society Series C*, 39(2), 205-217.
[doi:10.2307/2347760](https://doi.org/10.2307/2347760)

Williams, R. (2006). Generalized ordered logit/partial proportional odds
models for ordinal dependent variables. *The Stata Journal*, 6(1),
58-82.
[doi:10.1177/1536867X0600600104](https://doi.org/10.1177/1536867X0600600104)

O'Connell, M. and Laniyonu, A. (2025). Threshold-specific
cumulative-logit models for actuarial-risk audit. *Race and Justice*,
15(3), 428-453.

## See also

[`gocoef()`](https://rootcoder007.github.io/crimlab/reference/gocoef.md)
for the slopes by threshold.

## Examples

``` r
set.seed(7)
n <- 400
d <- data.frame(race = rbinom(n, 1, 0.4), age = rnorm(n))
e1 <- 0.5 + 0.8 * d$race - 0.5 * d$age
e2 <- -1 + 0.1 * d$race - 0.5 * d$age
u <- runif(n)
d$risk <- factor(1 + (u < plogis(e1)) + (u < plogis(e2)), 1:3,
                 c("low", "med", "high"))
r <- gologt(d, "risk", c("race", "age"), par = "age")
r
#> gologt: generalised ordered logit, 400 rows, 3 levels; log-likelihood -414.564
#>         term      thr        est         se          z            p         lo
#>  (Intercept)  low|med  0.5549532 0.14198273  3.9085966 9.283381e-05  0.2766722
#>  (Intercept) med|high -0.8280251 0.14776492 -5.6036652 2.098659e-08 -1.1176390
#>         race  low|med  0.5972905 0.22943154  2.6033496 9.231777e-03  0.1476129
#>         race med|high -0.1312269 0.22560745 -0.5816604 5.607954e-01 -0.5734094
#>          age      all -0.5709891 0.09757872 -5.8515741 4.869423e-09 -0.7622399
#>          hi
#>   0.8332342
#>  -0.5384112
#>   1.0469681
#>   0.3109556
#>  -0.3797383
#> against proportional odds: LR 8.875 on 1 df, p 0.00289
```
