# Exposure-response estimate with a bootstrap sensitivity interval

Fits the partially linear exposure-response model (the effect of an
exposure on an outcome, the confounders partialled out by cross-fitted
regressions), or any estimator with the same interface, then re-fits it
on `nboot` non-parametric resamples of the rows for a percentile
interval (Efron and Tibshirani 1993, chapter 13). The bootstrap interval
does not lean on the estimator's own standard error, so it is a check on
it: when the two disagree, the analytic one should not be trusted alone.

## Usage

``` r
exprspplr(data, y, expo, confs, nboot = 100L, rstate = 42L, estr = NULL)

# S3 method for class 'exprspplr'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit (an area and period, a person).

- y, expo, confs:

  Column names: the outcome and the exposure (single strings) and the
  confounders (one or more strings).

- nboot:

  Number of resamples, a non-negative whole number (default 100); at
  least 10 must succeed.

- rstate:

  Seed for the resampling (default 42), a whole number from -2147483647
  to \\2^{53} - 1\\.

- estr:

  `NULL` (the default) for the built-in cross-fitted partially linear
  regression, or a function (or the name of one)
  `(data, outcome, treatment, covariates)` returning a list with `ate`
  and `se` (and optionally `method`).

- x:

  A result of `exprspplr()`.

- ...:

  Not used.

## Value

A list of class `exprspplr` with `ate` (the full-sample estimate),
`seana` (the estimator's standard error), `seboot` (the standard
deviation of the replicates), `cilobs` and `ciupbs` (the 95\\ number of
successful resamples) and `method` (the estimator's method with
`" + percentile bootstrap"`).

## Details

In environmental criminology the exposure is a pollutant or heat and the
outcome a count of offences or victimisations by place and period, with
the neighbourhood's population and deprivation as confounders: the
estimate is the change in the outcome per unit of exposure.

The built-in estimator (`estr = NULL`) is the partially linear
regression \\y = \theta d + g(x) + e\\ fitted by double machine learning
with cross-fitting (the DML2 estimator of Chernozhukov et al. 2018): the
rows with every column present are split at random into five folds, the
outcome and the exposure are each regressed on the confounders by
ordinary least squares on four folds and predicted on the fifth, and
\\\theta\\ is the regression of the outcome's residuals on the
exposure's, with the heteroskedasticity-robust standard error. The folds
are drawn from R's generator seeded with 1 for every fit, the stream
then restored, so the estimate does not depend on the session. It
refuses a column that is not numeric, fewer than ten complete rows, and
an exposure with no variation left after the confounders are partialled
out.

Each resample draws `nrow(data)` rows with replacement and calls
`estr(resample, y, expo, confs)`; a resample on which the estimator
fails is skipped, and the run stops if fewer than 10 succeed. The
interval is the 2.5\\ `seboot` their standard deviation. Resampling uses
R's generator seeded with `rstate` and restores R's own random stream
afterwards. The cost is `nboot + 1` estimator calls; the built-in
estimator runs in the compiled core.

Refusals: column names that are not strings, an `nboot` that is not a
single non-negative whole number, an estimator that is not a function,
and the estimator's own errors on the full sample.

## References

Efron, B. and Tibshirani, R. J. (1993). An Introduction to the
Bootstrap. Chapman and Hall.

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
Newey, W. and Robins, J. (2018). Double/debiased machine learning for
treatment and structural parameters. The Econometrics Journal 21(1),
C1-C68. [doi:10.1111/ectj.12097](https://doi.org/10.1111/ectj.12097)

## See also

[`dmlx`](https://rootcoder007.github.io/crimlab/reference/dmlx.md) for
double machine learning with other learners and models.

## Examples

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
