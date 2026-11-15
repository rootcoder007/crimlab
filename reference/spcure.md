# Split-population (cure) model of recidivism

Asks two questions that an ordinary survival model runs together: what
share of people released will ever be reconvicted, and how soon those
who are reconvicted fail. A Kaplan-Meier curve of time to reconviction
levels off well above 0: many released people never fail, however long
they are followed. Schmidt and Witte (1989) model this with a split
population: each person would eventually fail with chance
`delta = 1 / (1 + exp(-x'gamma))` (the `"fail"` part), and the time to
failure of those who would fail follows a distribution with scale
`exp(x'beta)` (the `"time"` part): a Weibull with a common shape, an
exponential, or a lognormal (log time normal with mean `x'beta` and a
common standard deviation `sigma`). A person followed to time `t`
without failing contributes `1 - delta + delta S(t)`, one who failed at
`t` contributes `delta f(t)`. The parameters are estimated by maximum
likelihood (quasi-Newton steps with the analytic gradient) with standard
errors from the inverse of the Hessian.

## Usage

``` r
spcure(time, event, x = NULL, dist = "weibull", level = 0.95, maxit = 500)
```

## Arguments

- time:

  The time each person was followed to failure or to the end of
  follow-up: positive numbers.

- event:

  Whether each time is a failure (1 or `TRUE`) or the end of follow-up
  (0 or `FALSE`); failures and censored times both.

- x:

  `NULL` or the covariates of both parts: a numeric matrix, a data frame
  or list of numeric columns, or one numeric vector, one row per person.

- dist:

  The time to failure: `"weibull"`, `"exp"` or `"lognormal"`.

- level:

  The confidence level of the Wald intervals: one number strictly
  between 0 and 1.

- maxit:

  The quasi-Newton steps: one whole number from 1.

## Value

An object of class `spcure`, a list: `coef`, one row per parameter
(`part`: `"fail"`, `"time"` or `"shape"`; `term`, `est`, `se`, `lo`,
`hi`, `z`, `p`); `fit`, a named vector (`loglik`, `aic`, `iter`, `conv`
(1 when converged), `pfail`); `dist`; and `n`.

## Details

The two parts are told apart by the long right tail of the follow-up: if
everyone is followed only briefly, a low share failing and slow failure
look alike and the estimates are unstable (large standard errors).
`pfail` is the mean of `delta` over the people, the estimated share who
would ever fail.

## References

Schmidt, P. and Witte, A. D. (1989). Predicting criminal recidivism
using 'split population' survival time models. *Journal of
Econometrics*, 40(1), 141-159.
[doi:10.1016/0304-4076(89)90034-1](https://doi.org/10.1016/0304-4076%2889%2990034-1)

Maltz, M. D. (1984). *Recidivism*. Academic Press, Orlando.

## See also

[`kmrec()`](https://rootcoder007.github.io/crimlab/reference/kmrec.md)
for the Kaplan-Meier curve.

## Examples

``` r
set.seed(10)
n <- 600
prior <- rbinom(n, 1, 0.4)
fail <- rbinom(n, 1, plogis(-0.5 + 1.2 * prior))
tt <- ifelse(fail == 1, rweibull(n, 1.3, 300), Inf)
cens <- runif(n, 365, 1500)
r <- spcure(pmin(tt, cens), as.numeric(tt <= cens), data.frame(prior))
r
#> spcure: 600 people, weibull times; share who would ever fail 0.5268; log-likelihood -2343.2
#>   part        term        est         se         lo         hi          z
#>   fail (Intercept) -0.3366152 0.11268803 -0.5574797 -0.1157507 -2.9871424
#>   fail       prior  1.1306239 0.18545011  0.7671484  1.4940995  6.0966475
#>   time (Intercept)  5.7522378 0.06907445  5.6168543  5.8876212 83.2759080
#>   time       prior -0.0414838 0.09364132 -0.2250174  0.1420498 -0.4430074
#>  shape  log(shape)  0.3494171 0.04854205  0.2542765  0.4445578  7.1982364
#>             p
#>  2.815985e-03
#>  1.083160e-09
#>  0.000000e+00
#>  6.577604e-01
#>  6.099631e-13
```
