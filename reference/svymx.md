# Design-based estimates from a complex survey

Means, totals, proportions, ratios and variances from a stratified,
clustered, weighted sample (a victimisation survey, say), with standard
errors by linearisation that respect the design, by domain and with the
differences between domains; and a description of the design.

## Usage

``` r
svymx(
  y,
  w = NULL,
  strata = NULL,
  psu = NULL,
  fpc = NULL,
  stat = "mean",
  den = NULL,
  by = NULL,
  ref = NULL,
  level = 0.95
)
```

## Arguments

- y:

  The analysis variable: finite numbers, one per respondent (0 or 1 for
  `"prop"`).

- w:

  `NULL` or the sampling weights: positive numbers.

- strata:

  `NULL` or the stratum of each respondent: labels.

- psu:

  `NULL` (each respondent its own) or the first-stage unit of each
  respondent: labels.

- fpc:

  `NULL` or the number of first-stage units in the population of each
  respondent's stratum (constant within strata).

- stat:

  `"mean"`, `"total"`, `"prop"`, `"ratio"` or `"var"`.

- den:

  For `"ratio"`, the denominator variable.

- by:

  `NULL` or the domain of each respondent: labels.

- ref:

  The reference domain of the comparisons: a label or a position; the
  first by default.

- level:

  The confidence level of the intervals.

## Value

An object of class `svymx`, a list: `est`, one row per domain (`dom`,
`est`, `se`, `lo`, `hi`, `cv` the coefficient of variation in percent,
`n` the respondents, `df`); `cmp`, each domain against the reference
(`dom`, `diff`, `se`, `stat`, `p`; `NULL` without `by`); `des`, one row
(`n`, `strata`, `psu`, `sumw`, `df`, `effn`, `deffw`); `stat`; and
`level`.

## Details

Each estimate is a smooth function of weighted totals; its variance is
that of the total of its influence values `u_i` (for a mean
`w_i (y_i - m) / sum w`, a ratio `w_i (y_i - R x_i) / sum w x`, a total
`w_i y_i`) over the first-stage units within strata: \`sum_h (1 - f_h)
n_h / (n_h -

1.  sum_j (U_hj - Ubar_h)^2`, with `f_h = n_h /
    N_h`from`fpc`or 0 (first-stage units drawn with replacement). These are the estimators of the 'survey' package (Lumley 2004) with its default treatment of domains: a domain's influence values are 0 outside it and every first-stage unit stays in the variance, so the domain's own size is random, as it is.`"var"`is the population variance, the mean of`(y -
    ybar)^2 n / (n -
    1)` (`n`the units of the domain), as`survey::svyvar()\` gives it.
    Intervals use Student's t on the design's degrees of freedom,
    first-stage units less strata. A first-stage unit label is read
    within its stratum.

Without `w`, `strata` and `psu` the sample is simple random: the
variance is `s^2 / n` with `s^2` on `n - 1` (Bessel's correction, for a
sample from a large population); with `fpc` the population size,
`(1 - n / N) s^2 / n` (a sample without replacement from a finite
population), 0 when `n = N` (a census describes its population exactly).

`des` describes the design: units, strata, first-stage units, the sum of
the weights (the estimated population), the degrees of freedom, and
Kish's effective sample size `(sum w)^2 / sum w^2` and the design effect
of the weights `n sum w^2 / (sum w)^2`.

## References

Lumley, T. (2004). Analysis of complex survey samples. *Journal of
Statistical Software*, 9(8), 1-19.
[doi:10.18637/jss.v009.i08](https://doi.org/10.18637/jss.v009.i08)

Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
Springer.
[doi:10.1007/978-0-387-35099-8](https://doi.org/10.1007/978-0-387-35099-8)

Kish, L. (1965). *Survey Sampling*. Wiley, New York. ISBN
978-0-471-48900-9.

## Examples

``` r
set.seed(1)
d <- data.frame(region = rep(c("north", "south"), each = 60),
                block = rep(1:12, each = 10),
                victim = rbinom(120, 1, rep(c(0.2, 0.35), each = 60)),
                w = runif(120, 50, 150))
svymx(d$victim, d$w, d$region, d$block, stat = "prop")
#> svymx: prop, 120 respondents, 2 strata, 12 first-stage units
#>  dom       est         se        lo        hi      cv   n df
#>  all 0.2518814 0.02405017 0.1982943 0.3054685 9.54821 120 10
svymx(d$victim, d$w, d$region, d$block, stat = "prop", by = d$region)$cmp
#>     dom      diff         se     stat            p
#> 1 south 0.2234697 0.04419731 5.056184 0.0004947409
```
