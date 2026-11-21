# Confidence interval of a survey proportion

An interval for a proportion estimated from a complex survey, from the
estimate and its design-based standard error (those of
[`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md),
say). The usual interval, the estimate plus and minus a multiple of its
standard error (Wald), is too narrow for a small sample or a proportion
near 0 or 1, and can leave the range 0 to 1; statistical agencies ask
for one of the three others here instead.

## Usage

``` r
svypci(
  est,
  se,
  n,
  df = n - 1,
  method = c("beta", "wilson", "logit", "wald"),
  level = 0.95,
  trunc = TRUE
)
```

## Arguments

- est:

  The estimated proportion, from 0 to 1.

- se:

  Its standard error.

- n:

  The number of respondents of the sample (unweighted).

- df:

  The degrees of freedom of the design: first-stage units less strata,
  as
  [`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
  reports; `n - 1` by default (a simple random sample).

- method:

  `"beta"` (the default), `"wilson"`, `"logit"` or `"wald"`.

- level:

  The confidence level of the intervals.

- trunc:

  `TRUE` (the default) keeps the effective sample size at or below `n`.

## Value

A data frame with a row per estimate: `est`, `lo`, `hi` and `neff`, the
effective sample size `ne*` (for every method: it is what the design is
worth). `NA` where an input is missing, and for the logit interval at 0
or 1.

## Details

The first two put the design into an effective sample size, the size of
a simple random sample with the same variance, `ne = p (1 - p) / se^2`,
corrected for the degrees of freedom of the design (Korn and Graubard
1998): `ne* = ne (t(n - 1) / t(df))^2`, with
[`t()`](https://rdrr.io/r/base/t.html) the quantile of Student's t at
the level. With `trunc`, the default, `ne*` is never above `n`, as Korn
and Graubard and SAS have it; `survey::svyciprop()` does not truncate,
and `trunc = FALSE` gives its interval. An estimate of exactly 0 or 1,
or a standard error of 0, has no effective size of its own and takes
`n`.

`"beta"`: the Clopper-Pearson interval at `ne*` with `p ne*` events, by
beta quantiles (Korn and Graubard 1998): the "modified Clopper-Pearson"
interval. `"wilson"`: the Wilson score interval at `ne*` with the normal
quantile: the "modified Wilson" interval (Dean and Pagano 2015).
`"logit"`: the Wald interval of `log(p / (1 - p))`, with standard error
`se / (p (1 - p))` and Student's t on `df`, taken back to a proportion;
it does not exist at 0 or 1. `"wald"`: `p` plus and minus the t quantile
times `se`, the interval
[`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
gives.

Health Canada's guides to the Canadian surveys of alcohol and drug use
recommend, in this order, the modified Wilson, the modified
Clopper-Pearson and the logit interval, and advise against the Wald
interval (Health Canada 2022, section 7.3).

Errors: `est` outside 0 to 1; `se` negative; `n` not a whole number of
at least 2; `df` below 1; `level` outside (0, 1); `method` not one of
the four; `trunc` not `TRUE` or `FALSE`; arguments that are not numeric
or of different lengths (a length of 1 is repeated).

## References

Korn, E. L. and Graubard, B. I. (1998). Confidence intervals for
proportions with small expected number of positive counts estimated from
survey data. *Survey Methodology*, 24(2), 193-201.

Dean, N. and Pagano, M. (2015). Evaluating confidence interval methods
for binomial proportions in clustered surveys. *Journal of Survey
Statistics and Methodology*, 3(4), 484-503.
[doi:10.1093/jssam/smv024](https://doi.org/10.1093/jssam/smv024)

Health Canada (2022). *Canadian Postsecondary Education Alcohol and Drug
Use Survey, 2021-22: Public Use Microdata File user guide*. Section 7.3,
release guidelines.

## See also

[`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
for the estimate and its standard error;
[`svyqual()`](https://rootcoder007.github.io/crimlab/reference/svyqual.md)
for the quality level of an estimate.

## Examples

``` r
# 12% of 400 respondents, standard error 2.1 points, 38 design df
svypci(0.12, 0.021, 400, df = 38)
#>    est         lo      hi    neff
#> 1 0.12 0.08065552 0.16967 225.824
for (m in c("beta", "wilson", "logit", "wald")) {
  print(svypci(0.12, 0.021, 400, df = 38, method = m))
}
#>    est         lo      hi    neff
#> 1 0.12 0.08065552 0.16967 225.824
#>    est         lo        hi    neff
#> 1 0.12 0.08385073 0.1688613 225.824
#>    est         lo        hi    neff
#> 1 0.12 0.08355411 0.1694048 225.824
#>    est         lo        hi    neff
#> 1 0.12 0.07748772 0.1625123 225.824
# near 0 the Wald interval leaves the range; the others do not
svypci(0.01, 0.008, 400, df = 38, method = "wald")
#>    est           lo         hi     neff
#> 1 0.01 -0.006195153 0.02619515 145.8814
svypci(0.01, 0.008, 400, df = 38, method = "wilson")
#>    est          lo         hi     neff
#> 1 0.01 0.002272718 0.04287127 145.8814
try(svypci(1.2, 0.1, 400))
#> Error : svypci: `est` must be proportions, from 0 to 1.
```
