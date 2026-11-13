# Rates and proportions by group: ratios and differences

The rate (events over exposure) or the proportion (events over trials)
of each group, and its ratio to and difference from a reference group,
with approximate and exact intervals. Rows that share a group are
pooled: give one row per group, or one per case (0/1 outcomes with
`n = 1`, the default for proportions), or anything between.

## Usage

``` r
grprr(
  x,
  n = NULL,
  grp = NULL,
  ref = NULL,
  type = "rate",
  per = 1,
  level = 0.95
)
```

## Arguments

- x:

  The events: numbers from 0; or a `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md).

- n:

  The exposure (person-time, population) of each row for rates, or the
  trials for proportions: positive numbers; `NULL` for one trial per row
  (proportions only).

- grp:

  `NULL` (each row is a group, named by the names of `x`) or the group
  of each row: labels.

- ref:

  The reference group: a label or a position; the first group by
  default.

- type:

  `"rate"` (events over exposure) or `"prop"` (a proportion).

- per:

  The rates and differences are given per `per` units of `n` (1e5 for
  rates per 100,000, 100 for percentages): one positive number.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `grprr`, a list: `rates`, a data frame with one row
per group (`group`, `x`, `n`, `rate`, `lo`, `hi` the approximate
interval, `elo`, `ehi` the exact one, `rep`); `cmp`, one row per group
but the reference (`group`, `ratio`, `lo`, `hi`, `elo`, `ehi`, `diff`,
`dlo`, `dhi`); `ref`, `type`, `per` and `level`.

## Details

For a rate `x / n` the approximate interval is
`rate * exp(-+z / sqrt(x))` and the exact one the Poisson interval
through the gamma quantiles, as
[`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html);
the ratio's approximate interval is log-normal with variance
`1 / x + 1 / x0`, its exact interval the one conditional on the total
count (the Clopper-Pearson interval of `x / (x + x0)`, as the two-sample
[`stats::poisson.test()`](https://rdrr.io/r/stats/poisson.test.html)),
and the difference's is Wald. For a proportion the intervals are
Wilson's score interval and the exact Clopper-Pearson one
([`stats::binom.test()`](https://rdrr.io/r/stats/binom.test.html)); the
ratio's is log-normal (Katz et al.), with no exact interval (`NA`), and
the difference's is Newcombe's hybrid score interval. With a zero count
the log-normal intervals are `NA`: no continuity correction is made.

`rep` is a group's share of the events over its share of the exposure
(or the trials): over 1 when the group is over-represented among the
events, as in a comparison of who is placed in segregation with who is
in custody.

## References

Katz, D., Baptista, J., Azen, S. P. and Pike, M. C. (1978). Obtaining
confidence intervals for the risk ratio in cohort studies. *Biometrics*,
34(3), 469-474. [doi:10.2307/2530610](https://doi.org/10.2307/2530610)

Newcombe, R. G. (1998). Interval estimation for the difference between
independent proportions: comparison of eleven methods. *Statistics in
Medicine*, 17(8), 873-890.
[doi:10.1002/(SICI)1097-0258(19980430)17:8\<873::AID-SIM779\>3.0.CO;2-I](https://doi.org/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C873%3A%3AAID-SIM779%3E3.0.CO%3B2-I)

## Examples

``` r
# use-of-force incidents per 1,000 contacts by perceived race
grprr(c(white = 120, black = 95, other = 30),
      c(white = 60000, black = 21000, other = 15000), per = 1000)
#> grprr: rates per 1000, 95% intervals
#>  group   x     n    rate       lo       hi      elo      ehi       rep
#>  white 120 60000 2.00000 1.672347 2.391849 1.658199 2.391510 0.7836735
#>  black  95 21000 4.52381 3.699754 5.531409 3.660032 5.530123 1.7725948
#>  other  30 15000 2.00000 1.398372 2.860470 1.349392 2.855124 0.7836735
#> against white:
#>  group    ratio        lo       hi       elo      ehi    diff        dlo
#>  black 2.261905 1.7281381 2.960535 1.7096642 2.985434 2.52381  1.5462750
#>  other 1.000000 0.6702691 1.491938 0.6468078 1.502337 0.00000 -0.8001519
#>        dhi
#>  3.5013441
#>  0.8001519
# proportions from one row per case
grprr(c(1, 0, 0, 1, 1, 0, 1, 1), grp = c(1, 1, 1, 1, 2, 2, 2, 2),
      type = "prop")
#> grprr: proportions, 95% intervals
#>  group x n rate        lo        hi        elo       ehi rep
#>      1 2 4 0.50 0.1500390 0.8499610 0.06758599 0.9324140 0.8
#>      2 3 4 0.75 0.3006418 0.9544127 0.19412045 0.9936905 1.2
#> against 1:
#>  group ratio        lo       hi elo ehi diff        dlo       dhi
#>      2   1.5 0.4837821 4.650854  NA  NA 0.25 -0.3195573 0.6552867
```
