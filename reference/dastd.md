# Directly standardised rates

The rate each group would have with the stratum make-up (age bands, say)
of a standard population: the sum over strata of the standard's weights
times the group's stratum rates. The interval is Fay and Feuer's gamma
interval; the ratio of two standardised rates has a log-normal interval
and their difference a normal one.

## Usage

``` r
dastd(x, n, strat, grp = NULL, std, ref = NULL, per = 1, level = 0.95)
```

## Arguments

- x:

  The events: numbers from 0; or a `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md).

- n:

  The exposure of each row: positive numbers.

- strat:

  The stratum of each row: labels.

- grp:

  `NULL` (one group) or the group of each row: labels. Each group has
  one row per stratum.

- std:

  The standard population (or weights) of each stratum: numbers from 0
  with a positive sum, named by stratum or in the order of the strata (a
  factor's levels, or the sorted labels).

- ref:

  The reference group: a label or a position; the first group by
  default.

- per:

  The rates and differences are given per `per` units of `n` (1e5 for
  rates per 100,000, 100 for percentages): one positive number.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `dastd`, a list: `rates` (`group`, `x` the events,
`rate`, `lo`, `hi`, `var` the variance of the rate), `cmp` (one row per
group but the reference: `group`, `ratio`, `lo`, `hi`, `diff`, `dlo`,
`dhi`), `ref`, `per` and `level`.

## References

Fay, M. P. and Feuer, E. J. (1997). Confidence intervals for directly
standardized rates: a method based on the gamma distribution.
*Statistics in Medicine*, 16(7), 791-801.
[doi:10.1002/(SICI)1097-0258(19970415)16:7\<791::AID-SIM500\>3.0.CO;2-#](https://doi.org/10.1002/%28SICI%291097-0258%2819970415%2916%3A7%3C791%3A%3AAID-SIM500%3E3.0.CO%3B2-%23)

## Examples

``` r
d <- data.frame(age = rep(c("18-29", "30-49", "50+"), 2),
                grp = rep(c("A", "B"), each = 3),
                x = c(30, 22, 5, 41, 30, 12),
                n = c(5000, 9000, 6000, 3000, 8000, 7000))
dastd(d$x, d$n, d$age, d$grp, std = c(30, 45, 25), per = 1e4)
#> dastd: standardised rates per 10000
#>  group  x     rate       lo       hi      var
#>      A 57 31.08333 23.49743 40.39982 17.16806
#>      B 83 62.16071 48.83268 78.21030 52.02280
#> against A:
#>  group    ratio       lo      hi     diff     dlo      dhi
#>      B 1.999809 1.414352 2.82761 31.07738 14.7742 47.38057
```
