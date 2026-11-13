# Indirect standardisation: standardised ratios with shrinkage

The events observed in each group (an area, a police service) against
those expected from stratum rates, the standardised incidence ratio
(SIR) with its exact Poisson interval, funnel-plot limits for a ratio of
1 at its expected count, and the ratio shrunk towards the overall one by
the empirical Bayes estimate of Marshall (a gamma prior fitted by
moments), so that groups with little exposure are not ranked on noise.

## Usage

``` r
sirex(x, n, strat, grp, rates = NULL, level = 0.95)
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

  The group of each row: labels.

- rates:

  `NULL` for the stratum rates pooled over all groups, or the reference
  rate of each stratum: numbers from 0, named by stratum or in the order
  of the strata.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `sirex`, a list: `sir`, a data frame (`group`, `obs`,
`expd`, `sir`, `lo`, `hi`, `eb` the shrunk ratio, `shr` the share of the
prior in it, `flo`, `fhi` the funnel limits, `out` `TRUE` outside them);
`prior`, the prior's mean and variance (a variance of 0 when the ratios
vary no more than chance, and then every `eb` is the mean); and `level`.

## References

Marshall, R. J. (1991). Mapping disease and mortality rates using
empirical Bayes estimators. *Journal of the Royal Statistical Society,
Series C*, 40(2), 283-294.
[doi:10.2307/2347593](https://doi.org/10.2307/2347593)

## Examples

``` r
d <- data.frame(area = rep(c("N", "S", "E"), each = 2),
                age = rep(c("young", "old"), 3),
                x = c(12, 3, 30, 9, 4, 1),
                n = c(4000, 2000, 8000, 3000, 1500, 900))
sirex(d$x, d$n, d$age, d$area)
#> sirex: observed over expected, 95% intervals and funnel limits
#>  group obs      expd       sir        lo       hi eb shr       flo      fhi
#>      E   5  7.094162 0.7048049 0.2288482 1.644779  1   1 0.2819220 1.832493
#>      N  15 18.036409 0.8316511 0.4654688 1.371682  1   1 0.5544341 1.496972
#>      S  39 33.869429 1.1514809 0.8188155 1.574112  1   1 0.6790785 1.358157
#>    out
#>  FALSE
#>  FALSE
#>  FALSE
```
