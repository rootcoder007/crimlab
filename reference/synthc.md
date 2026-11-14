# Synthetic control for one treated unit

A weighted average of units never treated that tracks the treated unit
before the policy, as its counterfactual after; with the placebo test
that runs the same method on every comparison unit.

## Usage

``` r
synthc(y, unit, time, first, plac = TRUE)
```

## Arguments

- y:

  The outcome: finite numbers, one per row (counts from 0 for
  `fam = "pois"` or `"nb"`).

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; a balanced panel.

- first:

  The period the treated unit is first treated, on its rows; `NA`, 0 or
  `Inf` on the rows of the comparison units. Exactly one unit is
  treated, from a period after the first.

- plac:

  `TRUE` to run the placebo test.

## Value

An object of class `synthc`, a list: `weights` (`unit`, `w`); `path`
(`time`, `treated`, `synth`, `gap`); `fit`, one row: `rmspepre`,
`rmspepost`, their `ratio` and the placebo `p` (`NA` without the
placebos); `plac` (`unit`, `ratio`; `NULL` without); and `unit`, the
treated unit.

## Details

The weights solve `min_w sum_{t < T0} (y_1t - sum_j w_j y_jt)^2` over
`w_j >= 0`, `sum_j w_j = 1`, the periods before adoption `T0`. The
problem is the nearest point of a convex hull to the origin, solved
exactly by Wolfe's minimum-norm-point algorithm. The effect in each
period is the gap `y_1t - sum_j w_j y_jt`. For inference each comparison
unit in turn is treated as if it had adopted, with the others as its
donors; the ratio of the root mean squared gap after adoption to that
before ranks the treated unit among the placebos, and `p` is
`(1 + the placebos with a ratio at least the treated unit's) / (1 + the placebos)`
(Abadie, Diamond and Hainmueller 2010). The weights match the outcomes
before adoption only, not separate predictors.

## References

Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control
methods for comparative case studies: estimating the effect of
California's tobacco control program. *Journal of the American
Statistical Association*, 105(490), 493-505.
[doi:10.1198/jasa.2009.ap08746](https://doi.org/10.1198/jasa.2009.ap08746)

Wolfe, P. (1976). Finding the nearest point in a polytope. *Mathematical
Programming*, 11, 128-149.
[doi:10.1007/BF01580381](https://doi.org/10.1007/BF01580381)

## Examples

``` r
d <- expand.grid(year = 2010:2024, town = LETTERS[1:8])
i <- match(d$town, LETTERS)
d$first <- ifelse(d$town == "A", 2019, NA)
d$rate <- 20 + i + (d$year - 2010) * (0.2 + i / 20) +
  sin(seq_len(120)) - 3 * (d$town == "A" & d$year >= 2019)
s <- synthc(d$rate, d$town, d$year, d$first)
s
#> synthc: unit A, root mean squared gap 1.513816 before and 4.813741 after; placebo p 0.125
#>  unit          w
#>     B 0.97743517
#>     C 0.02256483
s$path
#>    time  treated    synth        gap
#> 1  2010 21.84147 21.73204  0.1094300
#> 2  2011 22.15930 21.39643  0.7628652
#> 3  2012 21.64112 21.91334 -0.2722228
#> 4  2013 20.99320 23.08438 -2.0911860
#> 5  2014 21.04108 24.10976 -3.0686850
#> 6  2015 21.97058 24.32360 -2.3530189
#> 7  2016 23.15699 23.80616 -0.6491748
#> 8  2017 23.73936 23.31002  0.4293338
#> 9  2018 23.41212 23.56819 -0.1560760
#> 10 2019 19.70598 24.62017 -4.9141882
#> 11 2020 19.50001 25.77562 -6.2756096
#> 12 2021 20.21343 26.24909 -6.0356627
#> 13 2022 21.42017 25.88213 -4.4619607
#> 14 2023 22.24061 25.28897 -3.0483649
#> 15 2024 22.15029 25.29182 -3.1415360
```
