# Rosenbaum sensitivity bounds for matched pairs

How strong would hidden bias have to be to explain away a matched-pairs
result? `gam` is Gamma, the most the odds of treatment can differ within
a pair because of an unobserved covariate (1: no hidden bias). `rswp()`
gives the bounds on the one-sided p-value of Wilcoxon's signed-rank test
at each Gamma (Rosenbaum 1987), from the normal approximation, and the
Gamma at which the upper bound reaches `alpha`. `rshl()` gives the
bounds on the Hodges-Lehmann estimate of the effect (Rosenbaum 1993),
computed exactly from the Walsh averages, with no grid.

## Usage

``` r
rswp(d, gam = 1, alpha = 0.05)

rshl(d, gam = 1)
```

## Arguments

- d:

  The treated-minus-control differences: finite numbers, one per pair.

- gam:

  Values of Gamma: numbers from 1.

- alpha:

  The level for the critical Gamma, strictly between 0 and 0.5.

## Value

`rswp()`: a list of `p` (a data frame of `gam`, `plo`, `pup`), `gcrit`
(the critical Gamma, below 1 when the test is not significant without
bias) and `t` (the signed-rank statistic). `rshl()`: a data frame of
`gam`, `lo`, `hi`. `NaN` when every difference is 0.

## References

Rosenbaum, P. R. (1987). Sensitivity analysis for certain permutation
inferences in matched observational studies. *Biometrika*, 74(1), 13-26.
[doi:10.1093/biomet/74.1.13](https://doi.org/10.1093/biomet/74.1.13)

Rosenbaum, P. R. (1993). Hodges-Lehmann point estimates of treatment
effect in observational studies. *Journal of the American Statistical
Association*, 88(424), 1250-1253.
[doi:10.1080/01621459.1993.10476405](https://doi.org/10.1080/01621459.1993.10476405)

## Examples

``` r
set.seed(1)
d <- rnorm(50, 0.5)
rswp(d, gam = c(1, 1.5, 2))
#> $p
#>   gam          plo          pup
#> 1 1.0 7.807040e-06 7.807040e-06
#> 2 1.5 7.347225e-09 8.087079e-04
#> 3 2.0 7.012883e-12 8.061278e-03
#> 
#> $gcrit
#> [1] 2.741017
#> 
#> $t
#> [1] 1085
#> 
rshl(d, gam = c(1, 1.5, 2))
#>   gam        lo        hi
#> 1 1.0 0.6526343 0.6526343
#> 2 1.5 0.4985346 0.7876529
#> 3 2.0 0.3952356 0.8803900
```
