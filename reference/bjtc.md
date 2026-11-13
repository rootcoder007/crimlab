# Bayesian journey to crime

Combines a journey-to-crime surface over origin zones with an
origin-destination matrix of solved offences (Levine 2009; Levine and
Block 2011): `P(Jtc)`, the surface; `P(O)`, the origins of all solved
offences; `P(O|Jtc)`, the origins of offences in the zones where the
series happened; their product, Bayes' numerator with the surface as
prior; and the risk, the product over `P(O)`, renormalised. Each block
sums to 1. The matrix describes past solved offences and inherits their
selection.

## Usage

``` r
bjtc(pj, od, dest)
```

## Arguments

- pj:

  The journey-to-crime value of each origin zone: numbers from 0, not
  all 0.

- od:

  The origin-destination counts: a matrix with one row per origin zone
  and one column per destination zone, numbers from 0.

- dest:

  The destination zone (column of `od`) of each incident of the series:
  whole numbers from 1; repeats count.

## Value

A data frame with one row per origin zone: `pj`, `po`, `pod`, `prod` and
`risk` (0 where `po` is 0).

## References

Levine, N. and Block, R. (2011). Bayesian journey-to-crime estimation:
an improvement in geographic profiling methodology. *The Professional
Geographer*, 63(2), 213-229.
[doi:10.1080/00330124.2010.547152](https://doi.org/10.1080/00330124.2010.547152)

## Examples

``` r
od <- matrix(c(5, 1, 0, 2, 4, 1, 0, 1, 6), 3)
bjtc(pj = c(0.5, 0.3, 0.2), od = od, dest = c(1, 2, 2))
#>    pj   po  pod       prod       risk
#> 1 0.5 0.35 0.45 0.59210526 0.55900621
#> 2 0.3 0.30 0.45 0.35526316 0.39130435
#> 3 0.2 0.35 0.10 0.05263158 0.04968944
```
