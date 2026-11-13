# Local G and G\* of Getis and Ord

The local statistics of Getis and Ord (1992) as z values (Ord and Getis
1995): `G_i` compares the weighted sum of the neighbours' values with
its expectation from the other `n - 1` values; `G*_i` (`star = TRUE`)
includes the unit itself with weight 1 and uses all `n`. Large positive
values mark concentrations of high values, large negative ones of low
values. With `star = TRUE` the weights are best binary
(`style = "raw"`).

## Usage

``` r
locgi(x, w, star = TRUE)
```

## Arguments

- x:

  The values: finite numbers, one per unit.

- w:

  The weights:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, or an spdep `listw` object, for as many units as `x`.

- star:

  `TRUE` for G\*, `FALSE` for G.

## Value

A numeric vector of z values, one per unit; `NaN` for a unit without
neighbours.

## References

Ord, J. K. and Getis, A. (1995). Local spatial autocorrelation
statistics: distributional issues and an application. *Geographical
Analysis*, 27(4), 286-306.
[doi:10.1111/j.1538-4632.1995.tb00912.x](https://doi.org/10.1111/j.1538-4632.1995.tb00912.x)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(60), y = runif(60))
v <- rpois(60, ifelse(p[, 1] > 0.7, 6, 2))
locgi(v, spwt(p, type = "band", d = 0.2, style = "raw"))
#>  [1]  0.23994257 -1.19393804 -0.85693776  4.16902467 -1.51580855  2.62650766
#>  [7]  3.78952137  0.87555456 -0.91284462 -1.74706652 -2.08423675 -0.89182354
#> [13]  1.19393804 -1.51232152  1.27896346 -0.37705262  2.28529782  2.46318889
#> [19] -0.23878761  2.08423675  4.04305111 -1.05904046  0.55717109 -2.22803818
#> [25] -0.07959587 -0.72460663 -1.06343180 -0.56842821  2.08423675 -0.73778305
#> [31] -0.72460663 -0.91284462 -0.72460663 -0.75790427  2.30828021  1.28199634
#> [37]  1.78732733 -1.97494476  2.30828021 -0.56071512  2.09530386  1.19393804
#> [43]  1.03474630  0.30383766 -1.72790811  1.03289627 -1.17365736 -1.62693646
#> [49]  1.44921325  0.78034560 -0.65203187  2.50225827 -1.44605478 -2.08423675
#> [55] -1.51918828 -1.17365736 -1.03474630 -0.91485098 -0.23878761 -0.55738971
```
