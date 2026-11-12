# Regression to the mean at selected hot spots

Places chosen because their count in a first period exceeds a cut show a
fall in the next period even with no treatment and no trend: if the two
periods are exchangeable, the selected places' mean change is never
positive, and places chosen for being low never fall on average. This
reports the change at the selected places next to the same statistic on
the data together with its period-swapped copy, which is exchangeable by
construction and so never positive; the excess of the first over the
second is the part of the fall that selection alone does not explain. A
treatment effect needs a control group selected the same way. See
`proof("rtmhs")`.

## Usage

``` r
rtmhs(x1, x2, cut, w = NULL)
```

## Arguments

- x1:

  Counts by place in the first period: finite numbers.

- x2:

  Counts by place in the second period: same length.

- cut:

  The selection cut: places with `x1 > cut` are selected.

- w:

  `NULL` for equal weights, or one non-negative weight per place.

## Value

A named numeric vector: `nsel` (the weight selected), `change` (the
weighted mean of `x2 - x1` at the selected places), `mirror` (the mean
of `x1 - x2` at places with `x2 > cut`), `sym` (the change in the data
with its swapped copy; never positive), `excess` (`change - sym`), `low`
(the mean change at places with `x1 < cut`) and `lowsym` (its
swapped-copy version; never negative). `NaN` where no place qualifies.

## References

Barnett, A. G., van der Pols, J. C. and Dobson, A. J. (2005). Regression
to the mean: what it is and how to deal with it. *International Journal
of Epidemiology*, 34(1), 215-220.
[doi:10.1093/ije/dyh299](https://doi.org/10.1093/ije/dyh299)

## Examples

``` r
set.seed(5)
mu <- rgamma(300, 2, 0.5)
x1 <- rpois(300, mu)
x2 <- rpois(300, mu)
rtmhs(x1, x2, cut = quantile(x1, 0.9))
#>       nsel     change     mirror        sym     excess        low     lowsym 
#> 27.0000000 -2.5555556 -1.8400000 -2.2115385 -0.3440171  0.1492537  0.2640449 
```
