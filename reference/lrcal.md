# Calibrated likelihood ratios

Recalibrates a set of likelihood ratios against the truth with the
pool-adjacent-violators algorithm, keeping their order. The result can
hold 0 and `Inf`: the data never contradicted those ratios. See
`proof("lrcal")`.

## Usage

``` r
lrcal(same, diff)
```

## Arguments

- same:

  Likelihood ratios from same-source comparisons: positive numbers.

- diff:

  Likelihood ratios from different-source comparisons: positive numbers.

## Value

A list: `same` and `diff`, the calibrated ratios in the order given.

## Examples

``` r
lrcal(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#> $same
#> [1] Inf Inf Inf
#> 
#> $diff
#> [1] 0 0 0
#> 
```
