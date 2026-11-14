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

## References

Ayer, M., Brunk, H. D., Ewing, G. M., Reid, W. T. and Silverman, E.
(1955). An empirical distribution function for sampling with incomplete
information. *The Annals of Mathematical Statistics*, 26(4), 641-647.
[doi:10.1214/aoms/1177728423](https://doi.org/10.1214/aoms/1177728423)

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
