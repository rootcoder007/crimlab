# Log-likelihood-ratio cost after recalibration

Splits the cost of a set of likelihood ratios into `min`, their cost
after recalibration with the pool-adjacent-violators algorithm, which
keeps their order and gives each block of pooled ratios the value that
costs that block least, and `cal`, the rest (`cllr - min`). That `min`
is the least cost of every order-keeping recalibration, so that `cal` is
never negative, is not proved. A system can separate the two kinds of
comparison well and still report ratios that are too large or too small.
See `proof("cllrmin")` and `proof("lrcal")`.

## Usage

``` r
cllrmin(same, diff)
```

## Arguments

- same:

  Likelihood ratios from same-source comparisons: positive numbers.

- diff:

  Likelihood ratios from different-source comparisons: positive numbers.

## Value

A named numeric vector: `cllr`, `min` and `cal`, with
`cllr = min + cal`.

## References

Brummer, N. and du Preez, J. (2006). Application-independent evaluation
of speaker detection. *Computer Speech and Language*, 20(2-3), 230-275.
[doi:10.1016/j.csl.2005.08.001](https://doi.org/10.1016/j.csl.2005.08.001)

## Examples

``` r
cllrmin(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#>      cllr       min       cal 
#> 0.3008346 0.0000000 0.3008346 
```
