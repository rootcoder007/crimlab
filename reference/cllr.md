# Log-likelihood-ratio cost

Scores a set of likelihood ratios against the truth. Ratios from
same-source comparisons should be large and ratios from different-source
comparisons small. The cost is 1 for a system that always answers 1. It
is above 0 for any set of ratios, and tends to 0 as the ratios of a
system that separates the two kinds perfectly grow without bound. See
`proof("cllr")` and
[`cllrmin()`](https://rootcoder007.github.io/crimlab/reference/cllrmin.md).

## Usage

``` r
cllr(same, diff)
```

## Arguments

- same:

  Likelihood ratios from same-source comparisons: positive numbers.

- diff:

  Likelihood ratios from different-source comparisons: positive numbers.

## Value

One number, the cost in bits.

## References

Brummer, N. and du Preez, J. (2006). Application-independent evaluation
of speaker detection. *Computer Speech and Language*, 20(2-3), 230-275.
[doi:10.1016/j.csl.2005.08.001](https://doi.org/10.1016/j.csl.2005.08.001)

## Examples

``` r
cllr(same = c(50, 200, 8), diff = c(0.01, 0.2, 1.5))
#> [1] 0.3008346
```
