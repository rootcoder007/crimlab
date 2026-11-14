# Indices of residential segregation

The evenness of groups across areal units (Massey and Denton 1988), from
a table of counts with one row per unit and one column per group:
Duncan's dissimilarity `D = sum |x_ik / X_k - x_il / X_l| / 2` between
each pair of groups (the share of either group that would have to move
for the two to be spread alike), and for each group against the rest the
segregation index `IS`, Theil's entropy index `H`, Atkinson's index with
shape `delta` and the Gini index. With weights `w` between the units,
Morrill's (1991) adjusted `D` subtracts the weighted mean absolute
difference of the pair's shares between neighbouring units, so that
segregation across a boundary counts less than segregation within. Units
with no one in them are left out. That `D` lies between 0 and 1 is
proved in Lean: see `proof("dissix")`.

## Usage

``` r
dissix(x, delta = 0.5, w = NULL)
```

## Arguments

- x:

  The counts: a matrix or data frame of numbers from 0, units by groups
  (named columns), every group with a positive total.

- delta:

  Atkinson's shape, above 0 and below 1.

- w:

  `NULL`, or weights between the units (as
  [`gmoran()`](https://rootcoder007.github.io/crimlab/reference/gmoran.md)
  takes them), with a positive total.

## Value

A list of class `dissix`: `d`, the matrix of `D`; `grp`, a data frame of
`group`, `is`, `h`, `atk` and `gini`; `dm`, `NULL` or the matrix of
Morrill's `D`; and `delta`.

## References

Duncan, O. D. and Duncan, B. (1955). A methodological analysis of
segregation indexes. *American Sociological Review*, 20(2), 210-217.
[doi:10.2307/2088328](https://doi.org/10.2307/2088328)

Massey, D. S. and Denton, N. A. (1988). The dimensions of residential
segregation. *Social Forces*, 67(2), 281-315.
[doi:10.2307/2579183](https://doi.org/10.2307/2579183)

Morrill, R. L. (1991). On the measure of geographic segregation.
*Geography Research Forum*, 11, 25-36.
[doi:10.64166/p1b4ca43](https://doi.org/10.64166/p1b4ca43)

## Examples

``` r
x <- cbind(a = c(90, 80, 10, 5, 50), b = c(10, 20, 85, 95, 50))
dissix(x)
#> dissix: dissimilarity between groups
#>           a         b
#> a 0.0000000 0.6284779
#> b 0.6284779 0.0000000
#>  group        is         h       atk      gini
#>      a 0.6284779 0.4052672 0.5217079 0.7680033
#>      b 0.6284779 0.4052672 0.5217079 0.7680033
dissix(x, w = spwt(nb = list(2, c(1, 3), c(2, 4), c(3, 5), 4),
                   type = "contig"))$dm
#>           a         b
#> a 0.0000000 0.3134779
#> b 0.3134779 0.0000000
```
