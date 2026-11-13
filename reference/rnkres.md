# Which pairs a noisy ranking orders

People (or areas, or courts) ranked by estimated scores are ordered
correctly whatever the noise only when their estimates differ by more
than their two interval half-widths together: below that, admissible
errors can reverse them, above it none can. `rnkres()` says which pairs
that holds for, and the share it fails for; the resolution of the
ranking is twice the largest half-width. A league table whose intervals
are wide orders few of its pairs (Goldstein and Spiegelhalter 1996). See
`proof("rnkres")`.

## Usage

``` r
rnkres(est, hw)
```

## Arguments

- est:

  The estimates: finite numbers, named or not, at least two.

- hw:

  The half-width of each estimate's interval: one number from 0, or one
  per estimate.

## Value

A list: `id`, a logical matrix, `TRUE` where the pair's order is
identified (`NA` on the diagonal), with the estimates' names; `npair`,
`nunid` (the pairs not identified), `share` (their share) and `res`.

## References

Goldstein, H. and Spiegelhalter, D. J. (1996). League tables and their
limitations: statistical issues in comparisons of institutional
performance. *Journal of the Royal Statistical Society, Series A*,
159(3), 385-443. [doi:10.2307/2983325](https://doi.org/10.2307/2983325)

## Examples

``` r
rnkres(c(a = 0.2, b = 0.35, c = 0.8), hw = c(0.1, 0.1, 0.05))
#> $id
#>       a     b    c
#> a    NA FALSE TRUE
#> b FALSE    NA TRUE
#> c  TRUE  TRUE   NA
#> 
#> $npair
#> [1] 3
#> 
#> $nunid
#> [1] 1
#> 
#> $share
#> [1] 0.3333333
#> 
#> $res
#> [1] 0.2
#> 
```
