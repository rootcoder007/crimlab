# Crosswalk of counts between two sets of areal units

Moves values (counts, populations) kept on one set of units (an old
boundary version, say) onto another, through the pairs of a source and a
target unit that overlap: each source's value is split among its targets
in proportion to `area`, the size of the overlap (its area, or its
population for a population-weighted crosswalk), which assumes the value
spread uniformly over the source. Where each source lies in one target
(a boundary version that only merges units) the result is exact
aggregation and `exact` is `TRUE`; splitting a source is an
approximation. A source with no pair is lost (`lost`), and each source's
main target (its largest share) maps an old identifier onto the new
scheme. The concentration of the result across the targets is reported
as the Gini coefficient and the share of the largest `top` targets.

## Usage

``` r
arwt(x, from, to, area = NULL, id = names(x), top = c(10, 20))
```

## Arguments

- x:

  The values of the source units: finite numbers.

- from, to:

  The overlapping pairs: the source unit (a label in `id`) and the
  target unit of each, as many as there are pairs.

- area:

  `NULL` to split each source evenly among its targets, or the size of
  each overlap, from 0.

- id:

  The labels of the source units, one per value of `x`, each once; by
  default `names(x)`.

- top:

  The numbers of largest targets whose share is reported: whole numbers
  from 1.

## Value

A list of class `arwt`: `out`, a data frame of `to` and `value` (every
target in a pair, in sorted order); `map`, a data frame of `from`, `to`
(its main target) and `share`; `exact`; `lost`; and `conc`, a named
vector: `gini` and `top<k>`.

## References

Lam, N. S.-N. (1983). Spatial interpolation methods: a review. *The
American Cartographer*, 10(2), 129-150.
[doi:10.1559/152304083783914958](https://doi.org/10.1559/152304083783914958)

## Examples

``` r
# four old units onto three new ones; unit d is split 70/30 by area
arwt(c(a = 120, b = 40, c = 60, d = 100),
     from = c("a", "b", "c", "d", "d"), to = c("N1", "N1", "N2", "N2", "N3"),
     area = c(1, 1, 1, 7, 3), top = 1)
#> arwt: 4 sources onto 3 targets, sources split by overlap
#>  to value
#>  N1   160
#>  N2   130
#>  N3    30
#>   gini   top1 
#> 0.2708 0.5000 
```
