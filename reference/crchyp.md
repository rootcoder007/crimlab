# Circle hypothesis of a crime series

Canter and Larkin's (1993) circle hypothesis: draw the circle whose
diameter joins the two crime sites of a series that lie farthest apart.
An offender whose home lies within it is a "marauder" (who travels out
from home), one whose home lies outside a "commuter". Gives the circle,
the share of the sites it holds (all of them, unless the series is far
from round) and, for each home given, whether it lies within.

## Usage

``` r
crchyp(x, y, home = NULL)
```

## Arguments

- x, y:

  The coordinates of the crime sites: finite numbers, at least two
  sites, in a planar projection (metres, feet).

- home:

  `NULL`, or the home of the offender: `c(x, y)`, or a matrix of two
  columns with one row per candidate home.

## Value

A list: `centre` (`x` and `y`), `r` the radius, `pair` the two sites
(row numbers) on the diameter, `share` the share of the sites within the
circle, and `mar`, `NULL` or for each home `TRUE` when it lies within (a
marauder).

## Details

Of several pairs equally far apart, the first in the order of the sites
is used. A point lies within the circle when its distance from the
centre is at most \\r + 10^{-12} \max(1, r)\\.

## References

Canter, D. and Larkin, P. (1993). The environmental range of serial
rapists. *Journal of Environmental Psychology*, 13(1), 63-69.
[doi:10.1016/S0272-4944(05)80215-4](https://doi.org/10.1016/S0272-4944%2805%2980215-4)

## Examples

``` r
crchyp(c(0, 4, 2), c(0, 0, 1), home = c(2, -1))
#> $centre
#> x y 
#> 2 0 
#> 
#> $r
#> [1] 2
#> 
#> $pair
#> [1] 1 2
#> 
#> $share
#> [1] 1
#> 
#> $mar
#> [1] TRUE
#> 
```
