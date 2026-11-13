# Head-bang smoothing

Smooths values at points or areal units (rates by area centroid, say) by
Hansen's (1991) head-bang rule, which removes isolated peaks and troughs
while keeping edges between regions of different level. For each unit
`i` the rule takes its `nn` nearest neighbours and every pair `(j, k)`
of them whose angle at `i` is at least `ang` degrees, so that `i` lies
roughly between `j` and `k` (a triple); with `ntrip > 0` only the
`ntrip` triples whose ends are nearest (smallest
`|p_j - p_i| + |p_k - p_i|`) are kept. Each triple gives a low value
`min(v_j, v_k)` and a high value `max(v_j, v_k)`; the low screen is the
median of the lows, the high screen the median of the highs, and the new
value is the median of the low screen, `v_i` and the high screen, so
`v_i` moves only if it lies outside both. With weights `w` (populations
at risk, say) the screens are weighted medians, each end weighted by its
own unit's weight, after Mungiole, Pickle and Simonson (1999): the
smallest value whose cumulative weight reaches half the total. A unit
with no triple (on the outer edge of the points) keeps its value;
Hansen's extrapolated end triples are not formed. Each of `iter` passes
starts from the previous pass's values.

## Usage

``` r
hdbang(
  x,
  y = NULL,
  v,
  w = NULL,
  nn = 8,
  ang = 135,
  ntrip = 0,
  iter = 1,
  lonlat = FALSE
)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  names them, or `lon` and `lat` with `lonlat = TRUE`) and `y` `NULL`;
  or `x` an sf or sfc object of POINT geometries, whose coordinate
  reference system must agree with `lonlat`. An unnamed two-column
  matrix is read as `x` then `y`; other columns are never taken by
  position. Vectors given separately (`x`, `y`, and `w` or `t` where a
  function takes them) are matched by name when every one is named and
  no name repeats, as elsewhere in the package; otherwise by position.
  No missing or infinite values.

- v:

  The values to smooth: finite numbers, one per point.

- w:

  `NULL`, or positive weights, one per point.

- nn:

  The number of nearest neighbours searched: from 2 to the number of
  points less 1.

- ang:

  The smallest angle at `i` of a triple, in degrees, from 90 to 180 (135
  by convention).

- ntrip:

  The most triples used per unit, nearest first; 0 for all.

- iter:

  The number of passes: a whole number from 0.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A numeric vector of smoothed values, one per point, in input order.

## References

Hansen, K. M. (1991). Head-banging: robust smoothing in the plane. *IEEE
Transactions on Geoscience and Remote Sensing*, 29(3), 369-378.

Mungiole, M., Pickle, L. W. and Simonson, K. H. (1999). Application of a
weighted head-banging algorithm to mortality data maps. *Statistics in
Medicine*, 18(23), 3201-3209.

## Examples

``` r
g <- expand.grid(x = 1:6, y = 1:6)
v <- rep(1, 36)
v[15] <- 50
hdbang(g, v = v)[15]
#> [1] 1
```
