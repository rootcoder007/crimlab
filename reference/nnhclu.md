# Nearest-neighbour hierarchical clustering

Groups points that lie closer together than a threshold distance, then
groups the groups, order by order (after Johnson 1967, as used for crime
hot spots). At the first order two points are linked when their distance
is below the threshold: `fixd`, or the mean nearest-neighbour distance
expected under complete spatial randomness in an area `area` plus `z`
standard errors, `0.5 sqrt(A / m) + z sqrt((4 - pi) A / (4 pi)) / m` for
`m` units, with `z = qnorm(p)` (`p = 0.5` gives the expectation itself,
smaller `p` a stricter threshold). Points with a linked neighbour are
taken in decreasing number of links; each one not yet clustered takes
every unclustered point linked to it; groups of at least `minpts` points
are kept. Their centres of minimum distance become the units of the next
order, with the threshold recomputed for their number, while at least
four remain, up to `maxord` orders. The iterative reassignment some
descriptions add after the first pass is not done: its rule is not
stated unambiguously.

## Usage

``` r
nnhclu(
  x,
  y = NULL,
  area = NULL,
  fixd = NULL,
  p = 0.5,
  minpts = 10,
  minhi = minpts,
  maxord = 3,
  k = 1,
  c = 1,
  count = NULL,
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

- area:

  The study area for the random-expectation threshold, in the square of
  the coordinates' unit (square metres when `lonlat = TRUE`). Give
  `area` or `fixd`, not both.

- fixd:

  A fixed threshold distance instead, positive.

- p:

  The one-tailed probability that a pair is as close by chance: one
  number strictly between 0 and 1.

- minpts, minhi:

  The fewest points (order 1) and member units (higher orders) a cluster
  keeps: whole numbers from 2.

- maxord:

  The highest order: one whole number from 1 to 20.

- k:

  The divisor offset: 0, 1 or 2; the total weight must exceed it.

- c:

  The multiplier of each cluster's standard deviational ellipse (see
  [`sdell()`](https://rootcoder007.github.io/crimlab/reference/sdell.md)).

- count:

  `NULL` for points, or for zones: one whole count from 0 per zone, the
  points `x`, `y` being the zones' centroids. A zone is taken as `count`
  incidents at its centroid, which these hot-spot rules treat exactly as
  coincident points. Zones with count 0 get `NA`; the random-expectation
  threshold counts incidents, not zones.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list. `cl`, an integer matrix with one row per point (input order) and
one column per order: the cluster (a row of `clusters`) the point
belongs to at that order, or `NA`. `clusters`, a data frame with one row
per cluster: `order`, `points`, `units` (member units), `mx`, `my` (mean
centre), `wx`, `wy` (centre of minimum distance; both in degrees when
`lonlat = TRUE`), `a`, `b`, `angle` (the ellipse of the member units,
`NaN` when too few for the divisor) and `thresh` (the threshold distance
used).

## References

Johnson, S. C. (1967). Hierarchical clustering schemes. *Psychometrika*,
32(3), 241-254.
[doi:10.1007/BF02289588](https://doi.org/10.1007/BF02289588)

Eck, J. E., Chainey, S., Cameron, J. G., Leitner, M. and Wilson, R. E.
(2005). *Mapping Crime: Understanding Hot Spots*. National Institute of
Justice, NCJ 209393.

## Examples

``` r
set.seed(1)
x <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
y <- c(rnorm(40, 0, 0.3), rnorm(40, 6, 0.3), runif(40, -5, 10))
nnhclu(x, y, area = 225, p = 0.05, minpts = 5)$clusters
#>   order points units         mx           my         wx          wy         a
#> 1     1     40    40 0.03845448 -0.003796574 0.06246601 -0.03941408 0.2470327
#> 2     1     36    36 6.00748914  5.922730801 5.99404607  5.93418792 0.3063012
#> 3     1      5     5 6.13610851  6.711097797 6.09427553  6.66090240 0.1958946
#>           b    angle    thresh
#> 1 0.2266867 2.820476 0.6309155
#> 2 0.2568727 2.060319 0.6309155
#> 3 0.0298542 0.769340 0.6309155
```
