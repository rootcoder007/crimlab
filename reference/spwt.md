# Spatial weights

Builds the neighbours and weights the autocorrelation statistics
([`gearc()`](https://rootcoder007.github.io/crimlab/reference/gearc.md),
[`gorg()`](https://rootcoder007.github.io/crimlab/reference/gorg.md),
[`lmorn()`](https://rootcoder007.github.io/crimlab/reference/lmorn.md),
[`locgi()`](https://rootcoder007.github.io/crimlab/reference/locgi.md))
use: the `k` nearest neighbours of each point (`type = "knn"`, ties
broken by the lower row number), the points within distance `d`
(`"band"`, the distance itself included), or a contiguity list `nb`
(`"contig"`, one element per unit giving its neighbours' numbers, as
spdep's `poly2nb()` returns). Weights are 1, or the inverse distance to
the power `idw`; `style = "row"` divides each unit's weights by their
sum (a unit without neighbours keeps none), `"raw"` leaves them.

## Usage

``` r
spwt(
  x = NULL,
  y = NULL,
  type = "knn",
  k = 4,
  d = NULL,
  nb = NULL,
  idw = 0,
  style = "row",
  metric = "euclid",
  lonlat = FALSE
)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as `sf::st_coordinates()` names them, or `lon` and `lat` with
  `lonlat = TRUE`) and `y` `NULL`; or `x` an sf or sfc object of POINT
  geometries, whose coordinate reference system must agree with
  `lonlat`. An unnamed two-column matrix is read as `x` then `y`; other
  columns are never taken by position. Vectors given separately (`x`,
  `y`, and `w` or `t` where a function takes them) are matched by name
  when every one is named and no name repeats, as elsewhere in the
  package; otherwise by position. No missing or infinite values.

- type:

  `"knn"`, `"band"` or `"contig"`.

- k:

  For `"knn"`: the number of neighbours, from 1 to the number of points
  less 1.

- d:

  For `"band"`: the distance, positive, in the unit of the coordinates
  (metres when `lonlat = TRUE`).

- nb:

  For `"contig"`: a list with one element per unit, each a vector of
  neighbour numbers (from 1; `0` alone for none).

- idw:

  The power of the inverse distance: 0 (binary weights) or more; above 0
  it needs the points.

- style:

  `"row"` or `"raw"`.

- metric:

  `"euclid"` or `"manhattan"` (`|dx| + |dy|`).

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list: `n` (the number of units), and one element per neighbour pair,
in order of unit then neighbour: `i`, `j` (integer unit numbers), `d`
(their distance, `NA` for contiguity without points) and `w` (the
weight). `x` may be `NULL` for contiguity with binary weights.

## References

Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models and
Applications*. Pion, London.

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(20), y = runif(20))
w <- spwt(p, k = 3)
str(w)
#> List of 5
#>  $ n: int 20
#>  $ i: int [1:60] 1 1 1 2 2 2 3 3 3 4 ...
#>  $ j: int [1:60] 12 16 19 5 11 14 13 16 17 6 ...
#>  $ d: num [1:60] 0.347 0.353 0.24 0.179 0.317 ...
#>  $ w: num [1:60] 0.333 0.333 0.333 0.333 0.333 ...
spwt(nb = list(2, c(1, 3), 2), type = "contig")$w
#> [1] 1.0 0.5 0.5 1.0
```
