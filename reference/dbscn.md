# DBSCAN clusters of incidents

Density-based clusters (Ester, Kriegel, Sander and Xu 1996): a core
point has at least `minpts` points within distance `eps`, itself
included; core points within `eps` of each other share a cluster, and a
point that is not core but lies within `eps` of a core point (a border
point) joins the first cluster that reaches it, clusters being grown
from the core points in input order. The rest is noise. With
`border = FALSE` border points are noise too (DBSCAN\*, Campello,
Moulavi and Sander 2013), which makes the result independent of the
order. These are the clusters of `dbscan::dbscan()`, numbered the same
way.

## Usage

``` r
dbscn(x, y = NULL, eps, minpts = 5, border = TRUE, lonlat = FALSE)
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

- eps:

  The neighbourhood radius: positive, in the unit of the coordinates
  (metres when `lonlat = TRUE`).

- minpts:

  The fewest points (the point itself included) within `eps` of a core
  point: a whole number from 1.

- border:

  `TRUE` to put border points in clusters.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list of class `dbscn`: `cl` (each point's cluster from 1, `NA` for
noise), `core` (whether each point is a core point) and `sizes` (the
points in each cluster).

## References

Ester, M., Kriegel, H.-P., Sander, J. and Xu, X. (1996). A density-based
algorithm for discovering clusters in large spatial databases with
noise. *Proceedings of the Second International Conference on Knowledge
Discovery and Data Mining*, 226-231.

Campello, R. J. G. B., Moulavi, D. and Sander, J. (2013). Density-based
clustering based on hierarchical density estimates. *Advances in
Knowledge Discovery and Data Mining*, 160-172.
[doi:10.1007/978-3-642-37456-2_14](https://doi.org/10.1007/978-3-642-37456-2_14)

## Examples

``` r
set.seed(1)
p <- rbind(cbind(x = rnorm(30, 0, 0.2), y = rnorm(30, 0, 0.2)),
           cbind(x = rnorm(30, 3, 0.2), y = rnorm(30, 1, 0.2)),
           cbind(x = runif(10, -2, 5), y = runif(10, -2, 3)))
d <- dbscn(p, eps = 0.3, minpts = 5)
d
#> dbscn: 2 clusters of 70 points, 10 noise
#>  cluster points
#>        1     30
#>        2     30
table(d$cl, useNA = "ifany")
#> 
#>    1    2 <NA> 
#>   30   30   10 
```
