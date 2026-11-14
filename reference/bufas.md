# Treated, comparison and excluded places around sites

The buffer design of place-based evaluations: a place (a block, an
address) within distance `r` of a site (a new facility, a treated
corner) is treated; one within `r` of a comparison point and farther
than `buf` from every site is a comparison; the rest are excluded, so
that spillover from the sites (displacement or a diffusion of benefits,
Bowers and Johnson 2003) does not reach the comparison. `buf` should be
several times `r`.

## Usage

``` r
bufas(x, y = NULL, site, comp, r, buf = 3 * r, lonlat = FALSE)
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

- site, comp:

  The sites and the comparison points, as coordinates in the form of `x`
  (longitude and latitude when `lonlat = TRUE`).

- r:

  The treatment radius: positive, in the unit of the coordinates (metres
  when `lonlat = TRUE`).

- buf:

  The distance from every site a comparison place keeps, from 0.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list of class `bufas`: `pts`, a data frame with one row per place:
`grp` (a factor: treated, comparison, excluded), `dsite` and `dcomp`
(the distances to the nearest site and comparison point); `n`, the
counts of the groups; `par`, `r` and `buf`.

## References

Bowers, K. J. and Johnson, S. D. (2003). Measuring the geographical
displacement and diffusion of benefit effects of crime prevention
activity. *Journal of Quantitative Criminology*, 19(3), 275-301.
[doi:10.1023/A:1024909009240](https://doi.org/10.1023/A%3A1024909009240)

## Examples

``` r
set.seed(4)
p <- cbind(x = runif(200, 0, 5000), y = runif(200, 0, 5000))
b <- bufas(p, site = cbind(c(1000, 3500), c(1200, 3000)),
           comp = cbind(c(2500, 4500), c(500, 4000)), r = 330)
b
#> bufas: radius 330, buffer 990
#>    treated comparison   excluded 
#>          5          4        191 
```
