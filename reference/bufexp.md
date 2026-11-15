# Exposure within a buffer

For each receptor (a home, a school), the number of sources (bars,
polluting sites, hot spots) within distance `r`, their summed weights
and the distance to the nearest: the buffer measure of exposure
assessment. Distances are planar, or great-circle metres when
`lonlat = TRUE` (the coordinates then longitude and latitude in
degrees).

## Usage

``` r
bufexp(x, y = NULL, src, r, w = NULL, lonlat = FALSE)
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

- src:

  The sources: a matrix or data frame of two coordinate columns, as `x`.

- r:

  The radius: one number from 0, in the unit of the coordinates (metres
  when `lonlat = TRUE`).

- w:

  `NULL` for a weight of 1 per source, or one weight from 0 per source.

- lonlat:

  `TRUE` when the receptors and sources are longitude and latitude in
  degrees.

## Value

A data frame with one row per receptor: `count`, `wsum` and `near`.

## References

Nuckols, J. R., Ward, M. H. and Jarup, L. (2004). Using geographic
information systems for exposure assessment in environmental
epidemiology studies. *Environmental Health Perspectives*, 112(9),
1007-1015. [doi:10.1289/ehp.6738](https://doi.org/10.1289/ehp.6738)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
rec <- cbind(x = c(0, 1, 5), y = c(0, 1, 5))
src <- cbind(x = c(0.5, 4), y = c(0.5, 4))
bufexp(rec, src = src, r = 1.5)
#>   count wsum      near
#> 1     1    1 0.7071068
#> 2     1    1 0.7071068
#> 3     1    1 1.4142136
```
