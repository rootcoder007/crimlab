# Standard deviational ellipse

The ellipse whose axes lie along the eigenvectors of the weighted
covariance matrix of the points about their mean centre, with semi-axes
`c` times the square roots of its eigenvalues (Wang, Shi and Miao 2015).
It is the ellipse of the rotated-axes formulas of Lefever (1926) and
Ebdon (1988) once the same divisor is used; texts differ only in the
divisor (`k`) and in the multiplier (some multiply by `sqrt(2)`). The
curve of the standard deviation by direction that Lefever drew is not
this ellipse (Gong 2002); the ellipse is returned. For bivariate normal
points the ellipse with `c = sqrt(-2 * log(1 - p))` covers a share `p`
of them: `c = 1` covers 39.3 per cent, `c = 2` 86.5 per cent.

## Usage

``` r
sdell(x, y = NULL, w = NULL, k = 1, c = 1, lonlat = FALSE)
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

- w:

  `NULL` for equal weights, or one non-negative weight per point (counts
  at a location, say), with a positive sum.

- k:

  The divisor offset: 0, 1 or 2; the total weight must exceed it.

- c:

  The multiplier of the standard deviations: one positive number.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A named numeric vector: `mx`, `my` (the mean centre, in degrees when
`lonlat = TRUE`), `a` and `b` (the semi-major and semi-minor axes,
`a >= b`, in the unit of the coordinates or metres), `angle` (of the
major axis, in radians anticlockwise from the x axis (east), in
`[0, pi)`; 0 when `a == b`) and `area` (`pi * a * b`).

## References

Wang, B., Shi, W. and Miao, Z. (2015). Confidence analysis of standard
deviational ellipse and its extension into higher dimensional Euclidean
space. *PLoS ONE*, 10(3), e0118537.
[doi:10.1371/journal.pone.0118537](https://doi.org/10.1371/journal.pone.0118537)

Gong, J. (2002). Clarifying the standard deviational ellipse.
*Geographical Analysis*, 34(2), 155-167.
[doi:10.1111/j.1538-4632.2002.tb01082.x](https://doi.org/10.1111/j.1538-4632.2002.tb01082.x)

## Examples

``` r
set.seed(2)
x <- rnorm(200)
sdell(x, 0.5 * x + rnorm(200, sd = 0.3), c = 2)
#>            mx            my             a             b         angle 
#> -0.0007417526  0.0681841044  2.4146550213  0.5262845280  0.4870207069 
#>          area 
#>  3.9923220528 
```
