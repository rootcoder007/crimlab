# Centre of minimum distance

The point `c` that minimises the weighted sum of distances
`sum(w * |p - c|)` to the points: the Weber point or spatial median,
which, unlike the median centre, does not depend on the axes. It is
found by Weiszfeld's iteration with the modification of Vardi and Zhang
(2000), which keeps every step defined and the sum non-increasing when
an iterate lands on a data point. The optimality gap certifies the
answer: it is 0 exactly at a minimiser. The minimiser is unique unless
all the points lie on one line with the weight split evenly.

## Usage

``` r
webpt(
  x,
  y = NULL,
  w = NULL,
  start = NULL,
  tol = 1e-10,
  maxit = 1000,
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

- w:

  `NULL` for equal weights, or one non-negative weight per point (counts
  at a location, say), with a positive sum.

- start:

  `NULL` to start at the mean centre, or `c(x, y)` (in degrees when
  `lonlat = TRUE`).

- tol:

  Stop when a step is at most `tol` times the root mean squared distance
  to the mean centre: one positive number.

- maxit:

  The largest number of steps: one whole number from 0.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A named numeric vector: `x`, `y` (the centre, in degrees when
`lonlat = TRUE`), `obj` (the weighted sum of distances), `gap` (the
optimality gap, 0 at a minimiser), `iter` (steps taken) and `code`: 2 at
an exact minimiser, 1 when the step fell below `tol`, 0 when `maxit` was
reached.

## References

Weiszfeld, E. and Plastria, F. (2009). On the point for which the sum of
the distances to n given points is minimum. *Annals of Operations
Research*, 167, 7-41 (translation of the 1937 paper).
[doi:10.1007/s10479-008-0352-z](https://doi.org/10.1007/s10479-008-0352-z)

Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and
associated data depth. *Proceedings of the National Academy of
Sciences*, 97(4), 1423-1426.
[doi:10.1073/pnas.97.4.1423](https://doi.org/10.1073/pnas.97.4.1423)

## Examples

``` r
webpt(c(0, 1, 0, 10), c(0, 0, 1, 10))
#>            x            y          obj          gap         iter         code 
#> 5.000000e-01 5.000000e-01 1.555635e+01 4.086973e-10 2.200000e+01 1.000000e+00 
```
