# Search-circle hot-spot areas

Hot-spot areas from search circles (Block and Block 1995): circles of
radius `r` are centred on the nodes of an `nx` by `ny` grid spanning the
rectangle `win` (corners included), the points in each are counted, the
circles holding at least two points are ranked by count and the first
`top` are kept; kept circles that share a point are merged, repeatedly,
and a hot-spot area is the set of points in a merged group, kept if it
has at least `minpts`. Descriptions differ on the circle's radius and on
the grid; both are the caller's choice here.

## Usage

``` r
stacel(
  x,
  y = NULL,
  win,
  nx = 20,
  ny = 20,
  r,
  top = 25,
  minpts = 5,
  k = 1,
  c = 1
)
```

## Arguments

- x, y:

  The points in planar coordinates (project longitude and latitude with
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  first): `x` and `y` numeric vectors of equal length, or `x` a matrix
  or data frame with columns named `x` and `y` (an unnamed two-column
  matrix is read as `x` then `y`). At least two points, all inside
  `win`; no missing values.

- win:

  The rectangle the grid spans, `c(xmin, xmax, ymin, ymax)`; points may
  lie outside it.

- nx, ny:

  The numbers of grid nodes across and up: whole numbers from 2 to
  10000.

- r:

  The search radius, positive, in the unit of the coordinates.

- top:

  The number of busiest circles kept: one whole number from 1.

- minpts:

  The fewest points a hot-spot area keeps: a whole number from 2.

- k:

  The divisor offset: 0, 1 or 2; the total weight must exceed it.

- c:

  The multiplier of each cluster's standard deviational ellipse (see
  [`sdell()`](https://rootcoder007.github.io/crimlab/reference/sdell.md)).

## Value

A list. `cl`, the hot-spot area of each point (a row of `clusters`) or
`NA`, in input order. `clusters`, a data frame with one row per area, by
size: `points`, `mx`, `my` (mean centre), `wx`, `wy` (centre of minimum
distance), `a`, `b` and `angle` (the standard deviational ellipse of its
points).

## References

Block, R. L. and Block, C. R. (1995). Space, place and crime: hot spot
areas and hot places of liquor-related crime. In J. E. Eck and D.
Weisburd (eds.), *Crime and Place*, Crime Prevention Studies 4, 145-183.
Criminal Justice Press.

## Examples

``` r
set.seed(1)
x <- c(rnorm(30, 2, 0.2), runif(60, 0, 10))
y <- c(rnorm(30, 2, 0.2), runif(60, 0, 10))
stacel(x, y, win = c(0, 10, 0, 10), r = 0.75, top = 5, minpts = 5)$clusters
#>   points       mx       my       wx       wy         a         b     angle
#> 1     35 2.027899 2.025648 2.055272 2.024903 0.3233517 0.2532551 0.1950933
```
