# A Gaussian-plume concentration field on a grid

Evaluates
[`gplume`](https://rootcoder007.github.io/crimlab/reference/pgsig.md) at
every point of an `x` by `y` ground-level grid (or at height `z`) and
returns the field as a matrix with the axes. It maps where the plume
from a single stack touches the ground, how far downwind the highest
concentration lies, and how wide the affected area is.

## Usage

``` r
plmfld(q, u, h, x, y, z = 0, stab = "D", set = "rural", mixh = NULL, nimg = 3)

# S3 method for class 'plmfld'
print(x, ...)
```

## Arguments

- q, u, h, stab, set, mixh, nimg:

  As in
  [`gplume`](https://rootcoder007.github.io/crimlab/reference/pgsig.md).

- x, y:

  Numeric vectors, strictly increasing: the downwind and crosswind
  coordinates of the grid, in metres from the source.

- z:

  Numeric scalar, the receptor height in metres; default `0` (ground
  level).

- ...:

  Not used.

## Value

A list of class `plmfld`: `field` (a numeric matrix of concentrations,
rows along `x`, columns along `y`), `x`, `y`, `dx` and `dy` (median
spacings, `NA` for a single coordinate), `units` (`"concentration"`) and
`title`.

## Details

Every grid point \\(x_i, y_j, z)\\ is passed to
[`gplume`](https://rootcoder007.github.io/crimlab/reference/pgsig.md) in
one call, with the source at the origin and the wind along \\+x\\; the
field is therefore zero (or vanishingly small) upwind and at \\x \le
0\\. The concentration is in the mass unit of `q` per cubic metre (grams
per cubic metre for `q` in grams per second). `dx` and `dy` are the
median grid spacings. The cost is one plume evaluation per grid point.

Errors: `x` or `y` not numeric or not strictly increasing; `z` negative;
and every check of
[`gplume`](https://rootcoder007.github.io/crimlab/reference/pgsig.md).

## See also

[`gplume`](https://rootcoder007.github.io/crimlab/reference/pgsig.md).

## Examples

``` r
f <- plmfld(100, 5, 50, x = seq(100, 3000, by = 100),
            y = seq(-500, 500, by = 50))
dim(f$field)
#> [1] 30 21
# where the ground-level maximum lies, downwind
f$x[which.max(apply(f$field, 1, max))]
#> [1] 800
# the coordinates must increase
try(plmfld(100, 5, 50, x = c(300, 100), y = c(-50, 0, 50)))
#> Error : plmfld: `x` and `y` must be strictly increasing receptor coordinates.
```
