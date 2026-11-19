# Thin-plate spline surface

`tpsurf()` fits a thin-plate spline to values measured at points in the
plane (rates at area centroids, say) and evaluates it on a grid over
their box: a smooth surface that passes through the values when `lam` is
0 and is smoothed towards a plane as `lam` grows. It complements
[`kdsurf()`](https://rootcoder007.github.io/crimlab/reference/kdsurf.md),
which smooths the points themselves rather than values at them.

## Usage

``` r
tpsurf(x, y, z, lam = 0, nx = 50L, ny = nx)

# S3 method for class 'tpsurf'
print(x, ...)

# S3 method for class 'tpsurf'
plot(x, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- x, y:

  The coordinates of the points. For
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), `x` is a
  `tpsurf` result.

- z:

  The value at each point.

- lam:

  The smoothing: one number from 0 (0 interpolates).

- nx, ny:

  The size of the grid: whole numbers from 2 to 1000.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html); for
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

- main:

  `NULL` (a default title) or one string.

- col:

  Not used: the surface is shaded by `pal`.

- pal:

  `NULL` (`"sequential"`) or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

## Value

An object of class `tpsurf`, a list: `x` and `y` (the grid), `z` (the
surface on it, `nx` by `ny`, as
[`graphics::image()`](https://rdrr.io/r/graphics/image.html) takes it),
`a` (the basis coefficients), `beta` (`const`, `x`, `y`), `fitted`,
`resid` and `lam`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
surface with its contours and returns `x`, `y` and `z` invisibly.

## Details

The surface is `f(s) = sum_i a_i phi(|s -` `s_i|) + b_0 + b_1 x + b_2 y`
with the radial basis `phi(r) = r^2 log r` (Duchon 1977); the
coefficients solve the linear system `(K + lam I) a + T b = z`,
`T'a = 0`, `K` the basis between the points and `T` the rows
`(1, x_i, y_i)` (Wahba 1990, ch. 2). The compiled kernel solves it by LU
decomposition and evaluates the grid. It refuses points that all lie on
one line (the plane is not determined) and, with `lam = 0`, points that
repeat (the system is singular). `fitted` is the surface at the points,
so `resid` is 0 for `lam = 0` up to rounding.

## References

Duchon, J. (1977). Splines minimizing rotation-invariant semi-norms in
Sobolev spaces. In W. Schempp and K. Zeller (eds.), *Constructive Theory
of Functions of Several Variables*, Lecture Notes in Mathematics 571,
85-100. Springer.
[doi:10.1007/BFb0086566](https://doi.org/10.1007/BFb0086566)

Wahba, G. (1990). *Spline Models for Observational Data*. SIAM.
[doi:10.1137/1.9781611970128](https://doi.org/10.1137/1.9781611970128)

## Examples

``` r
set.seed(2)
x <- runif(40)
y <- runif(40)
s <- tpsurf(x, y, sin(4 * x) + y^2 + rnorm(40, sd = 0.05), lam = 0.001,
            nx = 30)
s
#> tpsurf: thin-plate spline of 40 points, lambda 0.001, residual sum of squares 0.002046; grid 30 x 30
plot(s)
```
