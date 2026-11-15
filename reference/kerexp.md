# Exposure by kernel

The kernel-weighted density of sources at each receptor,
`sum_j w_j K(d_ij / h) / h^2`, with the bivariate normal kernel
`exp(-u^2 / 2) / (2 pi)` or the quartic kernel `3 / pi (1 - u^2)^2` (for
`u < 1`, 0 beyond) of Bailey and Gatrell (1995). Planar coordinates.

## Usage

``` r
kerexp(x, y = NULL, src, h, w = NULL, kern = "gaussian")
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

- h:

  The bandwidth: one positive number.

- w:

  `NULL` for a weight of 1 per source, or one weight from 0 per source.

- kern:

  `"gaussian"` or `"quartic"`.

## Value

A numeric vector, one value per receptor.

## References

Bailey, T. C. and Gatrell, A. C. (1995). *Interactive Spatial Data
Analysis*. Longman, Harlow.

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
kerexp(rbind(c(0.2, 0.3)), src = rbind(c(0.25, 0.35), c(0.6, 0.8)),
       h = 0.5, kern = "quartic")
#> [1] 3.668458
```
