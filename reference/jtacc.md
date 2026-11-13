# Accuracy of a journey-to-crime surface on a solved series

Given a surface over the grid of
[`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md)
and the offender's known base: the surface's value in the base's cell,
the hit score (the share of the grid searched, in decreasing order of
the surface, up to and including the base's cell), in the best and the
worst order of ties (Rossmo 2000; Canter et al. 2000), and the error
distance from the peak cell's centre to the base.

## Usage

``` r
jtacc(prob, win, base)
```

## Arguments

- prob:

  The surface: a matrix of values from 0, `nx` rows by `ny` columns, as
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md)
  returns.

- win:

  The grid's rectangle, as given to
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md).

- base:

  The known base, `c(x, y)`, inside `win`.

## Value

A named numeric vector: `p` (the value in the base's cell), `hitbest`
and `hitworst` (shares from 0 to 1), `err` (the error distance), `peak`
and `cell` (the peak's and the base's cell, numbered down the columns).

## References

Rossmo, D. K. (2000). *Geographic Profiling*. CRC Press.

## Examples

``` r
s <- jtsurf(c(2, 3, 4), c(2, 4, 3), win = c(0, 6, 0, 6), nx = 30,
            ny = 30, type = "negexp", par = c(1, 1), mode = "lik")
jtacc(s$prob, win = c(0, 6, 0, 6), base = c(3.2, 2.9))
#>            p      hitbest     hitworst          err         peak         cell 
#>   0.00932350   0.01000000   0.01111111   0.22360680 467.00000000 437.00000000 
```
