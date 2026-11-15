# The Leroux precision matrix

The precision of the area effects in the model of Leroux, Lei and
Breslow (2000), `tau (rho Q + (1 - rho) I)`, with `Q = D - W` the
Laplacian of the neighbour weights: independent effects at `rho = 0`,
the intrinsic conditional autoregression at `rho = 1`.

## Usage

``` r
lerprc(w, rho, tau = 1)
```

## Arguments

- w:

  The neighbours:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, an spdep `listw` or an `n x n` matrix of weights from 0 (its
  non-zero entries off the diagonal).

- rho:

  The spatial parameter: one number from 0 to 1.

- tau:

  The precision: one positive number.

## Value

The `n x n` precision matrix.

## References

Leroux, B. G., Lei, X. and Breslow, N. (2000). Estimation of disease
rates in small areas: a new mixed model for spatial dependence. In M. E.
Halloran and D. Berry (eds), *Statistical Models in Epidemiology, the
Environment, and Clinical Trials*, 179-191. Springer, New York.
[doi:10.1007/978-1-4612-1284-3_4](https://doi.org/10.1007/978-1-4612-1284-3_4)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
lerprc(rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0)), 0.5, 2)
#>      [,1] [,2] [,3]
#> [1,]    2   -1    0
#> [2,]   -1    3   -1
#> [3,]    0   -1    2
```
