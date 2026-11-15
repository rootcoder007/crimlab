# The BYM2 area structure

The scaled structure of the BYM2 model of Riebler et al. (2016): the
generalised inverse `G` of the Laplacian `Q` of the neighbour weights
(from its eigenvalues above `tol` times the largest), the scaling factor
`s`, the geometric mean of the diagonal of `G`, so that the structured
effect has a generalised variance of 1, and the covariance of the area
effects `((1 - phi) I + phi G / s) / tau`. The graph must be connected
and every area must have a neighbour.

## Usage

``` r
bym2sc(w, tau = 1, phi = 0.5, tol = 1e-10)
```

## Arguments

- w:

  The neighbours:
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)'s
  list, an spdep `listw` or an `n x n` matrix of weights from 0 (its
  non-zero entries off the diagonal).

- tau:

  The precision: one positive number.

- phi:

  The share of the variance that is structured: one number from 0 to 1.

- tol:

  The relative eigenvalue tolerance: one number strictly between 0 and
  1.

## Value

A list: `scale`, `ginv` and `cov`.

## References

Riebler, A., Sorbye, S. H., Simpson, D. and Rue, H. (2016). An intuitive
Bayesian spatial model for disease mapping that accounts for scaling.
*Statistical Methods in Medical Research*, 25(4), 1145-1165.
[doi:10.1177/0962280216660421](https://doi.org/10.1177/0962280216660421)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
w <- rbind(c(0, 1, 1, 0), c(1, 0, 1, 0), c(1, 1, 0, 1), c(0, 0, 1, 0))
bym2sc(w, 2, 0.4)$scale
#> [1] 0.3565926
```
