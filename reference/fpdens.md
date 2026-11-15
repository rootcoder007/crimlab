# Fokker-Planck evolution of a count density

The density of the Ornstein-Uhlenbeck process fitted as in
[`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md),
evolved by the Fokker-Planck equation \\\partial_t p = -\partial_x
(\theta (\mu - x) p) + \frac{1}{2} \sigma^2 \partial_x^2 p\\, here with
the drift taken outside the derivative as in the explicit scheme used:
central differences on `ngrid` points from 0 to 1.5 times the largest
count plus 1, `nstep` steps of 0.05, edges copied from their neighbours,
the density kept non-negative and of total 1. It starts as a normal
curve at the last count with standard deviation \\\sigma + 1\\. The
stationary variance is \\\sigma^2 / (2 \theta)\\.

## Usage

``` r
fpdens(y, ngrid = 64, nstep = 200)
```

## Arguments

- y:

  The counts, in time order, three at least.

- ngrid:

  The grid points.

- nstep:

  The steps.

## Value

A list of class `fpdens`: `theta`, `mu`, `sigma`, `grid`, `density`,
`peak` (the grid point of the highest density) and `statvar`.

## References

Risken, H. (1989). *The Fokker-Planck Equation*, 2nd ed. Springer.

## See also

Other counts over time:
[`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md),
[`catcor()`](https://rootcoder007.github.io/crimlab/reference/catcor.md),
[`lvfit()`](https://rootcoder007.github.io/crimlab/reference/lvfit.md),
[`ousim()`](https://rootcoder007.github.io/crimlab/reference/ousim.md),
[`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)

## Examples

``` r
set.seed(1)
f <- fpdens(rpois(200, 20))
f$peak
#> [1] 19.71429
```
