# Ornstein-Uhlenbeck (Langevin) fit and simulation of a count series

Treats a daily count series as a noisy mean-reverting intensity, \\dX_t
= \theta (\mu - X_t) dt + \sigma dW_t\\: \\\theta\\, \\\mu\\ and
\\\sigma\\ from the least-squares regression of the first differences on
the lagged values (\\\theta = -b\\, at least \\10^{-6}\\, \\\mu =
-a/b\\, \\\sigma\\ the residual standard deviation), then `npath`
Euler-Maruyama paths of `days / dt` steps from the last value. Draws
come from `seed`, as [`set.seed()`](https://rdrr.io/r/base/Random.html)
takes it.

## Usage

``` r
ousim(y, npath = 100, days = 365, dt = 1, seed = 42)
```

## Arguments

- y:

  The counts, in time order, three at least.

- npath:

  The paths.

- days:

  The horizon, in the time unit of `y`.

- dt:

  The step.

- seed:

  The seed.

## Value

A list of class `ousim`: `theta`, `mu`, `sigma`, `halflife` (\\\log 2 /
\theta\\), `q` (the 5%, 50% and 95% points of the last step), `paths`
(paths by steps) and `n`.

## References

Uhlenbeck, G. E. and Ornstein, L. S. (1930). On the theory of the
Brownian motion. *Physical Review*, 36(5), 823-841.
[doi:10.1103/PhysRev.36.823](https://doi.org/10.1103/PhysRev.36.823)

Kloeden, P. E. and Platen, E. (1992). *Numerical Solution of Stochastic
Differential Equations*. Springer.
[doi:10.1007/978-3-662-12616-5](https://doi.org/10.1007/978-3-662-12616-5)

## See also

Other counts over time:
[`arimfc()`](https://rootcoder007.github.io/crimlab/reference/arimfc.md),
[`catcor()`](https://rootcoder007.github.io/crimlab/reference/catcor.md),
[`fpdens()`](https://rootcoder007.github.io/crimlab/reference/fpdens.md),
[`lvfit()`](https://rootcoder007.github.io/crimlab/reference/lvfit.md),
[`sarfc()`](https://rootcoder007.github.io/crimlab/reference/sarfc.md)

## Examples

``` r
set.seed(1)
y <- rpois(200, 20)
ousim(y, npath = 20, days = 30)$q
#> [1] 11.14984 19.07234 24.60097
```
