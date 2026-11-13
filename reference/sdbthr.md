# Hot-spot threshold of the Short-D'Orsogna-Brantingham model

The continuum model of
[`sdbsim()`](https://rootcoder007.github.io/crimlab/reference/sdbsim.md)
has one uniform steady state, \\\bar b = \epsilon D \gamma / \omega\\,
\\\bar A = a_0 + \bar b\\, \\\bar\rho = \gamma / \bar A\\. Linearised
about it, a disturbance of wavenumber-squared \\k\\ evolves by a 2 by 2
matrix whose trace is negative at every \\k\\ and whose determinant is
\\P k^2 - Q k + R\\, \\R = \omega \bar A \> 0\\: a mode grows exactly
when the determinant is negative. So hot spots form when \\Q \> 0\\,
that is \\3 \epsilon D \bar\rho \> \omega + \eta \bar A\\ (repeat
victimisation outweighs decay and spreading), and \\Q^2 \> 4 P R\\; the
fastest-growing mode is near \\k^\* = Q / (2 P)\\. Without reinforcement
(`eps = 0`) the uniform state is stable at every wavelength
(`proof("sdbthr")`). `z` places the lattice's neighbour count on both
diffusion terms, as
[`sdbsim()`](https://rootcoder007.github.io/crimlab/reference/sdbsim.md)
does (`z = 1` is the continuum form).

## Usage

``` r
sdbthr(a0, omega, eps, dif, gam, eta, z = 1, k = NULL)
```

## Arguments

- a0, omega, eps, dif, gam, eta, z:

  The model's constants, as for
  [`sdbsim()`](https://rootcoder007.github.io/crimlab/reference/sdbsim.md):
  one number each, `eps` and `eta` from 0, the others positive.

- k:

  `NULL` or wavenumbers-squared (numbers from 0) at which to tabulate
  the trace, the determinant and the growth rate.

## Value

An object of class `sdbthr`, a list: `steady` (`b`, `a`, `rho`); `p`,
`q`, `r`, the determinant's coefficients; `thr`, the two sides of the
threshold (`lhs` \\3 \epsilon D \bar\rho\\, `rhs` \\\omega + \eta \bar
A\\); `unstable`; `kcrit`, `detcrit`, `grcrit`, the critical mode, its
determinant and growth rate (`NA` unless `p` and `q` are positive); and
`disp`, `NULL` or a data frame (`k`, `trace`, `det`, `growth`).

## References

Short, M. B., Bertozzi, A. L. and Brantingham, P. J. (2010). Nonlinear
patterns in urban crime: hotspots, bifurcations, and suppression. *SIAM
Journal on Applied Dynamical Systems*, 9(2), 462-483.
[doi:10.1137/090759069](https://doi.org/10.1137/090759069)

## Examples

``` r
sdbthr(a0 = 1, omega = 1, eps = 2, dif = 1, gam = 1, eta = 0.01)
#> sdbthr: steady b 2, rho 0.3333333; 3 eps dif rho = 2 against omega + eta A = 1.03
#> hot spots form; critical k 48.5
sdbthr(a0 = 1, omega = 1, eps = 0, dif = 1, gam = 1, eta = 0.01)$unstable
#> [1] FALSE
```
