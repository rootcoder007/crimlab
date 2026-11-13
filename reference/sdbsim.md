# Hot-spot dynamics: the Short-D'Orsogna-Brantingham model

Simulates the model of residential burglary of Short et al. (2008) on an
`n` by `n` grid of houses that wraps at its edges. Each site has an
attractiveness `a0 + b`; `b` rises by `theta` with every burglary there,
spreads to the four neighbours at rate `eta` and decays at rate `omega`,
so that burglaries draw further burglaries (repeat and near-repeat
victimisation). In the `"agent"` model offenders (`rho`, whole numbers
per site) burgle with probability \\1 - e^{-(a_0 + b) dt}\\ and leave,
or else move to a neighbouring site chosen in proportion to its
attractiveness; new offenders arrive at `gam` per site per unit time. In
the `"pde"` model the densities follow the continuum limit, \\b_t =
(\eta D / z) \Delta b - \omega b + \epsilon D \rho A\\ and \\\rho_t = (D
/ z)(\Delta \rho - 2 \nabla \cdot (\rho \nabla A / A)) - \rho A +
\gamma\\ with \\A = a_0 + b\\, by explicit Euler steps on a grid of
spacing `dx`; `sup` marks sites where policing suppresses crime (no
burglaries there). The step must satisfy
`dt <= z dx^2 / (4 dif max(1, eta))`.

## Usage

``` r
sdbsim(
  n,
  steps,
  mode = "agent",
  a0 = 1/30,
  omega = 1/15,
  eta = 0.03,
  gam = 0.002,
  theta = 1,
  dt = NULL,
  b = NULL,
  rho = NULL,
  dx = 1,
  dif = 1,
  eps = 1,
  z = 4,
  sup = NULL,
  seed = NULL
)
```

## Arguments

- n:

  The side of the grid: one whole number from 3 (at most 2000).

- steps:

  The number of steps: one whole number from 0.

- mode:

  `"agent"` or `"pde"`.

- a0:

  The static attractiveness: one positive number.

- omega:

  The decay rate of `b`: one number from 0 (`omega dt` at most 1).

- eta:

  The spread of `b` to neighbours: one number from 0 to 1.

- gam:

  The arrival rate of offenders per site: one number from 0.

- theta:

  The rise of `b` per burglary (`"agent"`): one number from 0.

- dt:

  The time step: one positive number; `NULL` for 1 (`"agent"`) or 0.01
  (`"pde"`).

- b, rho:

  `NULL` or `n` by `n` matrices of numbers from 0: the starting
  attractiveness (0 by default) and offenders (`"agent"`: whole numbers,
  none by default; `"pde"`: a density, `gam / a0` by default).

- dx, dif, eps, z:

  The grid spacing, the offenders' diffusion, the rise of `b` per unit
  of crime and the lattice's neighbour count of the `"pde"` model:
  positive numbers (`eps` from 0).

- sup:

  `NULL` or an `n` by `n` logical matrix of suppressed sites (`"pde"`).

- seed:

  `NULL` to draw the seed from R's generator, or one whole number from 0
  (`"agent"`).

## Value

An object of class `sdbsim`, a list: `b` and `rho`, the fields at the
end; `tot`, the burglaries at each site over the run (the integrated
crime rate for `"pde"`); `burg` and `meanb`, the burglaries and the mean
of `b` at each step; and `mode`.

## Details

Started at the uniform steady state
[`sdbthr()`](https://rootcoder007.github.io/crimlab/reference/sdbthr.md)
gives, the `"pde"` model stays there; started near it, a mode grows or
dies at the rate
[`sdbthr()`](https://rootcoder007.github.io/crimlab/reference/sdbthr.md)
predicts for its wavenumber. Without offenders the attractiveness only
decays: each step multiplies its largest value by at most `1 - omega dt`
(`proof("sdbsim")`) and its mean by exactly that.

## References

Short, M. B., D'Orsogna, M. R., Pasour, V. B., Tita, G. E., Brantingham,
P. J., Bertozzi, A. L. and Chayes, L. B. (2008). A statistical model of
criminal behavior. *Mathematical Models and Methods in Applied
Sciences*, 18(supp01), 1249-1267.
[doi:10.1142/S0218202508003029](https://doi.org/10.1142/S0218202508003029)

## Examples

``` r
s <- sdbsim(20, 50, a0 = 1 / 30, omega = 1 / 15, eta = 0.03,
            gam = 0.02, theta = 0.5, seed = 1)
s
#> sdbsim: agent model, 20 x 20 grid, 50 steps; burglaries 392; mean attractiveness at the end 0.1493598
p <- sdbsim(16, 200, mode = "pde", eps = 2, dt = 0.05,
            rho = matrix(0.07, 16, 16))
range(p$b)
#> [1] 0.05691214 0.05691214
```
