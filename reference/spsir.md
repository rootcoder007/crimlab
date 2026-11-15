# The metapopulation SIR model

The deterministic SIR model on patches (Keeling and Rohani 2008, sec.
7.3): the force of infection on patch `i` is `beta sum_j C_ij I_j / N_j`
with the coupling `C = (1 - kappa) I + kappa W~`, `W~` the
row-standardised neighbour weights, so that a share `kappa` of the
pressure comes from neighbours; recovery is at rate `gamma`. Integrated
by the classical fourth-order Runge-Kutta method in steps of `dt`. The
same equations describe the spread of an outbreak of violence between
neighbourhoods in contagion models.

## Usage

``` r
spsir(
  i0,
  w = NULL,
  pop = NULL,
  beta = 0.5,
  gamma = 0.2,
  kappa = 0.1,
  tmax = 100,
  dt = 0.1
)
```

## Arguments

- i0:

  The infected at time 0, one number from 0 per patch.

- w:

  `NULL` for isolated patches, or the neighbours, as in
  [`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md).

- pop:

  `NULL` for 1000 per patch, or the populations (positive, at least
  `i0`).

- beta, gamma:

  The transmission and recovery rates: numbers from 0.

- kappa:

  The share of the pressure from neighbours: from 0 to 1.

- tmax, dt:

  The horizon and the step: positive; at most 1e6 steps.

## Value

A list of class `spsir`: `attack` (the final attack rate,
`sum R / sum N`), `r0` (`beta / gamma`), `t` (the times) and `S`, `I`,
`R` (matrices: a row per time, a column per patch).

## References

Keeling, M. J. and Rohani, P. (2008). *Modeling Infectious Diseases in
Humans and Animals*. Princeton University Press, Princeton.
[doi:10.1515/9781400841035](https://doi.org/10.1515/9781400841035)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md)

## Examples

``` r
w <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
spsir(c(10, 0, 0), w, rep(1000, 3), beta = 0.6, gamma = 0.2,
      tmax = 60, dt = 0.5)
#> spsir: 3 patches to time 60, attack rate 0.9399, R0 3
```
