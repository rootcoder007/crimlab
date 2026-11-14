# Phase diagram of the inspection game

The inspection game of Helbing et al. (2010) with cooperators (C),
criminals (P) and inspectors (O), against C, P, O: C earns (1, 0, 1), P
earns (T, 0, T - 1) for the temptation `T`, O earns `1 - g` whatever it
meets for the inspection cost `g`. Discrete replicator steps
`x += rate * x * (f - mean f)` (negative shares set to 0, the shares
normalised) run from `x0` at each point of a grid of temptations and
costs; the final shares map the phases of the game, the criminals' share
standing for the crime rate.

## Usage

``` r
insphs(
  temp = seq(0.05, 1.8, length.out = 20),
  cost = seq(0.05, 1.2, length.out = 20),
  steps = 600,
  rate = 0.05,
  x0 = c(0.34, 0.33, 0.33)
)
```

## Arguments

- temp, cost:

  The temptations and inspection costs of the grid: finite numbers.

- steps:

  The replicator steps at each point: one whole number from 0.

- rate:

  The step size: one positive number.

- x0:

  The starting shares of C, P and O: three numbers from 0 with a
  positive sum.

## Value

A list: `crime`, `coop` and `insp`, the final shares of P, C and O
(matrices, one row per temptation and one column per cost); `temp` and
`cost`.

## References

Helbing, D., Szolnoki, A., Perc, M. and Szabo, G. (2010). Punish, but
not too hard: how costly punishment spreads in the spatial public goods
game. *New Journal of Physics*, 12, 083005.
[doi:10.1088/1367-2630/12/8/083005](https://doi.org/10.1088/1367-2630/12/8/083005)

## Examples

``` r
ph <- insphs(seq(0.1, 1.8, by = 0.3), seq(0.1, 1.1, by = 0.25),
             steps = 200)
round(ph$crime, 2)
#>      0.1 0.35  0.6 0.85  1.1
#> 0.1 0.00 0.00 0.00 0.00 0.00
#> 0.4 0.00 0.00 0.00 0.00 0.00
#> 0.7 0.00 0.01 0.02 0.03 0.04
#> 1   0.00 0.06 0.19 0.29 0.35
#> 1.3 0.03 0.25 0.52 0.64 0.69
#> 1.6 0.11 0.40 0.67 0.79 0.82
```
