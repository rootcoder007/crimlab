# The judge-leniency instrument

When cases are assigned to judges at random, the judge's leniency is an
instrument for the decision (detention, a custodial sentence): the
difference in outcomes between lenient and strict judges over the
difference in decisions is the Wald ratio (Imbens and Angrist 1994). It
is the average effect on the compliers, the cases whose decision depends
on the judge, only if no case defies (is detained by the lenient judge
and released by the strict one). With a share of defiers `s` whose
effect is `e`, the compliers number the first stage plus `s` and their
effect is `(wald * fs + s * e) / (fs + s)`; `sens` gives it over a grid,
so a reader sees how far the conclusion rests on monotonicity (de
Chaisemartin 2017). Code `z` so that the first stage is positive (1 for
the judge whose decision is more often 1); with a negative one the
complier share is negative and `ceff` is `NA`.

## Usage

``` r
jdgiv(z = NULL, d, y, w = NULL, dsh = c(0, 0.05, 0.1), deff = 0)
```

## Arguments

- z:

  `NULL`, or the instrument of each case: 0 for the strict judge, 1 for
  the lenient (or `FALSE`, `TRUE`).

- d:

  The decision of each case: 0 or 1; with `z = NULL`, two columns (a
  matrix or data frame), the decisions under the strict and under the
  lenient judge.

- y:

  The outcome of each case: numbers; with `z = NULL`, two columns, the
  outcomes without and with the decision.

- w:

  `NULL`, or a weight from 0 for each case.

- dsh:

  Defier shares for the sensitivity grid: shares from 0 to 1.

- deff:

  Defier effects for the sensitivity grid: numbers.

## Value

An object of class `jdgiv`, a list: `type` (`"obs"` or `"pop"`); `est`,
a named vector (observed: `fs` the first stage, `itt`, `wald`;
population: `itt`, `fs`, `wald`, `ate`); `shares` (observed: `complier`,
`always`, `never` if no one defies; population: also `defier`); and
`sens` (observed: `dsh`, `deff`, `cshare` the complier share, `ceff`
their effect) or `effects` (population: the mean effect of each type,
`NA` for an empty one).

## Details

Given the decisions and outcomes of a whole population under both judges
(`z = NULL`, `d` and `y` two columns each), the exact decomposition: the
shares and mean effects of compliers, defiers, always- and never-takers,
the intention-to-treat effect and the first stage, whose ratio is the
compliers' effect when there are no defiers. See `proof("jdgiv")`.

## References

Imbens, G. W. and Angrist, J. D. (1994). Identification and estimation
of local average treatment effects. *Econometrica*, 62(2), 467-475.
[doi:10.2307/2951620](https://doi.org/10.2307/2951620)

de Chaisemartin, C. (2017). Tolerating defiance? Local average treatment
effects without monotonicity. *Quantitative Economics*, 8(2), 367-396.
[doi:10.3982/QE601](https://doi.org/10.3982/QE601)

## Examples

``` r
set.seed(7)
z <- rbinom(1000, 1, 0.5)
d <- rbinom(1000, 1, ifelse(z == 1, 0.3, 0.5))
y <- rbinom(1000, 1, 0.2 + 0.1 * d)
jdgiv(z, d, y)
#> jdgiv: observed 
#>          fs         itt        wald 
#> -0.17071258  0.01969084 -0.11534496 
#>   complier     always      never 
#> -0.1707126  0.4723320  0.6983806 
#>   dsh deff      cshare ceff
#>  0.00    0 -0.17071258   NA
#>  0.05    0 -0.12071258   NA
#>  0.10    0 -0.07071258   NA
```
