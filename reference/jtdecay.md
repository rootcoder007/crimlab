# Distance-decay functions for journey to crime

The decay of the chance of offending with distance from an offender's
base, in the forms of Brantingham and Brantingham (1981), Rossmo (2000)
and Levine (2009): `"linear"`, `f = max(0, A + B d)`, `par = c(A, B)`;
`"negexp"`, `f = A exp(-b d)`, `par = c(A, b)` with `b >= 0`;
`"normal"`, `f = A dnorm(d, mu, s)`, `par = c(A, mu, s)`; `"lognormal"`,
`f = A dlnorm(d, mu, s)`, `par = c(A, mu, s)` (the lognormal density of
`d`; some texts use `log d^2` and `1 / d^2`, which is not a density);
`"trunc"`, the truncated negative exponential, `f = P d / dp` up to the
peak `dp` and `P exp(-C (d - dp))` beyond, `par = c(P, dp, C)`; and
`"table"`, values `tab` at distances `0, delta, 2 delta, ...`, linear in
between and 0 beyond, `par = delta` (as
[`jtcal()`](https://rootcoder007.github.io/crimlab/reference/jtcal.md)
tabulates).

## Usage

``` r
jtdecay(d, type, par, tab = NULL)
```

## Arguments

- d:

  Distances: numbers from 0.

- type:

  The form: one of the names above.

- par:

  The parameters of the form, of the length given above.

- tab:

  For `type = "table"`: at least two values from 0.

## Value

A numeric vector of `f(d)`, one per distance.

## References

Brantingham, P. J. and Brantingham, P. L. (eds.) (1981). *Environmental
Criminology*. Sage.

Levine, N. (2009). Introduction to the special issue on Bayesian
journey-to-crime modelling. *Journal of Investigative Psychology and
Offender Profiling*, 6(3), 167-185.
[doi:10.1002/jip.107](https://doi.org/10.1002/jip.107)

## Examples

``` r
jtdecay(c(0, 0.5, 1, 2, 5), type = "trunc", par = c(1, 0.5, 0.8))
#> [1] 0.00000000 1.00000000 0.67032005 0.30119421 0.02732372
```
