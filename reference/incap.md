# Incapacitation: crime rate, share prevented and the marginal year

An offender commits crimes at rate `lam` a year while free, and each
crime leads to a custodial sentence of `s` years with probability `q`.
Over a career the crime rate is `lam / (1 + lam q s)`: the share of
offending prevented, `lam q s / (1 + lam q s)`, never reaches 1, longer
or surer sentences lower the rate, and one more year of sentence
prevents `lam^2 q / ((1 + lam q s) (1 + lam q (s + 1)))` crimes a year,
less for each year added. A sentence uniform across offenders prevents a
larger share of the offending of those with the higher rate. The model
holds the rate fixed over the career and replaces nobody; see
[`incapr()`](https://rootcoder007.github.io/crimlab/reference/incapr.md)
for both. See `proof("incap")`.

## Usage

``` r
incap(lam, q, s, share = NULL)
```

## Arguments

- lam:

  Offending rates while free (crimes a year): positive numbers, one per
  group, as a numeric vector, a 1-d table, a list of single numbers or a
  one-row data frame.

- q:

  The chance a crime leads to a custodial sentence: numbers in the
  interval from 0 (excluded) to 1, recycled to the length of `lam`; one
  per group, in the forms of `lam` and matched to it by name when both
  are named (unnamed, by position: only names protect against a
  mismatch).

- s:

  Sentence lengths in years: non-negative numbers, recycled; one per
  group, in the forms of `lam`, matched by name.

- share:

  `NULL`, or one non-negative population share per group, summing to 1,
  in the forms of `lam`, matched by name.

## Value

A data frame with one row per group: `group`, `lam`, `q`, `s`, `rate`
(the long-run crime rate), `prev` (the share prevented) and `marg` (the
crimes a year one more year of sentence prevents). With `share`, a last
row `all` holds the weighted free rate in `lam`, the weighted rate in
`rate` and the share of all offending prevented in `prev`.

## References

Avi-Itzhak, B. and Shinnar, R. (1973). Quantitative models in crime
control. *Journal of Criminal Justice*, 1(3), 185-217.
[doi:10.1016/0047-2352(73)90061-5](https://doi.org/10.1016/0047-2352%2873%2990061-5)

## Examples

``` r
incap(lam = c(low = 2, high = 10), q = 0.1, s = 1, share = c(0.8, 0.2))
#>   group  lam   q  s     rate      prev      marg
#> 1   low  2.0 0.1  1 1.666667 0.1666667 0.2380952
#> 2  high 10.0 0.1  1 5.000000 0.5000000 1.6666667
#> 3   all  3.6  NA NA 2.333333 0.3518519        NA
```
