# Incapacitation under desistance and replacement

Crimes prevented by a sentence of `s` periods served from career age
`t0`, when the offending rate falls with career age and a share `r` of
the crimes prevented is committed by someone else instead. The count is
at most `s` times the rate at the start of the sentence and at least `s`
times the rate at its end; the same sentence served one period later
prevents no more; and, net of replacement, the estimate from the rate at
the start is still an upper bound. The usual constant-rate estimate
therefore overstates what a sentence buys under both assumptions. See
`proof("incapr")`.

## Usage

``` r
incapr(lam, t0, s, r = 0)
```

## Arguments

- lam:

  Offending rates by career age, `lam[1]` at age 0: non-negative numbers
  that never increase.

- t0:

  The career age at which the sentence starts: a whole number from 0,
  with `t0 + s` at most `length(lam) - 1`.

- s:

  The sentence in periods: a whole number from 1.

- r:

  The replacement share: one number from 0 to 1.

## Value

A named numeric vector: `prev` (crimes prevented), `net` (net of
replacement), `upper` (`s * lam` at age `t0`), `lower` (`s * lam` at age
`t0 + s`) and `later` (crimes prevented by the same sentence from age
`t0 + 1`).

## References

Shinnar, S. and Shinnar, R. (1975). The effects of the criminal justice
system on the control of crime: a quantitative approach. *Law and
Society Review*, 9(4), 581-611.
[doi:10.2307/3053340](https://doi.org/10.2307/3053340)

## Examples

``` r
incapr(lam = c(12, 10, 8, 6, 5, 4, 3, 2, 2, 1), t0 = 2, s = 3, r = 0.25)
#>  prev   net upper lower later 
#> 19.00 14.25 24.00 12.00 15.00 
```
