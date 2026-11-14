# Hit score of a geographic profile

A geographic profile ranks the cells of a search area; searched from the
highest score down, the hit score is the share of the area searched
before the offender's home is reached (Rossmo 2000), the smaller the
better. Cells scored the same as the home cell are counted both ways:
`hit` searches them after the home, `hitt` before it.

## Usage

``` r
hitscr(s, home)
```

## Arguments

- s:

  The scores of the cells: finite numbers, higher searched first (a
  profile from
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md),
  say).

- home:

  The cells holding the home: cell numbers from 1 to `length(s)`, one or
  several (several offenders, or candidate homes).

## Value

A data frame with one row per home: `home`, `above` (the cells scored
higher), `ties` (the cells scored the same, the home included), `rank`
(`above + 1`), and the hit scores in per cent, `hit`
(`100 * above / length(s)`) and `hitt`
(`100 * (above + ties) / length(s)`).

## References

Rossmo, D. K. (2000). *Geographic Profiling*. CRC Press.

## Examples

``` r
hitscr(c(0.9, 0.5, 0.7, 0.5, 0.1), home = 2)
#>   home above ties rank hit hitt
#> 1    2     2    2    3  40   80
```
