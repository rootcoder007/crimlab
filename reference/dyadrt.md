# Dyadic offending rates

Offences by one group against another, as a rate per member of the
offending group and as a rate per offender-victim pair exposed. The rate
per offender rises with the size of the victim group even when offenders
choose victims at random; the rate per pair does not. See
`proof("dyadrt")`.

## Usage

``` r
dyadrt(count, pop)
```

## Arguments

- count:

  A square matrix (or 2-d table) of offence counts, offending group in
  rows and victim group in columns: non-negative numbers. When it names
  both its rows and its columns, the columns are matched to the rows by
  name, and both to `pop` when `pop` is named.

- pop:

  Group sizes: positive numbers, as a numeric vector, a 1-d table, a
  list of single numbers or a one-row data frame. Unnamed, they are
  taken in the order of the rows: only names protect against a mismatch.

## Value

A data frame with one row per dyad: `off`, `vic`, `count`, `rtoff`,
`rtpair`, `nulloff` (the rate per offender under random mixing) and
`ratio` (`rtpair` over the pooled pair rate; 1 under random mixing).

## Examples

``` r
pop <- c(a = 800, b = 200)
dyadrt(count = outer(pop, pop) / 1000, pop = pop)
#>   off vic count rtoff rtpair nulloff ratio
#> 1   a   a   640   0.8      1     0.8     1
#> 2   b   a   160   0.8      1     0.8     1
#> 3   a   b   160   0.2      1     0.2     1
#> 4   b   b    40   0.2      1     0.2     1
```
