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

  A square matrix of offence counts, offending group in rows and victim
  group in columns: non-negative numbers.

- pop:

  Group sizes in the order of the rows: positive numbers, optionally
  named.

## Value

A data frame with one row per dyad: `off`, `vic`, `count`, `rate_off`,
`rate_pair`, `null_off` (the rate per offender under random mixing) and
`ratio` (`rate_pair` over the pooled pair rate; 1 under random mixing).

## Examples

``` r
pop <- c(a = 800, b = 200)
dyadrt(count = outer(pop, pop) / 1000, pop = pop)
#>   off vic count rate_off rate_pair null_off ratio
#> 1   a   a   640      0.8         1      0.8     1
#> 2   b   a   160      0.8         1      0.8     1
#> 3   a   b   160      0.2         1      0.2     1
#> 4   b   b    40      0.2         1      0.2     1
```
