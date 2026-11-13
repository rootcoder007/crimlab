# Actuarial life table of grouped times

How long cases take, from the counts of cases ending (and censored) in
each band of time: the conditional probability of ending in each band
with the censored counted as half a case at risk, the survival (the
share still open at the end of each band) with Greenwood's standard
error and a plain, log or log-log interval, the hazard, and the median,
linear within its band.

## Usage

``` r
lftab(breaks, events, cens = 0, level = 0.95, ci = "log-log")
```

## Arguments

- breaks:

  The band edges, increasing (the last may be `Inf`), or band labels as
  [`bandx()`](https://rootcoder007.github.io/crimlab/reference/bandx.md)
  reads them (`"0-3"`, `"Over 12"`), contiguous.

- events:

  The cases ending in each band: numbers from 0.

- cens:

  The cases censored in each band: one number or one per band.

- level:

  The level of the intervals: one number between 0 and 1.

- ci:

  `"log-log"`, `"log"` or `"plain"`.

## Value

An object of class `lftab`, a list: `tab` (`lower`, `upper`, `entering`,
`atrisk`, `events`, `cens`, `q`, `surv`, `se`, `lo`, `hi`, `hazard`),
`median` (`Inf` when survival stays above one half or crosses it in an
open band) and `ci`.

## Examples

``` r
lftab(c(0, 3, 6, 12, 24, Inf), events = c(40, 30, 20, 10, 5),
      cens = c(2, 2, 3, 3, 4))
#> lftab: life table, log-log intervals, median 4.851282 
#>  lower upper entering atrisk events cens         q       surv         se
#>      0     3      119  118.0     40    2 0.3389831 0.66101695 0.04357669
#>      3     6       77   76.0     30    2 0.3947368 0.40008921 0.04548923
#>      6    12       45   43.5     20    3 0.4597701 0.21614015 0.03896034
#>     12    24       22   20.5     10    3 0.4878049 0.11070593 0.03110609
#>     24   Inf        9    7.0      5    4 0.7142857 0.03163027 0.02088777
#>           lo         hi     hazard
#>  0.568007817 0.73860814 0.13605442
#>  0.310878228 0.48760071 0.16393443
#>  0.145266720 0.29631969 0.09950249
#>  0.059213092 0.18020436 0.05376344
#>  0.006579202 0.09308838         NA
```
