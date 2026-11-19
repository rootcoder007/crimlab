# Days on which an entitlement was missed

The law may give a person held in a restrictive unit something every
day: hours out of the cell, time with other people. `entcmp()` takes,
for each stay, the percentage of its days on which that was missed, and
gives the stays by how often it was missed, the share that missed it
every day, and the share of all days missed. With the hours received on
the days it was missed, it shows how far short those days fell; with the
percentage of days the person refused, how many stays missed more days
than refusal accounts for. These are the tables of Sprott and Doob
(2021, Tables 6, 7, 14, 17 and 21) on the four hours out of cell and the
two hours of human contact in Canada's Structured Intervention Units.

## Usage

``` r
entcmp(
  miss,
  days = NULL,
  by = NULL,
  hours = NULL,
  ref = NULL,
  cuts = c(20, 50),
  level = 0.95
)
```

## Arguments

- miss:

  The percentage of the days of each stay on which the entitlement was
  missed: numbers from 0 to 100.

- days:

  `NULL`, or the length of each stay in days (numbers from 0), for the
  share of all days missed.

- by:

  `NULL`, or the group of each stay: labels.

- hours:

  `NULL`, or the hours received on average on the days the entitlement
  was missed (numbers from 0; `NA` where none was missed).

- ref:

  `NULL`, or the percentage of the days of each stay on which the person
  refused: numbers from 0 to 100.

- cuts:

  The percentages that part the bands between never and always:
  increasing numbers above 0 and below 100.

- level:

  The level of the exact intervals: one number between 0 and 1.

## Value

A list: `tab`, a data frame with one row per group (then `"all"`) and
band (`group`, `band`, `n`, `share`); `sum`, one row per group (`group`,
`n`, `none` the stays that never missed, `all` those that missed every
day, `share` their share, `lo`, `hi`, `half` the stays that missed on
half their days or more, `dayshr` the share of all days missed, `NA`
without `days`); `hrs` (with `hours`), the stays that missed by the
hours they received (`band`, `n`, `share`); `ref` (with `ref`), one row
per group (`group`, `n`, `most` the stays that refused on more than half
their days, `net` those that missed on more days than they refused); and
`level`.

## References

Sprott, J. B. and Doob, A. N. (2021). *Solitary Confinement, Torture,
and Canada's Structured Intervention Units*. University of Toronto,
Centre for Criminology and Sociolegal Studies, 23 February 2021.

## Examples

``` r
miss <- c(100, 100, 80, 45, 10, 0, 100, 60)
entcmp(miss, days = c(20, 35, 16, 40, 18, 5, 62, 30),
       hours = c(0.4, 1.5, 2.2, 3, 3.5, NA, 0.8, 1),
       ref = c(10, 100, 5, 0, 0, 0, 30, 60))
#> $tab
#>   group                 band n share
#> 1   all                    0 1 0.125
#> 2   all         over 0 to 20 1 0.125
#> 3   all        over 20 to 50 1 0.125
#> 4   all over 50 to under 100 2 0.250
#> 5   all                  100 3 0.375
#> 
#> $sum
#>   group n none all share         lo        hi half    dayshr
#> 1   all 8    1   3 0.375 0.08523341 0.7551368    5 0.7415929
#> 
#> $hrs
#>            band n     share
#> 1     up to 0.5 1 0.1428571
#> 2 over 0.5 to 1 2 0.2857143
#> 3   over 1 to 2 1 0.1428571
#> 4   over 2 to 3 2 0.2857143
#> 5        over 3 1 0.1428571
#> 
#> $ref
#>   group n most net
#> 1   all 8    2   5
#> 
#> $level
#> [1] 0.95
#> 
```
