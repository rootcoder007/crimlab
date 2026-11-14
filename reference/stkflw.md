# Stock and flow of a population held over time

Person-days divided by the length of a period give the average daily
population (a stock: how many are held at once); divided by the people
admitted or released, the average length of stay (how long each stays).
The two share a numerator and `days = adp * t = people * alos` (Little's
law in the words of corrections planning; Lakner 1976), so any two of
people, days, average daily population and average length of stay give
the rest, except days with the average daily population, which carry the
same information. When the length of stay moves, the stock and the flow
can move in opposite directions; `pshare` is the share of the change in
days (on the log scale) that the change in people accounts for, the rest
coming from the length of stay. See `proof("stkflw")`.

## Usage

``` r
stkflw(
  people = NULL,
  days = NULL,
  adp = NULL,
  alos = NULL,
  t = 365,
  prd = NULL,
  from = NULL,
  to = NULL,
  expo = NULL,
  per = 1e+05,
  base = "first",
  rel = NULL,
  early = NULL
)
```

## Arguments

- people, days, adp, alos:

  Two of: the people admitted or released in each period (positive), the
  person-days (from 0), the average daily population (from 0; or a list
  with the periodic headcounts of each period, whose means are taken),
  and the average length of stay (positive). A `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  is taken too (as `people`): its elements fill the arguments by name.

- t:

  The length of each period in days: one positive number or one per
  period.

- prd:

  `NULL` (1, 2, ...) or the label of each period.

- from, to:

  `NULL`, or the first and last day of each period (dates or date
  strings); then `t` is the days from one to the other, counting both.

- expo:

  `NULL`, or the population at risk: one positive number or one per
  period, for rates per `per`.

- per:

  The rates are per `per` people: one positive number.

- base:

  Changes against the `"first"` period or the `"previous"` one.

- rel, early:

  `NULL`, or the releases of each period and the days by which every
  release is brought forward.

## Value

An object of class `stkflw`, a list: `flow`, a data frame with one row
per period (`period`, `t`, `people`, `days`, `adp`, `alos`, their
percent changes `peochg`, `aloschg`, `dayschg`, `adpchg`, `pshare`, with
`expo` `flowrate`, `stkrate`, `frchg`, `srchg`, and with `rel` `cut`,
`cutpct`); `mean`, the means over the periods of `people`, `days`, `adp`
and `alos`; `base` and `per`.

## Details

Given the releases of each period and a number of days by which every
release is brought forward, `cut` is the stock that removes in a steady
state, `rel * early / t` (releasing everyone six months early removes
six months of releases at once), and `cutpct` that share of the average
daily population.

## References

Lakner, E. (1976). *A Manual of Statistical Sampling Methods for
Corrections Planners*. University of Illinois at Urbana-Champaign.

Little, J. D. C. (1961). A proof for the queuing formula: L = lambda W.
*Operations Research*, 9(3), 383-387.
[doi:10.1287/opre.9.3.383](https://doi.org/10.1287/opre.9.3.383)

## Examples

``` r
stkflw(people = c(5200, 4800, 4300), days = c(310000, 318000, 330000),
       prd = 2021:2023)
#> stkflw: stock and flow, changes against the first period
#>  period   t people   days      adp     alos     peochg  aloschg  dayschg
#>    2021 365   5200 310000 849.3151 59.61538   0.000000  0.00000 0.000000
#>    2022 365   4800 318000 871.2329 66.25000  -7.692308 11.12903 2.580645
#>    2023 365   4300 330000 904.1096 76.74419 -17.307692 28.73218 6.451613
#>    adpchg    pshare
#>  0.000000        NA
#>  2.580645 -3.141506
#>  6.451613 -3.039708
# admissions implied by a population and a length of stay
stkflw(adp = 900, alos = 61)$flow$people
#> [1] 5385.246
```
