# Repeat victimisation: how often the same targets are victimised again

Asks how much crime falls on targets (people, households, addresses,
businesses) already victimised, and how soon a repeat follows. From one
record per incident with the target's identifier (and its date):

## Usage

``` r
rptvic(id, t = NULL, pop = NULL, brk = NULL)
```

## Arguments

- id:

  The target of each incident: labels.

- t:

  `NULL` or the time of each incident: numbers (days, say).

- pop:

  `NULL` or the number of targets at risk: one number, at least the
  number of victims.

- brk:

  `NULL` or the edges of the intervals of time since the previous
  incident: increasing numbers, two at least (needs `t`).

## Value

An object of class `rptvic`, a list: `tab`, one row per number of
victimisations `k` (from 0 with `pop`): `victims`, `incidents` and, with
`pop`, `expected`, the Poisson's victims; `stat`, a named vector:
`incidents`, `victims`, `repeats`, `rshare`, `conc`, `incid` and `prev`
(`NA` without `pop`), `rvict` (victims with two or more), `erep` (the
Poisson's number of them) and `ratio`; `gap`, `NULL` without `t`,
otherwise `n`, `mean` and `median` of the times between successive
incidents; and `lag`, `NULL` without `brk`, otherwise one row per
interval (`lo`, `hi`, `n`, `share`, `cum`).

## Details

- incidence, prevalence and concentration (Trickett, Osborn, Seymour and
  Pease 1992): incidents per target at risk, victims per target at risk,
  and incidents per victim; and the repeats, the incidents beyond each
  victim's first, with their share of all incidents;

- the victims by number of victimisations, against the numbers a Poisson
  count with the same mean over the `pop` targets at risk would give:
  more victims with two or more than the Poisson expects (`ratio`
  above 1) is victimisation concentrated beyond chance, the classic
  finding for burglary;

- the time between a target's successive incidents: their number, mean
  and median, and with `brk` the share of repeats falling in each
  interval since the previous incident, which for burglary is highest
  soon after (Johnson 2008).

A count of victims misses the targets never victimised: `pop` supplies
them. The Poisson comparison assumes every target has the same risk, so
it measures how far risk varies, and does not by itself tell a target's
lasting attractiveness ("flag") from a first incident raising the risk
of the next ("boost").

## References

Trickett, A., Osborn, D. R., Seymour, J. and Pease, K. (1992). What is
different about high crime areas? *The British Journal of Criminology*,
32(1), 81-89.
[doi:10.1093/oxfordjournals.bjc.a048181](https://doi.org/10.1093/oxfordjournals.bjc.a048181)

Johnson, S. D. (2008). Repeat burglary victimisation: a tale of two
theories. *Journal of Experimental Criminology*, 4(3), 215-240.
[doi:10.1007/s11292-008-9055-3](https://doi.org/10.1007/s11292-008-9055-3)

Farrell, G. and Pease, K. (1993). *Once Bitten, Twice Bitten: Repeat
Victimisation and its Implications for Crime Prevention*. Crime
Prevention Unit Series Paper 46. Home Office, London.

## Examples

``` r
set.seed(8)
hh <- sample(1:400, 300, TRUE, prob = rgamma(400, 0.5))
day <- sort(sample(1:365, 300, TRUE))
rptvic(hh, t = day, pop = 1000, brk = c(0, 7, 30, 90, 365))
#> rptvic: 300 incidents, 155 victims, 145 repeats (48.3%)
#>  k victims incidents     expected
#>  0     845         0 7.408182e+02
#>  1      77        77 2.222455e+02
#>  2      42        84 3.333682e+01
#>  3      19        57 3.333682e+00
#>  4       8        32 2.500261e-01
#>  5       6        30 1.500157e-02
#>  6       1         6 7.500784e-04
#>  7       2        14 3.214622e-05
#>  lo  hi  n      share        cum
#>   0   7  9 0.06206897 0.06206897
#>   7  30 26 0.17931034 0.24137931
#>  30  90 58 0.40000000 0.64137931
#>  90 365 52 0.35862069 1.00000000
```
