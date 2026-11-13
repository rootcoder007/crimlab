# Little's law on a court docket

From the arrival and disposition times of cases (a disposition missing
or after the horizon: still pending), the time-average backlog equals
the filing rate times the mean time to disposition when every case has
been disposed (Little 1961): `resid` is that identity's residual, 0 up
to rounding. The backlog counts each pending case to the horizon; the
mean time is over the disposed cases, which understates the truth when
cases are pending (see
[`losbnd()`](https://rootcoder007.github.io/crimlab/reference/losbnd.md)).
Given a target backlog, `need` is the mean time to disposition that
holds it at the present filing rate. See `proof("bklog")`.

## Usage

``` r
bklog(arrive, dispose, horizon = NULL, target = NULL)
```

## Arguments

- arrive:

  The arrival time of each case: numbers from 0.

- dispose:

  The disposition time of each case: numbers, `NA` for one pending.

- horizon:

  `NULL` (the latest arrival or disposition) or the end of the window:
  one positive number.

- target:

  `NULL`, or a target backlog: one number from 0.

## Value

An object of class `bklog`, a list: `est`, a named vector (`cases`,
`pending`, `rate` the filing rate, `wait` the mean time to disposition,
`backlog`, `occup` the case-time, `resid`); `horizon`; and `need` (with
`target`).

## References

Little, J. D. C. (1961). A proof for the queuing formula: L = lambda W.
*Operations Research*, 9(3), 383-387.
[doi:10.1287/opre.9.3.383](https://doi.org/10.1287/opre.9.3.383)

## Examples

``` r
bklog(arrive = c(0, 2, 3, 7, 9), dispose = c(5, 4, 12, NA, 15),
      horizon = 20, target = 2)
#> bklog: Little's law over a horizon of 20 
#>   cases pending    rate    wait backlog   occup   resid 
#>    5.00    1.00    0.25    5.50    1.75   35.00      NA 
```
