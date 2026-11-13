# Bounds implied by rounding and suppression in a published table

A published count stands for every true count that could have produced
it. Rounded to the nearest multiple of `rnd`, a true count lies within
`rnd / 2` of it (randomly rounded, within `rnd - 1`); suppressed below a
limit `k` (a cell written `"<k"`, or a marker such as `"x"` with
`suplim = k`), it lies from 0 to `k - 1`. Given the cells of an earlier
period, `chg` holds every change the two columns allow, and an interval
that spans 0 cannot tell a rise from a fall. See `proof("pbbnd")`.

## Usage

``` r
pbbnd(
  x,
  rnd = NULL,
  kind = "nearest",
  suplim = NULL,
  supmar = c("x", "X", "s", "..", "F", "suppressed", "n/a"),
  prev = NULL,
  pop = NULL,
  prepop = NULL,
  per = 1000
)
```

## Arguments

- x:

  The published cells: numbers from 0, or text holding numbers (with
  thousands separators), `"<k"` and suppression markers.

- rnd:

  `NULL`, or the rounding base: one number above 0.

- kind:

  `"nearest"` or `"random"` rounding.

- suplim:

  `NULL`, or the suppression limit of the cells marked as suppressed
  without one.

- supmar:

  The text markers of a suppressed cell.

- prev:

  `NULL`, or the cells of an earlier period, in the form of `x`.

- pop, prepop:

  `NULL`, or the populations now and before: one number or one per cell,
  above 0.

- per:

  The rates are per `per` people: one positive number.

## Value

An object of class `pbbnd`, a list: `bnd`, a data frame (`value` the
number read, `lower`, `upper`, `status` one of `"exact"`, `"rounded"`,
`"suppressed"` and `"missing"`), and `chg`, with `prev`, a data frame
(`chglo`, `chgup` the change's bounds, `pctlo`, `pctup` the percent
change's, and with populations `ratelo`, `rateup`, `preratelo`,
`prerateup`).

## Examples

``` r
pbbnd(c("1,240", "35", "<5", "x"), rnd = 5, suplim = 10)
#> pbbnd: what the published cells allow
#>  value  lower  upper     status
#>   1240 1237.5 1242.5    rounded
#>     35   32.5   37.5    rounded
#>     NA    0.0    4.0 suppressed
#>     NA    0.0    9.0 suppressed
pbbnd(c(120, 15), rnd = 5, prev = c(100, 20), pop = 5e4)
#> pbbnd: what the published cells allow
#>  value lower upper  status
#>    120 117.5 122.5 rounded
#>     15  12.5  17.5 rounded
#>  chglo chgup     pctlo    pctup ratelo rateup
#>     15    25  14.63415 25.64103   2.35   2.45
#>    -10     0 -44.44444  0.00000   0.25   0.35
```
