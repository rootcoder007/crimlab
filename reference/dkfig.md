# The dark figure of crime: bounds instead of a point

What recorded crime and a victimisation survey together say about the
crime that is never recorded, as intervals that hold under stated limits
on the errors (Manski 2003).

## Usage

``` r
dkfig(
  v = NULL,
  r = NULL,
  amax = 0,
  bmax = 0,
  thr = NULL,
  va = NULL,
  ra = NULL,
  rb = NULL,
  lam = 1,
  k = NULL
)
```

## Arguments

- v:

  `NULL`, or survey (true) victimisation rates: numbers from 0 to 1.

- r:

  The recorded rates: one number or one per `v`, from 0 to 1.

- amax, bmax:

  The largest false-positive and false-negative rates of the survey:
  numbers from 0, with a sum below 1.

- thr:

  `NULL`, or a threshold above every `v`.

- va, ra, rb:

  `NULL`, or the survey's offences of the covered type, the recorded
  ones of that type (both above 0, `ra` at most `va`) and the recorded
  offences of the uncovered type (from 0).

- lam:

  The least ratio of the uncovered type's recording rate to the covered
  type's: above 0, at most 1.

- k:

  `NULL`, or the number of offences in each incident: whole numbers from
  1.

## Value

An object of class `dkfig`, a list: `bnd` (with `v`: `v`, `r`, `lo`,
`hi`, `dlo`, `dhi` the dark figure's bounds, `ratio` the upper bound
over the recorded rate, and with `thr` `brk`), `cmp` (with `va`: `rate`
A's recording rate, `rateup`, the total's bounds `totlo`, `tothi`, the
dark figure's `dklo`, `dkhi` and `width`) and `hier` (with `k`: `incid`,
`offen`, `ratio`, `extra`, `lower`, `upper`); the parts not asked for
are `NULL`.

## Details

- With a survey rate `v` and a recorded rate `r`, and survey false
  positives at most `amax` and false negatives at most `bmax`, the true
  rate lies in
  `[max(r, (v - amax) / (1 - amax)), min(1, v / (1 - bmax))]` and the
  dark figure (true minus recorded) in that interval shifted by `r`.
  Given `thr`, `brk` is the largest `bmax` under which the conclusion
  "the true rate is below `thr`" survives: `1 - v / thr`.

- When the survey covers only some offences (type A, `va` of them, `ra`
  recorded) and the recorded count of the rest is `rb`, recorded at no
  less than `lam` times A's rate, the total and the dark figure have the
  bounds in `cmp`.

- With `k`, the number of offences in each recorded incident, the
  hierarchy rule's arithmetic: one incident in the counts stands for
  between 1 and `max(k)` offences.

See `proof("dkfig")`.

## References

Manski, C. F. (2003). *Partial Identification of Probability
Distributions*. Springer.
[doi:10.1007/b97478](https://doi.org/10.1007/b97478)

## Examples

``` r
dkfig(v = c(0.12, 0.08), r = 0.03, amax = 0.01, bmax = 0.2, thr = 0.15)
#> dkfig: bounds on the dark figure
#>      v    r         lo   hi        dlo  dhi    ratio       brk
#> 1 0.12 0.03 0.11111111 0.15 0.08111111 0.12 5.000000 0.2000000
#> 2 0.08 0.03 0.07070707 0.10 0.04070707 0.07 3.333333 0.4666667
dkfig(va = 900, ra = 300, rb = 450, lam = 0.5)
#> dkfig: bounds on the dark figure
#>         rate       rateup        totlo        tothi         dklo         dkhi 
#>    0.3333333    0.5555556 1350.0000000 3600.0000000  600.0000000 2850.0000000 
#>        width 
#> 2250.0000000 
dkfig(k = c(1, 1, 2, 1, 3))
#> dkfig: bounds on the dark figure
#> incid offen ratio extra lower upper 
#>   5.0   8.0   1.6   0.6   5.0  15.0 
```
