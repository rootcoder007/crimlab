# Rates and disparities when exposure is measured by a proxy

Recorded counts (stops, searches, uses of force) are rates per unit of
exposure to police, and exposure is seldom observed. If each group's
true exposure lies within a factor `gamma` of an observed proxy
(residents, drivers on the road, calls for service), each group's true
rate lies within a factor `gamma` of its proxy rate, and the ratio of
two groups' true rates within a factor `gamma^2` of the proxy ratio. The
direction of a disparity is settled only when the proxy ratio lies
outside that band. See `proof("expbnd")`.

## Usage

``` r
expbnd(y, m, gamma, ref = 1)
```

## Arguments

- y:

  Recorded counts by group: positive numbers, as a numeric vector, a 1-d
  table, a list of single numbers or a one-row data frame.

- m:

  The exposure proxy by group: positive numbers, in the forms of `y`;
  named inputs are matched to `y` by name, unnamed ones by position
  (only names protect against a mismatch).

- gamma:

  The factor within which true exposure lies around the proxy: one
  number, at least 1.

- ref:

  The reference group, by name (matched exactly, case included; a factor
  by its label) or by position in `y` as given.

## Value

A data frame with one row per group, in the order of `y`: `group`, `ref`
(`TRUE` for the reference group), `rate` (`y / m`), `lo` and `hi` (the
bounds on the true rate), `ratio` (the proxy rate against the reference
group), `rlo` and `rhi` (the bounds on the true ratio; every end is
reached by some admissible exposure), and `sign`: 1 when the true ratio
exceeds 1 for every admissible exposure, -1 when it is below 1 for every
one, 0 when the data leave it open.

## References

Knox, D., Lowe, W. and Mummolo, J. (2020). Administrative records mask
racially biased policing. *American Political Science Review*, 114(3),
619-637.
[doi:10.1017/S0003055420000039](https://doi.org/10.1017/S0003055420000039)

## Examples

``` r
expbnd(y = c(white = 100, black = 300), m = c(1000, 1000), gamma = 1.5)
#>   group   ref rate         lo   hi ratio       rlo  rhi sign
#> 1 white  TRUE  0.1 0.06666667 0.15     1 0.4444444 2.25    0
#> 2 black FALSE  0.3 0.20000000 0.45     3 1.3333333 6.75    1
```
