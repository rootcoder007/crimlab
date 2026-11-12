# Benchmark decomposition of a disparity

Splits each group's disparity in force per resident, against a reference
group, into a disparity in contact per resident times a disparity in
force per contact. The split is a product, never a sum. See
`proof("dspbm")`.

## Usage

``` r
dspbm(pop, contact, force, ref = 1)
```

## Arguments

- pop:

  Residents by group: positive numbers, as a numeric vector, a 1-d
  table, a list of single numbers or a one-row data frame.

- contact:

  Police contacts by group: positive numbers, in the forms of `pop`;
  named inputs are matched to `pop` by name, unnamed ones by position
  (only names protect against a mismatch).

- force:

  Uses of force by group: positive numbers, in the forms of `pop`,
  matched by name.

- ref:

  The reference group, by name (matched exactly, case included; a factor
  by its label) or by position in `pop` as given.

## Value

A data frame with one row per group, in the order of `pop`: `group`,
`ref` (`TRUE` for the reference group), `cd` (contact per resident),
`fd` (force per contact) and `rd` (force per resident), each against the
reference group, with `rd = cd * fd`.

## Examples

``` r
dspbm(pop = c(a = 9000, b = 1000), contact = c(a = 900, b = 300),
      force = c(a = 9, b = 6))
#>   group   ref cd fd rd
#> 1     a  TRUE  1  1  1
#> 2     b FALSE  3  2  6
```
