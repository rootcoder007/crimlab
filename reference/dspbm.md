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

  Residents by group: positive numbers, optionally named.

- contact:

  Police contacts by group: positive numbers, same length.

- force:

  Uses of force by group: positive numbers, same length.

- ref:

  The reference group, by name (matched exactly, case included) or
  position.

## Value

A data frame with one row per group: `group`, `cd` (contact per
resident), `fd` (force per contact) and `rd` (force per resident), each
against the reference group, with `rd = cd * fd`.

## Examples

``` r
dspbm(pop = c(a = 9000, b = 1000), contact = c(a = 900, b = 300),
      force = c(a = 9, b = 6))
#>   group cd fd rd
#> 1     a  1  1  1
#> 2     b  3  2  6
```
