# Placement volatility of each person-year

Counts how often a person's placement region changes within a period (a
fiscal year), from records that each carry the region at the time of
placement `ra` and the region of the most recent placement `rb`. For a
person-year with records `r = 1, ..., R` in a given order the count is
\$\$v = \sum_r 1(A_r \ne B_r) + \sum\_{r \> 1} 1(A_r \ne A\_{r-1}),\$\$
a change within a record (`within`) plus a change between consecutive
records (`across`), as Ruhela (2026) defines placement volatility for
Ontario's restrictive confinement records (Appendix H).

## Usage

``` r
plcvol(ra, rb, id, ord = c("region", "given"))
```

## Arguments

- ra, rb:

  The region at the time of placement and the region of the most recent
  placement of each record: labels (text, factors or numbers), compared
  as text.

- id:

  The person-year of each record: one label, or a list of labels (person
  and year) combined.

- ord:

  `"region"` or `"given"`, the order of the records within a
  person-year.

## Value

A data frame of class `plcvol`, one row per person-year in the order of
first appearance: `id` (the labels, one column each when `id` is a
list), `n` (records), `within`, `across` and `vm` (their sum).

## Details

The first sum does not depend on the order of the records, the second
does. A release without placement dates fixes no order, so `ord` says
which is used: `"region"` sorts each person-year's records by `ra` (the
order of the published analysis code), under which `across` is the
number of distinct regions at placement less one and no longer depends
on how the rows arrived; `"given"` keeps the rows in the order given.
See `proof("plcvol")`: `0 <= v <= 2R - 1`, with equality possible;
`v = 0` when every record stays in one region; the count is unchanged by
an injective relabelling of the regions and can fall when two regions
are merged.

Rows missing `id`, `ra` or `rb` are dropped. The work is one ordering of
the rows, `O(n log n)` for `n` rows, then one pass.

## References

Ruhela, V. S. (2026). Alert complexity and placement volatility in
Ontario restrictive confinement data. Major research paper, Centre for
Criminology and Sociolegal Studies, University of Toronto.

## Examples

``` r
ra <- c("East", "East", "West", "North", "North")
rb <- c("East", "West", "West", "North", "East")
plcvol(ra, rb, id = c(1, 1, 1, 2, 2))
#> plcvol: 2 person-years, records in region order; mean volatility 1.5 (within 1, across 0.5)
#>  id n within across vm
#>   1 3      1      1  2
#>   2 2      1      0  1
plcvol(ra, rb, id = c(1, 1, 1, 2, 2), ord = "given")
#> plcvol: 2 person-years, records in given order; mean volatility 1.5 (within 1, across 0.5)
#>  id n within across vm
#>   1 3      1      1  2
#>   2 2      1      0  1
```
