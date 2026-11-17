# Alert complexity of each person-year

Counts the distinct alert configurations among a person-year's records.
Each record carries flags (in Ontario's records the mental health,
suicide risk and suicide watch alerts); its configuration is the tuple
of the flags, `(1, 0, 0)` for a mental health alert alone, coded as the
text `"100"`. Alert complexity is the number of distinct configurations
the person-year's records span, and `cx`, complexity of two or more, is
the treatment of Ruhela (2026, Appendix H).

## Usage

``` r
alrtcx(flags, id, keep = NULL)
```

## Arguments

- flags:

  The alert flags of each record: a data frame, list or matrix with one
  column per flag; each value 1 / 0, `TRUE` / `FALSE` or yes / no text
  (case and surrounding space ignored).

- id:

  The person-year of each record, as for
  [`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md).

- keep:

  `NULL` to count every configuration, or the configurations to count,
  as codes such as `"100"` (one digit per flag).

## Value

A data frame of class `alrtcx`, one row per person-year in order of
first appearance: `id`, `n` (records), `ac` (distinct configurations
among those counted), `cx` (1 when `ac >= 2`), and `anyal` (1 when any
record carries an alert).

## Details

`cx = 1` holds exactly when two records are in different configurations:
a record without alerts and one with a single alert give `cx = 1`, while
one record with all three alerts gives `cx = 0`. The count measures
change of alert state over the records, not alerts held at once
(`proof("alrtcx")`). It does not depend on the order of the records or
on an injective relabelling of the configurations, never falls when a
record is added, and is at most the number of configurations allowed.
`keep` restricts the count to the configurations given (the published
analysis counts the six seen in the release and ignores any other). The
work is `O(n log n)`.

## References

Ruhela, V. S. (2026). Alert complexity and placement volatility in
Ontario restrictive confinement data. Major research paper, Centre for
Criminology and Sociolegal Studies, University of Toronto.

## Examples

``` r
f <- data.frame(mh = c("Yes", "No", "No", "Yes"),
                sr = c("No", "No", "Yes", "Yes"),
                sw = c("No", "No", "No", "Yes"))
alrtcx(f, id = c(1, 1, 2, 3))
#> alrtcx: 3 person-years; 1 (33.3%) span two or more alert configurations
#>  id n ac cx anyal
#>   1 2  2  1     1
#>   2 1  1  0     1
#>   3 1  1  0     1
```
