# Person-year table of a restrictive confinement release

Reduces a placement-level release of Ontario's restrictive confinement
records (table a01 of the open data, one row per person, fiscal year,
pair of regions and alert configuration, with a count of placements) to
one row per person and fiscal year, carrying placement volatility
([`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md)),
alert complexity
([`alrtcx()`](https://rootcoder007.github.io/crimlab/reference/alrtcx.md))
and the covariates of Ruhela (2026): gender, age group, fiscal year and
regional cluster (the regions the person-year's records name, sorted and
joined by `" + "`).

## Usage

``` r
otpyr(
  x,
  cols = NULL,
  years = NULL,
  ord = "region",
  keep = c("100", "010", "110", "011", "111", "000"),
  base = getOption("crimlab.data", "https://data.crimlab.dev")
)
```

## Arguments

- x:

  A data frame, or the key of a table of the open data service.

- cols:

  The column names: a named character vector with any of `year`, `id`,
  `ra`, `rb`, `gender`, `age`, `mh`, `sr`, `sw`, `np` (the placements of
  the row), overriding the defaults.

- years:

  `NULL` for every fiscal year in `x`, or the years to keep.

- ord:

  The order of the records for
  [`plcvol()`](https://rootcoder007.github.io/crimlab/reference/plcvol.md).

- keep:

  The alert configurations counted by
  [`alrtcx()`](https://rootcoder007.github.io/crimlab/reference/alrtcx.md).

- base:

  The data service, for a key (see
  [`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)).

## Value

A data frame of class `otpyr`, one row per person-year: `id`, `year`,
`n` (rows), `np` (placements), `within`, `across`, `vm`, `ac`, `cx`,
`anyal`, `gender`, `age`, `rc`; attributes `years`, `drop`, `ord`.

## Details

`x` is the table, or a key of the open data service read by
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md).
Columns are found by `cols`, whose defaults are the names the release
has always used, matched without regard to case or underscores, so the
publisher's `end_fiscal_year` finds `EndFiscalYear`. The fiscal years
are those in the data; `years` keeps some of them. Identifiers are
reassigned every fiscal year in the release, so the person-year is the
unit, never the person.

Rows missing any column used are dropped and counted (attribute `drop`);
an alert flag that is neither yes nor no (nor 1 or 0) is a missing one.
Alert complexity counts the six configurations the release shows
(`keep`), as the published analysis does.

## References

Ruhela, V. S. (2026). Alert complexity and placement volatility in
Ontario restrictive confinement data. Major research paper, Centre for
Criminology and Sociolegal Studies, University of Toronto.

## Examples

``` r
d <- data.frame(
  EndFiscalYear = c(2025, 2025, 2025, 2026, 2026),
  UniqueIndividual_ID = c("a", "a", "b", "a", "c"),
  Region_AtTimeOfPlacement = c("East", "West", "East", "North", "East"),
  Region_MostRecentPlacement = c("East", "West", "East", "East", "East"),
  Gender = "Male", Age_Category = "25 to 49",
  MentalHealth_Alert = c("Yes", "No", "No", "No", "Yes"),
  SuicideRisk_Alert = "No", SuicideWatch_Alert = "No",
  Number_Of_Placements = c(3, 1, 2, 4, 1)
)
otpyr(d)
#> otpyr: 4 person-years, fiscal years 2025-2026 (2)
#>  year py  cx  vm
#>  2025  2 0.5 0.5
#>  2026  2 0.0 0.5
```
