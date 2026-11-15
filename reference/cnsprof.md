# Read the 2021 Census Indigenous Population Profile

Reads a table of the Indigenous Population Profile, 2021 Census of
Population (Statistics Canada catalogue 98-510-X2021001 to
98-510-X2021016), as published, into one row per geography, gender, age
group, characteristic and group: First Nations, Metis, Inuit, the other
identity groups and non-Indigenous people (or, in the household tables,
Indigenous and non-Indigenous households), each in its own rows.

## Usage

``` r
cnsprof(x)
```

## Arguments

- x:

  A profile table as published (a data frame read from its CSV).

## Value

A data frame of class `cnsprof`: `geo`, `dguid`, `geolvl`, `residence`
(on or off reserve, in the tables that split it), `gender`, `age`,
`chaid`, `characteristic`, `level`, `parid`, `rootid`, `sample`,
`group`, `value`, `quality`, `tnrsf`, `tnrlf`.

## Details

The profile gives each characteristic in a hierarchy shown only by
indentation (“Total - Highest certificate, diploma or degree”, “High
(secondary) school diploma or equivalency certificate”, ...): the level,
the parent and the topic's total (`rootid`) are worked out from it, so a
share can be taken of its parent. A topic marked "25% sample data" comes
from the long-form sample (`sample = "25%"`); the others from the full
count (`"100%"`). Counts are randomly rounded to multiples of 5 (of 10
below 10) by Statistics Canada, and a cell suppressed to protect
confidentiality (`x`), too unreliable (`F`), not available (`..`) or not
applicable (`...`) has no value here; `E` values are kept and marked.
The global non-response rates of the geography (`tnrsf`, `tnrlf`, in
percent) are kept: Statistics Canada suppresses an area whose long-form
rate is 50% or more, and a high rate means a larger error than the
sample alone.

## References

Statistics Canada (2023). Indigenous Population Profile, 2021 Census of
Population. Statistics Canada Catalogue no. 98-510-X2021001. Released
June 21, 2023. Published at
www12.statcan.gc.ca/census-recensement/2021/dp-pd/ipp-ppa/index.cfm?Lang=E
(the site answers scripts with 403); the Wayback Machine holds a copy of
that address captured on 8 March 2026 (web.archive.org, snapshot
20260308015010).

Statistics Canada (2022). Guide to the Census of Population, 2021.
Catalogue no. 98-304-X.

## See also

[`profshr`](https://rootcoder007.github.io/crimlab/reference/profshr.md),
[`cnsstd`](https://rootcoder007.github.io/crimlab/reference/cnsstd.md),
[`cnscwb`](https://rootcoder007.github.io/crimlab/reference/cnscwb.md).

## Examples

``` r
# made-up numbers in the layout of a population profile table, for one
# place and two groups (the published tables hold many more rows)
nm <- c(
  "Total - Age groups of the population - 100% data",
  paste("Number of total income recipients aged 15 years and over in",
        "private households in 2020 - 25% sample data"),
  "Average total income in 2020 among recipients ($)",
  paste("Total - Highest certificate, diploma or degree for the population",
        "aged 25 years and over in private households - 25% sample data"),
  "  No certificate, diploma or degree",
  "  Bachelor's degree or higher",
  "Participation rate",
  "Unemployment rate"
)
x <- data.frame(
  GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
  AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3),
                rep("25 to 64 years", 2)),
  RESIDENCE = "Total - Residence on or off reserve",
  CHARACTERISTIC_ID = 1:8, CHARACTERISTIC_NAME = nm,
  `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14),
  SYMBOL_IND_IDENT = "",
  NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6),
  SYMBOL_NON_IND = "",
  check.names = FALSE
)
p <- cnsprof(x)
p[p$chaid == 5, c("group", "characteristic", "level", "parid", "value")]
#> Indigenous Population Profile, tidied
#> 
#>               group                     characteristic  level  parid  value
#> Indigenous identity  No certificate, diploma or degree      1      4    150
#>      Non-Indigenous  No certificate, diploma or degree      1      4    400
```
