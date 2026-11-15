# Carbon intensity of electricity by location, offline, anywhere

The grams of CO2-equivalent emitted per kilowatt-hour of electricity at
a location, with the year and source of the figure. The table is
compiled into the package, so it works offline in any country: every
country at the latest year published by Our World in Data from Ember's
yearly electricity data (the latest year each grid reports: 2025 or 2024
for nearly all, older for a few grids, Ukraine (2022) among them, named
in `year`), sub-national zones from Electricity Maps' 2024 yearly data
as redistributed in the Green Algorithms data release v3.1 (Canadian
provinces and territories, US balancing authorities, Australian states,
Indian and Japanese regions, Brazilian and Chilean systems, ...), and
the world average. All figures are lifecycle intensities (generation
plus upstream), the basis both sources publish.

## Usage

``` r
co2int(loc = "WORLD")

co2inttab()

# S3 method for class 'co2int'
print(x, ...)
```

## Arguments

- loc:

  A location code, alpha-3 code or name; a number in gCO2e/kWh, returned
  as a user-supplied value; or `NULL` to detect the location from the
  environment, time zone or locale.

- x:

  A result of `co2int()`.

- ...:

  Not used.

## Value

`co2int()`: a list of class `co2int` with `g_per_kwh`, `location`,
`name`, `year`, `basis`, `source` and `note` (a character vector, empty
unless a fallback or detection was involved). `co2inttab()`: the full
table as a data frame with the columns `location`, `iso3`, `name`,
`g_per_kwh`, `year`, `basis` and `source`.

## Details

A location is matched by its code (`"CA"`, `"CA-ON"`, `"US-CAL-CISO"`),
its ISO 3166-1 alpha-3 code (`"CAN"`) or its name (`"Canada"`),
case-insensitively. An unknown sub-national code falls back to its
country and says so in `note`. `NULL` detects the location with
[`dtctloc`](https://rootcoder007.github.io/crimlab/reference/dtctloc.md).

Computing has a footprint of its own, and research that reports one (as
journals and funders increasingly ask) needs the grid's intensity:
[`co2fp`](https://rootcoder007.github.io/crimlab/reference/co2fp.md)
multiplies a computation's energy by it.

## Provenance of the table

The country rows are Ember's yearly electricity data as published by Our
World in Data; the sub-national zones and the Canadian provinces come
from the Electricity Maps yearly files redistributed in the Green
Algorithms data repository (v3.1). Every row names its source and year.
Values are reproduced unchanged: Ember reports three small diesel grids
(Montserrat, Saint Helena, Uzbekistan) at exactly 1000 gCO2e/kWh and
Turkmenistan's gas fleet at 1306, and the table does not round or cap
them. Kosovo appears once, under the code `XK`.

Errors: a location the table does not know; a number that is not a
single non-negative value; text that is not valid in its encoding.

## References

Ember (2026). Yearly Electricity Data, via Our World in Data.

Electricity Maps (2025). 2024 Yearly Carbon Intensity Data, via Green
Algorithms data v3.1.

Lannelongue, L., Grealey, J. and Inouye, M. (2021). Green Algorithms:
quantifying the carbon footprint of computation. Advanced Science 8(12),
2100707.
[doi:10.1002/advs.202100707](https://doi.org/10.1002/advs.202100707)

## See also

[`co2fp`](https://rootcoder007.github.io/crimlab/reference/co2fp.md),
[`dtctloc`](https://rootcoder007.github.io/crimlab/reference/dtctloc.md).

## Examples

``` r
co2int("CA-ON")
#> Carbon intensity of electricity: 90.97 gCO2e/kWh
#>   Location: CA-ON (Canada, Ontario), 2024
#>   Basis: lifecycle
#>   Source: Electricity Maps (2025). Canada 2024 Yearly Carbon Intensity Data (Version January 27, 2025). Electricity Maps (via Green Algorithms data v3.1)
co2int("Canada")$year
#> [1] 2025
co2int("CA-XX")$note
#> [1] "no row for zone CA-XX; its country CA used"
co2int(120)$source
#> [1] "user-supplied carbon intensity"
nrow(co2inttab())
#> [1] 378
# the Canadian provinces in the table, lowest intensity first
tab <- co2inttab()
ca <- tab[
  startsWith(tab$location, "CA-"),
  c("location", "name", "g_per_kwh", "year")
]
head(ca[order(ca$g_per_kwh), ])
#>     location                              name g_per_kwh year
#> 236    CA-NL Canada, Newfoundland and Labrador     27.60 2024
#> 227    CA-QC                    Canada, Quebec     35.04 2024
#> 322    CA-MB                  Canada, Manitoba     52.89 2024
#> 367    CA-BC          Canada, British Columbia     70.34 2024
#> 241    CA-PE      Canada, Prince Edward Island     70.94 2024
#> 313    CA-ON                   Canada, Ontario     90.97 2024
# a location the table does not know is an error, not a silent world
# average
try(co2int("Atlantis"))
#> Error : co2int: Unknown location 'Atlantis'. Use a code, alpha-3 code or name from co2inttab(), or pass a number in gCO2e/kWh.
```
