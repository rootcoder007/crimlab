# Per-area pollution burden

Applies
[`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
to every row of a table with one row per area (forward sortation area,
census tract, neighbourhood, ...) and returns the rows sorted by
attributable cases, worst first: where in a city or region the pollution
burden falls.

## Usage

``` r
polbrdarea(
  areatab,
  areacol = "fsa",
  expcol = "exposure",
  popcol = "population",
  basratecol = "brate",
  poll = "NO2",
  y = "allcause",
  rare = FALSE,
  hirate = FALSE
)

# S3 method for class 'polbrdarea'
print(x, ...)
```

## Arguments

- areatab:

  A data frame with one row per area: a label, the exposure in
  micrograms per cubic metre, the population in persons and the baseline
  rate per person-year.

- areacol, expcol, popcol, basratecol:

  Column names (single strings) of the label, exposure, population and
  rate; defaults `"fsa"`, `"exposure"`, `"population"` and `"brate"`.

- poll, y, rare, hirate:

  Passed to
  [`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
  for every area.

- x:

  An object of class `polbrdarea`.

- ...:

  Not used.

## Value

A data frame of class `polbrdarea` with one row per area, sorted by
`attrcases`, largest first: the four input columns plus `rr`, `paf`,
`attrcases` and `basecases`.

## Details

Each area is treated as wholly exposed at its own concentration
(`expoprev = 1`), with the pollutant's default reference concentration.
Every area is checked as
[`polbrdn`](https://rootcoder007.github.io/crimlab/reference/polbrdn.md)
checks a single population, and the first area that fails stops the call
with that area's error. Refusals: a column name that is not a single
string or not in the table, an empty table, and a non-numeric exposure,
population or rate column. The area labels can be of any type.

## See also

[`expconcidx`](https://rootcoder007.github.io/crimlab/reference/expconcidx.md)
for how the burden is spread across income.

## Examples

``` r
# simulated areas
t <- data.frame(fsa = c("A1", "B2"), exposure = c(28, 18),
                population = c(40000, 60000), brate = 0.008)
polbrdarea(t)
#> Pollution burden by area: 2 areas, 18.75 attributable of 800 cases (2.34%)
#>  fsa exposure population brate    rr     paf attrcases basecases
#>   A1       28      40000 0.008 1.036 0.03502    11.205       320
#>   B2       18      60000 0.008 1.016 0.01572     7.544       480
try(polbrdarea(t[, c("fsa", "exposure")]))
#> Error : polbrdarea: areatab missing columns: population, brate.
```
