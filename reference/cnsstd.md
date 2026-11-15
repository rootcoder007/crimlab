# An Indigenous standard population from the census

The age structure of a group in a geography of the 2021 Census
Indigenous Population Profile, as the standard of
[`dastd`](https://rootcoder007.github.io/crimlab/reference/dastd.md) or
[`stdprev`](https://rootcoder007.github.io/crimlab/reference/stdprev.md).

## Usage

``` r
cnsstd(
  x,
  grp = "Indigenous identity",
  geo = "Canada",
  ages = c("0 to 4 years", "5 to 9 years", "10 to 14 years", "15 to 19 years",
    "20 to 24 years", "25 to 34 years", "35 to 44 years", "45 to 54 years",
    "55 to 64 years", "65 to 74 years", "75 years and over"),
  resid = "Total - Residence on or off reserve"
)
```

## Arguments

- x:

  A
  [`cnsprof`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md).

- grp:

  The group (default `"Indigenous identity"`).

- geo:

  The geography: its name (default `"Canada"`) or, where names repeat,
  its DGUID.

- ages:

  The age groups, as the profile names them.

- resid:

  On or off reserve, in the tables that split it (default the total).

## Value

A data frame: `age`, `population`, `weight` (shares summing to 1).

## Details

Comparisons of First Nations, Metis and Inuit health are standardised to
the Indigenous population's own, younger, age structure (Kumar and
Tjepkema 2019): a standard as old as Canada's gives weight to ages where
few Indigenous people are. This takes the population in private
households of each age group (the profile's first characteristic) for
the group and geography given. Age groups that do not overlap are chosen
by default (0 to 4, 5 to 9, 10 to 14, 15 to 19, 20 to 24, 25 to 34, 35
to 44, 45 to 54, 55 to 64, 65 to 74, 75 and over); collapse them to the
age groups of the rates being standardised.

## References

Kumar MB, Tjepkema M (2019). Suicide among First Nations people, Metis
and Inuit (2011-2016): findings from the 2011 Canadian Census Health and
Environment Cohort (CanCHEC). Statistics Canada catalogue no.
99-011-X2019001.

## See also

[`dastd`](https://rootcoder007.github.io/crimlab/reference/dastd.md),
[`stdprev`](https://rootcoder007.github.io/crimlab/reference/stdprev.md).

## Examples

``` r
# made-up population counts by age group for two groups, in the layout of a
# profile table
ages <- c("0 to 14 years", "15 to 64 years", "65 years and over")
x <- data.frame(
  GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
  AGE_GROUP = c("Total - Age groups", ages),
  RESIDENCE = "Total - Residence on or off reserve",
  CHARACTERISTIC_ID = 1,
  CHARACTERISTIC_NAME = "Total - Age groups of the population - 100% data",
  `IND_IDENT (2)` = c(1200, 400, 700, 100), SYMBOL_IND_IDENT = "",
  NON_IND = c(5000, 800, 3300, 900), SYMBOL_NON_IND = "",
  check.names = FALSE
)
cnsstd(cnsprof(x), geo = "Lakeside", ages = ages)
#>                 age population     weight
#> 1     0 to 14 years        400 0.33333333
#> 2    15 to 64 years        700 0.58333333
#> 3 65 years and over        100 0.08333333
```
