# Synthetic example data in the shape of published tables

Draws a data frame with the columns of a published administrative table
and values of the kinds it holds, for examples, tests and teaching where
the real table cannot be shipped or fetched. Each column is drawn on its
own, uniformly from its values: nothing in the result is a record of
anyone, and no relation between columns is kept.

## Usage

``` r
synfix(kind = NULL, n = 100, seed = 1)
```

## Arguments

- kind:

  The name of a preset or of a family of presets, a named list of the
  values of each column, or `NULL` for the list of presets.

- n:

  The number of rows: one whole number from 1 to 1e7.

- seed:

  The seed: one whole number from 0 to 2147483647.

## Value

A data frame of `n` rows; for a family, a named list of them (named
`"a01"`, ..., or `"incidents"`, ...); for `kind = NULL`, a data frame
with one row per column of every preset: `kind`, `column`, `type`
(`"labels"`, `"range"` or `"id"`) and `values` (the labels separated by
`"|"`, or the ends of the range).

## Details

A preset gives the columns of one table: the Ontario corrections tables
on segregation, restrictive confinement and deaths in custody
(`"otis-a01"`, `"otis-b01"` to `"otis-b09"`, `"otis-c01"` to
`"otis-c12"`, `"otis-d01"` to `"otis-d07"`) and on use of force in
correctional institutions (`"cuof-incidents"`, `"cuof-incident_type"`,
`"cuof-indigenous"`, `"cuof-inmate_incident"`,
`"cuof-inmate_participant"`, `"cuof-institution_summary"`,
`"cuof-location_summary"`, `"cuof-race"`, `"cuof-ethnic_origin"`,
`"cuof-religion"`, `"cuof-select_incident_summary"`,
`"cuof-staff_incident"`). A column of a preset takes one of a set of
labels, a whole number from a range, or a synthetic identifier
(`"syn-00001"` to `"syn-"` and `max(2, n %/% 4)` in five digits, so that
identifiers repeat as people's do). `"otis"` or `"cuof"` alone gives the
list of every table of the family, the `i`-th drawn with the seed
`(seed + i) %% 2^31`. `kind = NULL` lists the presets.

Any other table is a named list, one element per column, each a vector
of the values the column takes (labels or numbers), drawn with
replacement; a value given once is a constant column.

The draws are those of R's `sample.int(m, n, replace = TRUE)` after
`set.seed(seed)`, made by the package's own copy of R's generator, so
the result is the same on every system and in every language, and the
session's random numbers are not touched.

## Examples

``` r
head(synfix("otis-b02", n = 20, seed = 3))
#>   EndFiscalYear UniqueIndividual_ID Gender Region_MostRecentPlacement
#> 1          2022           syn-00002 Female                    Toronto
#> 2          2019           syn-00001 Female                    Central
#> 3          2021           syn-00004 Female                    Western
#> 4          2024           syn-00004 Female                    Eastern
#> 5          2021           syn-00003   Male                    Western
#> 6          2019           syn-00002   Male                    Eastern
#>   Age_Category TotalAggregatedDays_Segregation
#> 1          50+                              57
#> 2     18 to 24                              76
#> 3     18 to 24                               8
#> 4     25 to 49                              49
#> 5          50+                              54
#> 6     18 to 24                              44
x <- synfix("cuof", n = 5)
names(x)
#>  [1] "incidents"               "incident_type"          
#>  [3] "indigenous"              "inmate_incident"        
#>  [5] "inmate_participant"      "institution_summary"    
#>  [7] "location_summary"        "race"                   
#>  [9] "ethnic_origin"           "religion"               
#> [11] "select_incident_summary" "staff_incident"         
synfix(list(year = 2019:2023, sex = c("F", "M")), n = 4, seed = 9)
#>   year sex
#> 1 2021   F
#> 2 2023   M
#> 3 2021   M
#> 4 2021   F
head(synfix())
#>       kind                     column   type
#> 1 otis-a01              EndFiscalYear  range
#> 2 otis-a01        UniqueIndividual_ID     id
#> 3 otis-a01   Region_AtTimeOfPlacement labels
#> 4 otis-a01 Region_MostRecentPlacement labels
#> 5 otis-a01                     Gender labels
#> 6 otis-a01               Age_Category labels
#>                                     values
#> 1                                2018|2024
#> 2                                         
#> 3 Central|Eastern|Northern|Toronto|Western
#> 4 Central|Eastern|Northern|Toronto|Western
#> 5                              Female|Male
#> 6                    18 to 24|25 to 49|50+
```
