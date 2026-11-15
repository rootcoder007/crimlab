# Check a table against a preset before use

`prechk()` reads a table the way
[`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
would under a preset and reports, in one error, every problem that would
stop it or quietly change its result: the columns the preset reads
(directly, in a filter or in a recipe) that the table lacks, the labels
of a column outside the codes the publisher lists for it, and the labels
a mapping recode with no `*` rule does not cover. A table that passes is
returned unchanged, so the check can sit in front of the call:
`clprep(prechk(d, p), p)`.

## Usage

``` r
prechk(data, preset, codes = NULL)
```

## Arguments

- data:

  The table: a data frame.

- preset:

  The preset's name, as
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  lists them, or a preset of one's own laid out as the rows of
  `presets.tsv` (see
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)).

- codes:

  The labels each column may hold: a data frame with columns `column`
  and `value`, one row per label. `NULL` (the default) takes the rows of
  `codes.tsv` whose `agency` is the preset's name up to its first `-`;
  none for a preset of one's own unless given here.

## Value

`data`, invisibly, when it passes; otherwise an error naming each
missing column and each unknown label with its column.

## Details

The codes of the shipped presets are rows of `codes.tsv` (in the
installed package's `presets` folder): `agency`, `column`, `value` and
`source`, one row per label. A preset named `otis-...` is checked
against the rows of agency `otis`: the labels of the data dictionary of
the Ontario Ministry of the Solicitor General for its segregation,
restrictive confinement and deaths in custody tables, one list per
column over all the tables that carry it (so `Gender` allows `All`,
which only the tables of medians and modes use). Labels are compared
without case or surrounding space, as
[`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
compares them, and missing values are not codes
([`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
drops their rows). Columns with no codes listed, counts and days among
them, are not checked for labels.

The presets describe the layouts of published tables, aggregate counts
and de-identified records; crimlab ships none of the tables.

## Examples

``` r
# a table laid out as the published OTIS b03 table, simulated
d <- data.frame(
  EndFiscalYear = rep(2021:2022, each = 4),
  Region_AtTimeOfPlacement = rep(c("Central", "Toronto"), 4),
  Institution_AtTimeOfPlacement = rep(c("Quinte Detention Centre",
                                        "Toronto South Detention Centre"),
                                      4),
  Alert_Type = "Mental Health Alert",
  Alert_Presence = rep(c("Yes", "Yes", "No", "No"), 2),
  Number_SegregationPlacements = c(31, 52, 140, 260, 28, 61, 133, 248)
)
prechk(d, "otis-b03-cnt")
# a region the publisher does not list, and a column missing
bad <- d[-6]
bad$Region_AtTimeOfPlacement[1] <- "Southern"
try(prechk(bad, "otis-b03-cnt"))
#> Error : prechk: the table does not fit preset otis-b03-cnt:
#> * `data` has no column Number_SegregationPlacements.
#> * column Region_AtTimeOfPlacement has labels not among its codes: Southern.
# a preset of one's own, and a label its recode does not cover
mine <- data.frame(preset = "my-cnt", generic = "cntreg",
                   arg = c("y", "grp"), column = c("count", "sex"),
                   recode = c("", "f=1;m=0"))
e <- data.frame(count = c(4, 30, 6, 41), sex = c("f", "m", "f", "x"))
try(prechk(e, mine))
#> Error : prechk: the table does not fit preset my-cnt:
#> * column sex has labels with no recode: x.
```
