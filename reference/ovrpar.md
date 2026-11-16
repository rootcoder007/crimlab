# Code police oversight reports into fields

Reads the fields of a schema from each report: the dates, the counts of
officials and witnesses, the service, the outcome and whatever else the
schema names, each by its rule (see
[`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md)).
The rules are deterministic and offline, so the same report always gives
the same row, and a folder of saved reports becomes a data set whose
every value can be traced to the words that gave it.

## Usage

``` r
ovrpar(x, schema = ovrsch(), html = NA)
```

## Arguments

- x:

  The reports: their pages (HTML) or texts, a character vector.

- schema:

  The fields: a data frame with the character columns `name`, `rule`,
  `cue`, `alt` and `sect`, as
  [`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md)
  gives them.

- html:

  `NA` (the default) to read each element as HTML when it holds markup
  (a `<` followed by a letter, `/` or `!`), `TRUE` or `FALSE` to say.

## Value

A data frame with one row per report and one column per field, in the
order of `schema`; a row of `NA` for a report that is `NA`.

## Details

Each report is first made text as
[`ovrtxt()`](https://rootcoder007.github.io/crimlab/reference/ovrtxt.md)
makes it, then each field is read by its rule; a field the report does
not state is `NA`. A field of the rules `num`, `tag` and `before` is a
whole number, one of the rule `flag` `TRUE` or `FALSE`, one of the rule
`date` a `YYYY-MM-DD` date (as text), and the others text. Nothing is
fetched.

## Examples

``` r
r <- paste(
  "<h1>Director's Report - Case # 23-OFD-001</h1>",
  "<p>On January 6, 2023, the Barrie Police Service contacted the SIU.",
  "In the night, a 34-year-old man was injured.</p>",
  "<p>Number of SIU Investigators assigned: 3</p>",
  "<p>WO #1 Interviewed</p><p>WO #2 Interviewed</p>",
  "<p>Joseph Martino</p><p>Director</p>"
)
p <- ovrpar(r)
p[, c("case_number", "police_service", "date_siu_notified_iso",
      "siu_investigators", "number_of_witness_officials",
      "age_affected", "directors_name")]
#>   case_number        police_service date_siu_notified_iso siu_investigators
#> 1  23-OFD-001 Barrie Police Service            2023-01-06                 3
#>   number_of_witness_officials age_affected directors_name
#> 1                           2           34 Joseph Martino
```
