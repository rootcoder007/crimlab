# Check a coded table of oversight reports

Flags the cells of a table coded from oversight reports (by
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md),
by hand, or by any other means) that are not of the form their field
calls for, so that a parsing slip or a typing error is found before the
table is analysed. Nothing is changed.

## Usage

``` r
ovrchk(data, schema = ovrsch(), shape = NULL)
```

## Arguments

- data:

  The coded table: a data frame, one row per report.

- schema:

  The fields, as
  [`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md)
  gives them; only `name` and `rule` are read.

- shape:

  `NULL`, or patterns named by column: a named character vector.

## Value

A data frame with one row per row of `data`: `row`, `issues` (how many
cells are flagged) and `what` (each flagged column with its form,
`"date_of_incident_iso: date"`, separated by `"; "`).

## Details

The form of a column comes from the rule of its field in `schema`: a
`date` field must hold real dates as `YYYY-MM-DD`; a `num`, `tag` or
`before` field whole numbers from 0; a `flag` field `TRUE`, `FALSE`, yes
or no (in any case, or y and n). `shape` adds columns that must follow a
pattern, such as case numbers: in a pattern `9` stands for a digit, `A`
for a capital letter, `a` for a small letter and any other character for
itself, with alternatives separated by `"|"` (`"99-AAA-999|99-AA-999"`).
Columns of `schema` that the table lacks are passed over. A missing
value (`NA` or `""`) is not an error.

## Examples

``` r
d <- data.frame(case_number = c("23-OFD-001", "23-OFD-1", NA),
                date_of_incident_iso = c("2023-01-05", "2023-02-30", ""),
                number_of_subject_officials = c(1, 2.5, NA),
                charges_recommended = c(FALSE, NA, TRUE))
ovrchk(d, shape = c(case_number = "99-AAA-999|99-AA-999"))
#>   row issues
#> 1   1      0
#> 2   2      3
#> 3   3      0
#>                                                                                 what
#> 1                                                                                   
#> 2 date_of_incident_iso: date; number_of_subject_officials: count; case_number: shape
#> 3                                                                                   
```
