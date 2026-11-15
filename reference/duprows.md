# Duplicated rows, with their groups

Returns the rows that share a combination of `cols` with at least one
other row, grouped so the duplicates sit together. A duplicate is
usually a join that fanned out or a re-release appended rather than
replaced, both of which are visible only once the offending rows are in
front of you.

## Usage

``` r
duprows(data, cols = NULL, quiet = FALSE)
```

## Arguments

- data:

  A data frame whose rows are checked for duplicates.

- cols:

  Character vector of the columns defining a duplicate, or `NULL` (the
  default) for all of them.

- quiet:

  Logical: `TRUE` skips the count of duplicate rows and groups, which is
  otherwise sent through
  [`message()`](https://rdrr.io/r/base/message.html); default `FALSE`.

## Value

A data frame: the duplicated rows (every occurrence, including the
first), ordered so each group sits together, with an integer
`dupe_count` column giving each group's size and row names reset. Zero
rows (with the same columns) when there are none.

## Details

Two rows are duplicates when their values, converted to text, are equal
in every column of `cols`; two `NA`s are equal. Columns are taken by
position when `cols = NULL`, so a table with repeated column names can
still be inspected. Errors when `data` is not a data frame, a named
column does not exist, or a named column's name is shared by two
columns.

## See also

[`misssmry`](https://rootcoder007.github.io/crimlab/reference/misssmry.md)

## Examples

``` r
df <- data.frame(
  id = c(1, 2, 2, 3, 3, 3),
  value = c("a", "b", "b", "c", "d", "c")
)

# Duplicated on every column.
duprows(df)
#> -- Duplicated rows
#>    4 row(s) of 6
#>   id value dupe_count
#> 1  2     b          2
#> 2  2     b          2
#> 3  3     c          2
#> 4  3     c          2

# Duplicated on the identifier alone, which catches more.
duprows(df, "id")
#> -- Duplicated rows
#>    5 row(s) of 6
#>   id value dupe_count
#> 1  2     b          2
#> 2  2     b          2
#> 3  3     c          3
#> 4  3     d          3
#> 5  3     c          3

# No duplicates gives zero rows, not an error.
duprows(data.frame(x = 1:3))
#> -- Duplicated rows
#>    0 row(s) of 3
#> [1] x          dupe_count
#> <0 rows> (or 0-length row.names)

# a column that does not exist is an error
try(duprows(df, "code"))
#> Error : duprows: no such column(s): code.
```
