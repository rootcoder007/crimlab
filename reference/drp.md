# Drop empty or constant columns and rows

`drpempty()` removes rows or columns that are entirely missing.
`drpconst()` removes columns that hold a single distinct value. Use them
after reading a spreadsheet export, which often carries blank padding
rows, empty columns and columns that repeat one value (a year, a
jurisdiction) on every row.

## Usage

``` r
drpempty(data, which = c("rows", "cols"), quiet = FALSE)

drpconst(data, naasval = FALSE, quiet = FALSE)
```

## Arguments

- data:

  A data frame whose empty rows or columns are dropped.

- which:

  `"rows"`, `"cols"`, or both (the default, `c("rows", "cols")`).

- quiet:

  Logical: `TRUE` skips the summary of what was dropped (rows and
  columns before and after, the columns removed), which is otherwise
  sent through [`message()`](https://rdrr.io/r/base/message.html);
  default `FALSE`.

- naasval:

  Logical: treat `NA` as a distinct value, so a column of `NA`s plus one
  real value counts as two (default `FALSE`).

## Value

The data frame, with the offending rows or columns removed. The names of
the columns that were dropped are attached as the `"dropped"` attribute
(a character vector, empty when none); dropped rows are not listed.

## Details

A constant column carries no information and breaks anything that scales
by variance, so it is worth removing – but it is also a FINDING. A
column that was informative in an earlier release and is constant in a
fresh one means the source changed, so check that before deleting it and
moving on.

`drpempty()` removes all-`NA` columns first, then all-`NA` rows among
the remaining columns, and resets the row names. An empty string is a
value, not missing. `drpconst()` counts the distinct non-missing values
of each column (with `naasval = TRUE`, `NA` counts as a value), so an
all-`NA` column is also dropped by default. Neither function modifies
anything in place.

Errors: `data` not a data frame (or with missing or duplicated names);
`which` other than `"rows"` and/or `"cols"`.

## See also

[`profcols`](https://rootcoder007.github.io/crimlab/reference/profcols.md),
which reports `n_distinct` without removing anything.

## Examples

``` r
df <- data.frame(
  keep = c(1, 2, NA),
  allna = c(NA, NA, NA),
  constant = c(7, 7, 7)
)

drpempty(df)
#> -- Empty rows and columns dropped
#>    rows 3 -> 3, columns 3 -> 2
#>    columns dropped (all missing): allna
#>   keep constant
#> 1    1        7
#> 2    2        7
#> 3   NA        7
attr(drpempty(df), "dropped")
#> -- Empty rows and columns dropped
#>    rows 3 -> 3, columns 3 -> 2
#>    columns dropped (all missing): allna
#> [1] "allna"

drpconst(df)
#> -- Constant columns dropped
#>    columns 3 -> 1
#>    columns dropped (one value): allna, constant
#>   keep
#> 1    1
#> 2    2
#> 3   NA
attr(drpconst(df), "dropped")
#> -- Constant columns dropped
#>    columns 3 -> 1
#>    columns dropped (one value): allna, constant
#> [1] "allna"    "constant"

# Rows only.
drpempty(data.frame(a = c(1, NA), b = c(2, NA)), which = "rows")
#> -- Empty rows and columns dropped
#>    rows 2 -> 1, columns 2 -> 2
#>   a b
#> 1 1 2

# A column of NAs plus one value is constant by default, and not when
# NA is treated as a value of its own.
x <- data.frame(v = c(NA, NA, 5))
ncol(drpconst(x))
#> -- Constant columns dropped
#>    columns 1 -> 0
#>    columns dropped (one value): v
#> [1] 0
ncol(drpconst(x, naasval = TRUE))
#> -- Constant columns dropped
#>    columns 1 -> 1
#> [1] 1

# only data frames
try(drpempty(c(1, NA)))
#> Error : drpempty: `data` must be a data frame.
```
