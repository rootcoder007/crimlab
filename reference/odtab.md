# Tables of the open data service

`odls()` lists the tables of the project's open data service and
`odtab()` loads one by its key, as the service serves it at the time of
the call. Everything about a table (its address, rows, columns and
SHA-256) is read from the service's manifest at run time, so a table
that its publisher has replaced is loaded as it is now, with no new
release of crimlab.

## Usage

``` r
odtab(
  key,
  cols = NULL,
  refresh = FALSE,
  base = getOption("crimlab.data", "https://data.crimlab.dev")
)

odls(base = getOption("crimlab.data", "https://data.crimlab.dev"))
```

## Arguments

- key:

  A table's key, as `odls()` lists it.

- cols:

  The names of the columns to keep, or `NULL` (all).

- refresh:

  `TRUE` downloads the table again even when the cache holds a copy that
  matches.

- base:

  The address of the service, or of a mirror of it (a folder holding
  `manifest.json` and `<db>/<table>.csv.gz`, which may be a `file://`
  address). The default is `getOption("crimlab.data")`, or the project's
  service when that is not set.

## Value

`odls()` a data frame with a row per table: `key`, `rows`, `bytes` (of
the gzip file), `sha256` (`NA` for a table that cannot be checked) and
`columns`, a list of the column names. `odtab()` the table as a data
frame.

## Details

A downloaded file is used only when its SHA-256 is the one the manifest
lists. When it is not, the manifest is read once more (the table may
have been replaced between the two requests) and, if the file still does
not match, it is deleted and `odtab()` stops. A table the manifest lists
without a SHA-256 cannot be checked and is not loaded.

Tables are kept in `tools::R_user_dir("crimlab", "cache")`. A kept file
is used again only while its SHA-256 is the one the manifest lists now;
otherwise, or with `refresh = TRUE`, the table is downloaded again.

The columns are those of the file, named as served. `cols` picks some of
them, matched without regard to case or underscores (`"case_id"` finds
`"CaseID"`), so that a publisher's respelling of a name does not break a
call; an exact name is preferred.

The values are those the file holds. Every column is read as text; an
empty field and the text `NA` are missing values. A column becomes
numbers only when every value in it is a plain decimal number that a
number holds exactly as written: digits with at most a leading minus, a
decimal point and an exponent, no leading zero, at most 15 significant
digits. So identifiers such as `"007"` (not the same as `"7"`), long
identifiers, dates and codes such as `"F"` stay text, and no column
becomes logical. A column of whole numbers within R's integer range is
integer, any other numeric column double.

Errors: a key the manifest does not list; a table listed without a
SHA-256; a file that does not match it; a line with more or fewer fields
than the header; a table whose number of rows is not the one the
manifest lists; a column the table does not have, or a name that matches
more than one; a service that cannot be reached, answers with an error
or does not answer in time (60 seconds for the manifest, 600 for a
table, or `getOption("timeout")` if longer).

## Examples

``` r
if (FALSE) { # \dontrun{
# needs the internet
tabs <- odls()
head(tabs[, c("key", "rows")])
x <- odtab(tabs$key[!is.na(tabs$sha256)][1L])
str(x)
} # }
```
