# The tables of the data service about Indigenous peoples

Lists the tables of the project's open data service drawn from
Statistics Canada's tables on Indigenous peoples and from the 2021
Census Indigenous Population Profile, with their descriptions, sizes and
licences, optionally those whose description matches a pattern.
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
loads one by its key.

## Usage

``` r
indtabs(
  pat = NULL,
  base = getOption("crimlab.data", "https://data.crimlab.dev")
)
```

## Arguments

- pat:

  A regular expression matched against the descriptions without regard
  to case, or `NULL` (all tables).

- base:

  The address of the service, or of a mirror of it (a folder holding
  `manifest.json` and `<db>/<table>.csv.gz`, which may be a `file://`
  address). The default is `getOption("crimlab.data")`, or the project's
  service when that is not set.

## Value

A data frame with a row per table: `key`, `rows`, `description`,
`licence` and `source`.

## Details

Only openly licensed tables of aggregates are mirrored. The data of the
First Nations Information Governance Centre (the First Nations Regional
Health Survey and its companions) are governed by the OCAP principles of
ownership, control, access and possession, and are not; they are asked
for from the Centre.

Errors: `pat` not one non-empty string or not a valid regular
expression; those of
[`odls()`](https://rootcoder007.github.io/crimlab/reference/odtab.md).

## References

Statistics Canada. *Indigenous Population Profile, 2021 Census of
Population*. Ottawa.

Statistics Canada. *Statistics Canada Open Licence*.
<https://www.statcan.gc.ca/en/reference/licence>

## See also

[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
and
[`odls()`](https://rootcoder007.github.io/crimlab/reference/odtab.md);
[`indiex()`](https://rootcoder007.github.io/crimlab/reference/indiex.md)
for the tables by subject.

## Examples

``` r
if (FALSE) { # \dontrun{
# needs the internet
t <- indtabs("mental health")
t[, c("key", "rows")]
} # }
try(indtabs(pat = ""))
#> Error : indtabs: `pat` must be one non-empty string.
```
