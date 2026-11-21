# Fetch one police oversight report

Downloads a director's report of Ontario's Special Investigations Unit
by its report number, as the page the Unit serves, to a file.
[`ovrtxt()`](https://rootcoder007.github.io/crimlab/reference/ovrtxt.md)
then gives its text and
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
codes it into fields.

## Usage

``` r
ovrget(drid, dest, lang = c("en", "fr"), over = FALSE, wait = 120)
```

## Arguments

- drid:

  The report number: one positive whole number.

- dest:

  The file to write.

- lang:

  `"en"` (the default) or `"fr"`.

- over:

  `TRUE` replaces `dest` when it exists.

- wait:

  The seconds to wait for the file (`getOption("timeout")` when that is
  longer).

## Value

As
[`srcget()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md):
invisibly, a list of `path`, `bytes`, `sha256` and `url`.

## Details

The report number is the `drid` in the address of a report's page. The
Unit publishes each report in English and in French.

Errors: `drid` not one positive whole number; `lang` not `"en"` or
`"fr"`; those of
[`srcget()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md).

## See also

[`ovrtxt()`](https://rootcoder007.github.io/crimlab/reference/ovrtxt.md),
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md);
[`srcls()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
for other sources.

## Examples

``` r
if (FALSE) { # \dontrun{
# needs the internet
f <- ovrget(1234, tempfile(fileext = ".html"))
substr(ovrtxt(readLines(f$path, warn = FALSE)), 1, 200)
} # }
try(ovrget(0, tempfile()))
#> Error : ovrget: `drid` must be one positive report number (the drid of the report's address).
```
