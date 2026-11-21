# Ask whether a newer release exists

Asks the project's release page once, now, for the number of the latest
release, and says whether it is newer than the installed package.
Nothing is sent but that request (a limit of two seconds), and nothing
is kept.

## Usage

``` r
updchk()
```

## Value

Invisibly, a list: `current` (the installed version), `latest` (the
latest release, `NA` when the page could not be reached) and `newer`
(`TRUE` when `latest` is above `current`). One line is printed.

## Details

The package never asks on its own: no request is made when it is
attached unless a notice has been turned on (see the section "Notice of
a newer release" of
[crimlab-package](https://rootcoder007.github.io/crimlab/reference/crimlab-package.md)).

## Examples

``` r
if (FALSE) { # \dontrun{
# needs the internet
updchk()
} # }
```
