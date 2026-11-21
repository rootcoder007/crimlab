# Statistics Canada tables on Indigenous peoples, by subject

The tables of the project's open data service that hold what is most
often asked for about First Nations people, Metis and Inuit: health
indicators, perceived mental health, wellbeing, suicidal thoughts, the
census counts of population and of households, life expectancy and
infant mortality. `indiex()` finds the tables of one subject among those
[`indtabs()`](https://rootcoder007.github.io/crimlab/reference/indtabs.md)
lists, and with `load = TRUE` loads the first.

## Usage

``` r
indiex(
  name = c("health", "mentalhealth", "wellbeing", "suicide", "censuspop", "censushh",
    "lifeexp", "infantmort"),
  load = FALSE,
  base = getOption("crimlab.data", "https://data.crimlab.dev")
)
```

## Arguments

- name:

  The subject: `"health"`, `"mentalhealth"`, `"wellbeing"`, `"suicide"`,
  `"censuspop"`, `"censushh"`, `"lifeexp"` or `"infantmort"`.

- load:

  `TRUE` loads the first table found, with
  [`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md).

- base:

  The address of the service, or of a mirror of it (a folder holding
  `manifest.json` and `<db>/<table>.csv.gz`, which may be a `file://`
  address). The default is `getOption("crimlab.data")`, or the project's
  service when that is not set.

## Value

The rows of
[`indtabs()`](https://rootcoder007.github.io/crimlab/reference/indtabs.md)
for the subject; with `load = TRUE`, the first of those tables as a data
frame.

## Details

The package holds no data: the tables are Statistics Canada's, mirrored
by the service under the Statistics Canada Open Licence and loaded and
checked by
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md).
A subject is matched by words in a table's description, so a table the
service adds later is found without a new release; a subject with no
table there gives a data frame of no rows.

Errors: `name` not one of the subjects; `load` not `TRUE` or `FALSE`;
with `load`, a subject with no table; those of
[`indtabs()`](https://rootcoder007.github.io/crimlab/reference/indtabs.md)
and
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md).

## See also

[`indtabs()`](https://rootcoder007.github.io/crimlab/reference/indtabs.md),
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md);
[`sctidy()`](https://rootcoder007.github.io/crimlab/reference/sctidy.md)
and
[`scests()`](https://rootcoder007.github.io/crimlab/reference/scests.md)
to read a Statistics Canada table once loaded.

## Examples

``` r
if (FALSE) { # \dontrun{
# needs the internet
indiex("mentalhealth")[, c("key", "rows")]
x <- indiex("lifeexp", load = TRUE)
} # }
try(indiex("weather"))
#> Error : indiex: `name` must be one of "health", "mentalhealth", "wellbeing", "suicide", "censuspop", "censushh", "lifeexp", "infantmort", not "weather".
```
