# Presets: a published dataset's columns onto a generic function's arguments

A preset says how the columns of one published dataset (an Ontario
segregation table, the Special Investigations Unit's cases, a
use-of-force file) feed the arguments of one of crimlab's generic
functions, so that `generic(clprep(data, preset))` gives the published
analysis while the same generic runs on anyone's data given as plain
vectors. A preset describes a layout and ships no data: the user
supplies the table.

## Usage

``` r
clprep(data = NULL, preset = NULL)
```

## Arguments

- data:

  The dataset: a data frame with the columns the preset names. `NULL`
  (with `preset` `NULL`) lists the presets.

- preset:

  The preset's name, as `clprep()` lists them; `NULL` for the list.

## Value

With a preset, an object of class `clpre`: a list of the generic's
arguments, which the generic takes as its first argument, with
attributes `preset`, `generic`, `source`, `units`, `note` and `drop`.
Without one, a data frame of the presets: `preset`, `generic`, `args`
and `source`.

## Details

The presets are rows of `presets.tsv` (in the installed package's
`presets` folder), one per preset and argument: `preset`, `generic`,
`arg`, `column` (the dataset's column), `recode`, `filter`, `units`,
`note`, `source`. Python reads the same file. A `recode` is one of:
empty (the column as it is); `fac` (a factor, one effect per value);
`num` (numbers); `one` (a 1 for every row); `sub:a,b` (characters a to
b, as a year from an ISO date); `const:v` (the argument is v, no
column); or a mapping `x|y=1;has:abc=1;*=0` from labels, compared
without case or surrounding space (`has:` a label containing the text),
to values, `*` covering every other label and `*=NA` dropping those
rows. A `filter` `col=value` keeps the rows whose column has that value.
An `arg` of `recipe` derives columns first: `long` stacks the two count
columns named in `recode` as `count` with `with` 1 for the first; `acvm`
and `pyr` reduce one row per placement to one per person and year (see
Details of the returned attributes).

Rows with a missing value in any column the preset reads are dropped, as
the source analyses drop them; their number is the attribute `drop`.

## Examples

``` r
head(clprep())
#>         preset generic               args
#> 1 otis-b03-cnt  cntreg y grp adj clus gee
#> 2 otis-b04-cnt  cntreg          y grp adj
#> 3 otis-b06-cnt  cntreg y grp adj clus gee
#> 4 otis-b07-cnt  cntreg   recipe y grp adj
#> 5 otis-b08-cnt  cntreg y grp adj clus gee
#> 6 otis-b09-cnt  cntreg          y grp adj
#>                                                                                                                        source
#> 1 Ontario Ministry of the Solicitor General, segregation and restrictive confinement data (Ontario Data Catalogue), table b03
#> 2 Ontario Ministry of the Solicitor General, segregation and restrictive confinement data (Ontario Data Catalogue), table b04
#> 3 Ontario Ministry of the Solicitor General, segregation and restrictive confinement data (Ontario Data Catalogue), table b06
#> 4 Ontario Ministry of the Solicitor General, segregation and restrictive confinement data (Ontario Data Catalogue), table b07
#> 5 Ontario Ministry of the Solicitor General, segregation and restrictive confinement data (Ontario Data Catalogue), table b08
#> 6 Ontario Ministry of the Solicitor General, segregation and restrictive confinement data (Ontario Data Catalogue), table b09
# a table laid out as the published OTIS c01 table (numbers made up)
d <- data.frame(EndFiscalYear = rep(2020:2022, each = 2),
                Gender = rep(c("Female", "Male"), 3),
                NumberIndividuals_RestrictiveConfinement =
                  c(40, 310, 52, 335, 47, 360))
p <- clprep(d, "otis-c01-cnt")
cntreg(p)
#> cntreg: 6 counts; reference 0
#> rate ratios (95% intervals):
#>   fit term       irr        lo        hi             p   sep
#>  pois    1 0.1383085 0.1158296 0.1651497 6.153279e-106 FALSE
#>    nb    1 0.1383085 0.1158295 0.1651499 6.174131e-106 FALSE
#> negative binomial theta 7655874; AIC Poisson 48.8636, negative binomial 50.86374
# the same model on plain vectors
cntreg(d$NumberIndividuals_RestrictiveConfinement, d$Gender == "Female",
       adj = factor(d$EndFiscalYear))
#> cntreg: 6 counts; reference FALSE
#> rate ratios (95% intervals):
#>   fit term       irr        lo        hi             p   sep
#>  pois TRUE 0.1383085 0.1158296 0.1651497 6.153279e-106 FALSE
#>    nb TRUE 0.1383085 0.1158295 0.1651499 6.174131e-106 FALSE
#> negative binomial theta 7655874; AIC Poisson 48.8636, negative binomial 50.86374
```
