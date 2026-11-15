# Internal soundness of a region map

Checks that a region map (a crosswalk from units, such as correctional
facilities or police detachments, to the statistical regions that
contain them) lists each unit once, gives every unit a region, and names
only regions the reference geography knows. Run it on any crosswalk
before tables are built on it.

## Usage

``` r
rgnmapintg(map, unit, reg, regs = NULL)
```

## Arguments

- map:

  The region map, a data frame, one row per unit.

- unit, reg:

  The names of its unit and region columns.

- regs:

  `NULL`, or every valid region code (the identifier column of the
  population table); with `NULL` the third check is not run.

## Value

A data frame of class `rgnmapintg`, one row per check: `check`,
`observed`, `expected` (always 0) and `pass`.

## Details

These are the failures a comparison against published output cannot see,
because they would be present on both sides: a unit matched into two
regions, a unit matched into none, a region code that is a typo or
belongs to a neighbouring jurisdiction. None needs the geometry.

The checks as computed: "units assigned more than one region" counts
rows whose (non-blank) unit repeats an earlier row's, whether or not the
region differs, so an exact duplicate row also fails; "units assigned no
region" counts rows whose region is `NA` or blank; "region codes not in
the reference geography" counts non-blank region codes not in `regs`,
compared as text after trimming spaces (so `"3506"` and `3506` match,
`"03506"` does not).

## See also

[`rgnmapcmp()`](https://rootcoder007.github.io/crimlab/reference/rgnmapcmp.md),
[`rgnmapalt()`](https://rootcoder007.github.io/crimlab/reference/rgnmapalt.md).

## Examples

``` r
# a simulated crosswalk of facilities to census divisions
cw <- data.frame(inst = c("North Jail", "South Jail", "East Jail"),
                 cd = c("3557", "3520", "3506"))
rgnmapintg(cw, "inst", "cd", regs = c("3557", "3520", "3506"))
#> rgnmapintg: every check passed 
#>                                        check observed expected pass
#>          units assigned more than one region        0        0 TRUE
#>                     units assigned no region        0        0 TRUE
#>  region codes not in the reference geography        0        0 TRUE
cw$cd[3] <- "2406"
rgnmapintg(cw, "inst", "cd", regs = c("3557", "3520", "3506"))
#> rgnmapintg: 1 of 3 checks failed 
#>                                        check observed expected  pass
#>          units assigned more than one region        0        0  TRUE
#>                     units assigned no region        0        0  TRUE
#>  region codes not in the reference geography        1        0 FALSE
cw2 <- rbind(cw, data.frame(inst = c("North Jail", "West Jail"),
                            cd = c("3560", "")))
rgnmapintg(cw2, "inst", "cd")
#> rgnmapintg: 2 of 2 checks failed 
#>                                check observed expected  pass
#>  units assigned more than one region        1        0 FALSE
#>             units assigned no region        1        0 FALSE
try(rgnmapintg(cw, "institution", "cd"))
#> Error : rgnmapintg: `map` has no column `institution`.
```
