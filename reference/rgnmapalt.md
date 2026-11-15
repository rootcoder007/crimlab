# Check a region map against an independently derived assignment

Compares the region each unit was assigned with the region a different
method assigns it (a name lookup, a postal geography, an administrative
list), and marks the disagreements that are already known and explained.
This is the only check here that can catch an error in the method that
built the map, because it does not use that method.

## Usage

``` r
rgnmapalt(map, unit, reg, route, known = character())
```

## Arguments

- map:

  The region map, a data frame.

- unit, reg:

  The names of its unit and region columns.

- route:

  The second method's assignment: a vector named by unit, or a data
  frame with the columns `unit` and `reg`.

- known:

  The units whose disagreement is expected and documented.

## Value

A data frame of the disagreements, unexplained first: `unit`, `primary`
(the map's region), `second` (the route's) and `known`. `sum(!x$known)`
unexplained disagreements is the result, zero the pass; `nrow(x)` should
equal the number of documented cases, as a documented case that stops
disagreeing means the route changed.

## Details

A second route usually has understood weaknesses: a place name that is a
community rather than a municipality, an amalgamated city, a township
absorbed into a neighbour. Listing them in `known` by name keeps the
check sharp instead of loosening a tolerance until everything passes,
and an unexpected disagreement cannot hide inside an allowance. Units
and regions are compared as text with spaces trimmed; units the route
does not cover, and route entries that are `NA` or blank, are skipped.

## See also

[`rgnmapcmp()`](https://rootcoder007.github.io/crimlab/reference/rgnmapcmp.md),
[`rgnmapintg()`](https://rootcoder007.github.io/crimlab/reference/rgnmapintg.md).

## Examples

``` r
cw <- data.frame(inst = c("North Jail", "South Jail", "Hill Jail"),
                 cd = c("3557", "3520", "3506"))
route <- c("North Jail" = "3557", "South Jail" = "3521",
           "Hill Jail" = "3519")
d <- rgnmapalt(cw, "inst", "cd", route, known = "Hill Jail")
d
#>         unit primary second known
#> 1 South Jail    3520   3521 FALSE
#> 2  Hill Jail    3506   3519  TRUE
sum(!d$known)
#> [1] 1
try(rgnmapalt(cw, "inst", "cd", c("3557", "3520")))
#> Error : rgnmapalt: `route` must be named by unit, or be a data frame.
```
