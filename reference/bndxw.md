# Crosswalk between boundary versions

Carries values kept on one set of areas (the 158 neighbourhoods of a
city's new boundaries, a census vintage's tracts) to another (the 140
old neighbourhoods) through a correspondence table: each link says what
share `w` of a source area's value falls in a target area. Counts of a
source wholly inside one target move whole (weight 1, exact); a source
split between targets is shared by area or population (uniform density
within the source, an assumption). The same table read the other way
disaggregates. `map` also gives each source's main target, the column to
add to a table keyed by the source areas.

## Usage

``` r
bndxw(from, to, x, unit, w = NULL)
```

## Arguments

- from, to:

  The correspondence: the source and target area of each link, labels. A
  `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  is taken too (as `from`): its elements fill the arguments by name.

- x:

  The values to carry: finite numbers, one per record.

- unit:

  The source area of each record: labels (as in `from`).

- w:

  `NULL` (every link weight 1), or the share of each link: numbers from
  0.

## Value

A list: `tab`, a data frame with one row per target (`to`, `value`);
`map`, one row per source (`from`; `to`, the target of its largest
share, the first of ties; `share`, that share; `wsum`, the sum of its
shares, 1 when nothing is lost or counted twice); `lost`, the total of
the records whose area is not in `from`; `unmatched`, those areas.

## Examples

``` r
# old areas a, b merge into A; c splits between A and B by area
bndxw(from = c("a", "b", "c", "c"), to = c("A", "A", "A", "B"),
      x = c(3, 4, 5, 8), unit = c("a", "b", "b", "c"),
      w = c(1, 1, 0.25, 0.75))
#> $tab
#>   to value
#> 1  A    14
#> 2  B     6
#> 
#> $map
#>   from to share wsum
#> 1    a  A  1.00    1
#> 2    b  A  1.00    1
#> 3    c  B  0.75    1
#> 
#> $lost
#> [1] 0
#> 
#> $unmatched
#> character(0)
#> 
```
