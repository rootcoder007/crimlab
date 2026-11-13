# Points in polygons, and the population the units cover

Finds the region holding each point (a facility, a station, an incident)
by the even-odd rule, a point on a boundary counting as inside; a region
given as several rings (parts and holes) is one region, and a hole is
any ring inside another of the same region. A point outside every region
gets `NA`; one inside overlapping regions takes the first by label and
`nreg` says how many hold it. Coordinates are taken as planar (longitude
and latitude work for regions small enough that their edges are straight
in degrees).

## Usage

``` r
pinpl(x, y = NULL, poly, pop = NULL)
```

## Arguments

- x, y:

  The points in planar coordinates (project longitude and latitude with
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  first): `x` and `y` numeric vectors of equal length, or `x` a matrix
  or data frame with columns named `x` and `y`. At least two points; no
  missing values.

- poly:

  The regions: a data frame with columns `x` and `y` (or `lon` and
  `lat`), the vertices in order around each ring, `id` (the region's
  label) and, for regions of several rings, `ring`.

- pop:

  `NULL`, or the population of each region (from 0), in the sorted order
  of the labels.

## Value

A list of class `pinpl`: `reg` (each point's region label, `NA`
outside), `nreg` (the regions holding it), `regs` (a data frame of
`region`, `units`, `pop`, `hasunit` and `popshr`), and `cov`, a named
vector: `regions`, `withunit`, `units`, `totpop`, `covpop`, `uncpop` and
`covshr` (`NA` without populations).

## Details

With the regions' populations, the coverage: the population of the
regions holding a unit, and its share of the total. That share is not a
rate denominator. A region without a unit is not unserved: units serve
catchments, which a location does not state, and a rate built on the
regions that hold units pairs a partial denominator with a whole
numerator, inflating it unevenly.

Polygons from a spatial file convert with, for instance,
[`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html),
whose columns `X`, `Y`, `L1` (ring) and `L2` or `L3` (feature) become
`x`, `y`, `ring` and `id`.

## References

Shimrat, M. (1962). Algorithm 112: position of point relative to
polygon. *Communications of the ACM*, 5(8), 434.
[doi:10.1145/368637.368653](https://doi.org/10.1145/368637.368653)

## Examples

``` r
sq <- function(x0, id) {
  data.frame(x = x0 + c(0, 1, 1, 0), y = c(0, 0, 1, 1), id = id)
}
poly <- rbind(sq(0, "A"), sq(1, "B"), sq(2, "C"))
u <- pinpl(cbind(x = c(0.5, 0.2, 2.5, 9), y = c(0.5, 0.7, 0.5, 9)),
           poly = poly, pop = c(1200, 800, 450))
u
#> pinpl: 3 of 4 points in 2 of 3 regions; they cover 67.3% of the population, which is not a rate denominator
#>  region units  pop hasunit    popshr
#>       A     2 1200    TRUE 0.4897959
#>       B     0  800   FALSE 0.3265306
#>       C     1  450    TRUE 0.1836735
u$reg
#> [1] "A" "A" "C" NA 
```
