# Hot spots and how long they last

Counts incidents in square grid cells of side `cell` (from the corner
`org`) for each period and marks as hot, in each period, the
`k = ceiling(top * m)` busiest of the `m` cells that ever hold an
incident, with every cell tied with the k-th; when a tie would make more
than `2 k` cells hot the counts do not single them out and the tied
cells are left out (counted in `tiedout`). Each period gets its
concentration (the fewest cells holding half its incidents), the share
of its incidents in hot cells, and the Jaccard overlap of its hot set
with the period before (`NA` when either is empty); each cell gets the
number of periods it was hot and its longest run of successive hot
periods, the persistence that tells chronic hot spots from flare-ups
(Weisburd, Bushway, Lum and Yang 2004).

## Usage

``` r
hsper(x, y = NULL, prd, cell, top = 0.05, org = NULL, lonlat = FALSE)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  names them, or `lon` and `lat` with `lonlat = TRUE`) and `y` `NULL`;
  or `x` an sf or sfc object of POINT geometries, whose coordinate
  reference system must agree with `lonlat`. An unnamed two-column
  matrix is read as `x` then `y`; other columns are never taken by
  position. Vectors given separately (`x`, `y`, and `w` or `t` where a
  function takes them) are matched by name when every one is named and
  no name repeats, as elsewhere in the package; otherwise by position.
  No missing or infinite values.

- prd:

  The period of each incident: labels (a year, a month, a quarter),
  taken in sorted order (a factor in the order of its levels), with no
  missing values.

- cell:

  The side of a cell: positive, in the unit of the coordinates (metres
  when `lonlat = TRUE`).

- top:

  The share of the active cells counted hot, above 0 and at most 1.

- org:

  `NULL` for the lowest `x` and `y`, or the grid's corner `c(x, y)`
  (projected, as the coordinates are, in metres when `lonlat = TRUE`).

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list of class `hsper`: `periods`, a data frame with one row per
period: `period`, `n` (incidents), `active` (cells with any), `nhalf`
and `shalf` (the fewest cells holding half the incidents, and that
number over the active cells of all periods), `hot`, `tiedin` (hot cells
beyond `k` from a tie), `tiedout`, `hotshr` (the share of incidents in
hot cells) and `jacc`; `cells`, a data frame with one row per active
cell: its centre `x`, `y` (in degrees when `lonlat = TRUE`), `n`, `nhot`
and `run`; `hot`, a logical matrix of cells by periods; `cell`, each
incident's cell (a row of `cells`); and `k`.

## References

Weisburd, D., Bushway, S., Lum, C. and Yang, S.-M. (2004). Trajectories
of crime at places: a longitudinal study of street segments in the city
of Seattle. *Criminology*, 42(2), 283-322.
[doi:10.1111/j.1745-9125.2004.tb00521.x](https://doi.org/10.1111/j.1745-9125.2004.tb00521.x)

## Examples

``` r
set.seed(1)
n <- 600
p <- cbind(x = c(rnorm(300, 2, 0.3), runif(300, 0, 10)),
           y = c(rnorm(300, 2, 0.3), runif(300, 0, 10)))
h <- hsper(p, prd = rep(2021:2024, length.out = n), cell = 1, top = 0.1)
h
#> hsper: 99 active cells, the 10 busiest hot in each of 4 periods
#>  period   n active nhalf      shalf hot tiedin tiedout    hotshr      jacc
#>    2021 150     57     4 0.04040404  20     10       0 0.7533333       NaN
#>    2022 150     57     4 0.04040404  19      9       0 0.7466667 0.2580645
#>    2023 150     54     4 0.04040404  11      1       0 0.6666667 0.2000000
#>    2024 150     62     4 0.04040404  17      7       0 0.7000000 0.2173913
subset(h$cells, nhot == 4)
#>           x        y  n nhot run
#> 11 1.513147 1.533045 88    4   4
#> 12 1.513147 2.533045 76    4   4
#> 21 2.513147 1.533045 90    4   4
#> 22 2.513147 2.533045 54    4   4
```
