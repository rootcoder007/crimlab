# Square-root law of response distance

Kolesar and Blum's square-root law: with available units spread like a
random (Poisson) pattern over a region of area `area`, the mean distance
from a call to the nearest available unit is \\c \sqrt{A / (u (1 -
b))}\\ for `units` u of which a share `busy` b is busy; \\c = 1/2\\ in
straight lines and \\\sqrt{2 \pi} / 4 \approx 0.627\\ on a rectangular
street grid. Quadrupling the units halves the distance and scaling area
and units together leaves it unchanged (`proof("sqrtlw")`). The busy
share can come from
[`hypq()`](https://rootcoder007.github.io/crimlab/reference/hypq.md)'s
workloads.

## Usage

``` r
sqrtlw(area, units, busy = 0, metric = "euclid", speed = NULL)
```

## Arguments

- area:

  The area served: one positive number (squared distance units).

- units:

  The units deployed: one positive number.

- busy:

  The share of units busy: one number from 0, below 1.

- metric:

  `"euclid"` (straight lines) or `"grid"` (right-angle streets).

- speed:

  `NULL` or the travel speed (distance per unit time): one positive
  number.

## Value

A list: `dist`, the mean distance; `const`, the constant `c`; `avail`,
the units available; `time`, the mean travel time (`NA` without
`speed`).

## References

Kolesar, P. and Blum, E. H. (1973). Square root laws for fire engine
response distances. *Management Science*, 19(12), 1368-1378.
[doi:10.1287/mnsc.19.12.1368](https://doi.org/10.1287/mnsc.19.12.1368)

## Examples

``` r
sqrtlw(100, 4)
#> $dist
#> [1] 2.5
#> 
#> $const
#> [1] 0.5
#> 
#> $avail
#> [1] 4
#> 
#> $time
#> [1] NA
#> 
sqrtlw(100, 8, busy = 0.3, metric = "grid", speed = 0.5)
#> $dist
#> [1] 2.648109
#> 
#> $const
#> [1] 0.6266571
#> 
#> $avail
#> [1] 5.6
#> 
#> $time
#> [1] 5.296219
#> 
```
