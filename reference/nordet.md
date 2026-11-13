# Noisy-OR detection by patrols

The chance that an incident is seen by patrols when each officer within
`radius` of it sees it independently with probability `p`: with `k`
officers in range, `1 - (1 - p)^k` (the detection model of Barman and
Barman 2026). With `draw = TRUE`, one detection outcome is drawn for
each incident. A model of who is seen, under stated assumptions; how
patrols are actually placed and how much they see are matters of fact
the model does not settle.

## Usage

``` r
nordet(x, off, radius, p = 0.85, draw = FALSE, seed = NULL)
```

## Arguments

- x:

  The incidents: a numeric matrix of two columns, or a data frame with
  columns `x` and `y` (or two columns), planar coordinates.

- off:

  The officers, in the same forms; `NULL` or no rows for none.

- radius:

  The range: one positive number, in the coordinates' units (distances
  are Euclidean; an officer exactly at `radius` counts).

- p:

  The chance one officer in range sees an incident: one number above 0
  and at most 1.

- draw:

  `TRUE` to draw whether each incident is detected.

- seed:

  `NULL` to draw the seed of the detections from R's generator, or one
  whole number from 0; read only with `draw`.

## Value

A data frame with one row per incident: `k` (officers in range), `prob`
and, with `draw`, `det` (1 detected, 0 not).

## Examples

``` r
set.seed(1)
inc <- matrix(runif(20), ncol = 2)
off <- matrix(runif(10), ncol = 2)
nordet(inc, off, radius = 0.3, p = 0.7)
#>    k prob
#> 1  2 0.91
#> 2  2 0.91
#> 3  0 0.00
#> 4  2 0.91
#> 5  1 0.70
#> 6  2 0.91
#> 7  0 0.00
#> 8  0 0.00
#> 9  1 0.70
#> 10 1 0.70
```
