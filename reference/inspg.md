# The criminal, citizen and inspector game

Replicator dynamics of a rock-paper-scissors game of three strategies
(criminals, citizens, inspectors), each beating the next with stake `a`:
`x_i' = a x_i (x_{i+1} - x_{i+2})`. The centroid is the only interior
rest point and `x1 x2 x3` is conserved, so the shares circle the
centroid for ever and more inspection does not steadily lower crime
(D'Orsogna and Perc 2015). Integrated by fourth-order Runge-Kutta; the
drift of the invariant measures the integration error.

## Usage

``` r
inspg(a, x, steps = 1000, dt = 0.01)
```

## Arguments

- a:

  The stake: one finite number, not 0.

- x:

  The starting shares: three numbers from 0 summing to 1.

- steps:

  The Runge-Kutta steps: one whole number from 0, at most 1e7.

- dt:

  The step size: one positive number.

## Value

A list: `traj`, the shares after each step (`steps + 1` rows); `inv`,
the invariant `x1 x2 x3` along the path; `drift`, its largest departure
from the start; `dvdt`, its exact derivative at the start (0); `field`,
the vector field at the start; `fixed`, the centroid.

## References

D'Orsogna, M. R. and Perc, M. (2015). Statistical physics of crime: a
review. *Physics of Life Reviews*, 12, 1-21.
[doi:10.1016/j.plrev.2014.11.001](https://doi.org/10.1016/j.plrev.2014.11.001)

## Examples

``` r
g <- inspg(1, c(0.5, 0.3, 0.2), steps = 200)
g$drift
#> [1] 3.334832e-14
```
