# A true base rate behind a noisily recorded one

A recorded outcome is often a proxy: a rearrest for a reoffence, a
recorded crime for one committed. With false-positive noise `a`, the
chance a case without the outcome is recorded with it, and
false-negative noise `b`, the chance a case with it is not, a true base
rate `p` is recorded as `pobs = p (1 - b) + (1 - p) a`. Given only the
boxes the noise rates lie in, `nsbnd()` gives the interval of true rates
consistent with each recorded rate, and says for each pair of groups
whether the order of their true rates follows from the records whatever
the noise.

## Usage

``` r
nsbnd(pobs, a, b)
```

## Arguments

- pobs:

  Recorded rates: numbers from 0 to 1, named or not.

- a:

  The box of the false-positive noise rate: one number `a` (the box from
  0 to `a`) or two, its ends; from 0.

- b:

  The box of the false-negative noise rate, in the same way. The upper
  ends of `a` and `b` must add to less than 1.

## Value

A list: `bnd`, a data frame with one row per rate: `rate` (its name, or
its position), `pobs`, `lo`, `hi` and `width`; `pairs`, one row per pair
`i` before `j`: `i`, `j`, `ord` (`"lower"`: `i`'s true rate is below
`j`'s for every admissible noise; `"higher"`; `"undecided"`; `NA` when
either has no interval) and `bd` (`NA` when the rates are equal or no
box keeps them ordered).

## Details

The true rate `(pobs - a) / (1 - a - b)` falls as `a` grows and rises as
`b` grows, so its ends are reached at corners of the boxes: the lower at
the largest `a` (capped at `pobs`) and smallest `b`, the upper at the
smallest `a` and largest `b` (capped at `1 - pobs`). Both ends are
attained, so no narrower interval follows from the boxes alone; with
known noise rates (each box one point) the interval is the one true
rate. A recorded rate that no admissible noise can produce (below the
smallest `a`, or above `1 - ` the smallest `b`) has no interval (`NaN`).

Two groups' true rates are ordered whatever the noise exactly when their
intervals do not overlap; when they overlap there is admissible noise
under which either is the larger. `bd` is the largest upper end of the
`b` box, `a` held, at which the lower recorded rate's interval stays
below the other's: how much under-recording the order survives. The
boxes are the analyst's claim about how the records were made; the
bounds are only as good as that claim. See `proof("nsbnd")`.

## References

Fogliato, R., Chouldechova, A. and G'Sell, M. (2020). Fairness
evaluation in presence of biased noisy labels. *Proceedings of the 23rd
International Conference on Artificial Intelligence and Statistics*,
PMLR 108, 2325-2336.
[doi:10.48550/arXiv.2003.13808](https://doi.org/10.48550/arXiv.2003.13808)

## Examples

``` r
# rearrest rates of two groups; up to 10% false positives, 20% missed
nsbnd(c(a = 0.35, b = 0.55), a = 0.10, b = 0.20)
#> $bnd
#>   rate pobs        lo     hi     width
#> 1    a 0.35 0.2777778 0.4375 0.1597222
#> 2    b 0.55 0.5000000 0.6875 0.1875000
#> 
#> $pairs
#>   i j   ord  bd
#> 1 a b lower 0.3
#> 
# with up to 40% missed, the order is no longer identified
nsbnd(c(a = 0.35, b = 0.55), a = 0.05, b = 0.40)$pairs
#>   i j       ord    bd
#> 1 a b undecided 0.335
# known noise rates give the true rate itself
nsbnd(0.4 * 0.8 + 0.6 * 0.1, a = c(0.1, 0.1), b = c(0.2, 0.2))$bnd
#>   rate pobs  lo  hi width
#> 1    1 0.38 0.4 0.4     0
```
