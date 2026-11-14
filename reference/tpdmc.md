# Decay chain of a body

A discrete-time Markov chain of decomposition: from each transient stage
(fresh, bloat, active and advanced decay, or the user's own stages) a
body advances to the next stage with probability `decay * (1 - pres)`
(after the last, to skeletal), mummifies with probability `mum * pres`,
and otherwise stays. `pres` is the strength of what preserves the body
(aridity, quicklime, a sealed casket), 0 for none. From the fundamental
matrix come the probability of ending skeletal or mummified from each
stage and the expected steps to get there. With no preservation nothing
mummifies, so the column `mummified` of `absorb` is also the change in
the chance of mummifying that the preservation makes. With `steps` above
0 one path is drawn from `start`.

## Usage

``` r
tpdmc(
  pres = 0,
  decay = 0.5,
  mum = 0.5,
  states = c("fresh", "bloat", "active", "advanced"),
  start = 1,
  steps = 0,
  seed = NULL
)
```

## Arguments

- pres:

  The preservation: one number from 0 to 1.

- decay:

  The chance to advance a stage without preservation: one number above
  0, at most 1.

- mum:

  The chance to mummify at full preservation: one number from 0 to 1.

- states:

  The names of the transient stages, in order.

- start:

  The stage a drawn path starts from: a name of `states` or its number.

- steps:

  The longest path drawn (it stops when absorbed): one whole number from
  0, 0 for none.

- seed:

  `NULL` (the draws follow
  [`set.seed()`](https://rdrr.io/r/base/Random.html)) or a whole number
  from 0.

## Value

A list: `P`, the transition matrix (stages, then `skeletal` and
`mummified`); `N`, the fundamental matrix of the transient stages;
`absorb`, the probability of ending `skeletal` or `mummified` from each
stage; `time`, the expected steps to absorption from each stage;
`start`; `path`, `NULL` or the states of the drawn path.

## References

Kemeny, J. G. and Snell, J. L. (1960). *Finite Markov Chains*. Van
Nostrand.

## Examples

``` r
tpdmc(pres = 0.4)$absorb
#>          skeletal mummified
#> fresh      0.1296    0.8704
#> bloat      0.2160    0.7840
#> active     0.3600    0.6400
#> advanced   0.6000    0.4000
tpdmc(pres = 0.4, steps = 50, seed = 1)$path
#> [1] "fresh"     "fresh"     "mummified"
```
