# Hypercube queue of patrol units

Larson's hypercube model of patrol cars: calls arrive from each atom
(beat, reporting area) at rate `lam` and go to the first free unit in
the atom's order of preference; each unit serves at rate `mu`. The state
is the set of busy units, \\2^n\\ states for `n` units, and the
stationary distribution is found exactly. From it come each unit's
workload (the share of time busy), the share of calls each unit answers
from each atom, the probability that every unit is busy (a lost call, or
with `queue = "wait"` a call that waits; Larson's infinite-line model,
the queue lumped into the all-busy state), the mean queue and wait, and
with travel times the mean response travel time. Whatever the orders,
the number of busy units follows Erlang's loss or delay system. The
noisy-OR chance that patrols detect an event is
[`nordet()`](https://rootcoder007.github.io/crimlab/reference/nordet.md).

## Usage

``` r
hypq(lam, pref, mu = 1, queue = "loss", tt = NULL)
```

## Arguments

- lam:

  The call rate of each atom: positive numbers.

- pref:

  The dispatch orders: a matrix with one row per atom, each a
  permutation of the unit numbers 1 to `n` (at most 10 units), or a list
  of such vectors.

- mu:

  The service rate of a unit (1 / mean service time): one positive
  number.

- queue:

  `"loss"` (a call finding every unit busy is lost) or `"wait"` (it
  waits; the calls must be fewer than the units serve).

- tt:

  `NULL` or the travel times: a matrix with one row per unit and one
  column per atom.

## Value

An object of class `hypq`, a list: `state`, each set of busy units
(`busy`, one character per unit, `"1"` busy) with its probability `p`;
`wl`, the workloads; `disp`, the share of all calls answered at once by
each unit (rows) from each atom (columns); `queue`; `pfull`, the
probability that every unit is busy; `lq` and `wq`, the mean queue and
wait (0 for `"loss"`); `rho`, the calls over what the units can serve;
`travel`, `NULL` or the mean travel time overall (`mean`) and by atom
(`atom`).

## References

Larson, R. C. (1974). A hypercube queuing model for facility location
and redistricting in urban emergency services. *Computers and Operations
Research*, 1(1), 67-95.
[doi:10.1016/0305-0548(74)90076-8](https://doi.org/10.1016/0305-0548%2874%2990076-8)

## Examples

``` r
pref <- rbind(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1), c(3, 1, 2))
tt <- matrix(c(2, 5, 8, 6, 3, 7, 9, 6, 2, 7, 8, 3), 3)
hypq(c(0.4, 0.7, 0.3, 0.6), pref, mu = 1.5, tt = tt)
#> hypq: 3 units, rho 0.4444444; P(all busy) 0.109215
#> workloads: 0.3722186 0.4025169 0.4129779
#> mean travel time 3.944009
```
