# Transition matrix between states

The moves of people between states (regions, institutions, security
levels) from their sequence of records: the counts and probabilities of
each move between consecutive records of the same person, the stationary
distribution of the chain (where the population would settle if the
probabilities held), the share of moves that stay in their state, and
the Theil index of the moves between states (0 when they spread evenly
over the pairs of states, higher when a few pairs carry them).

## Usage

``` r
trnmat(state, id, time = NULL)
```

## Arguments

- state:

  The state of each record: labels.

- id:

  The person (or sequence) of each record: labels. Give
  `paste(id, year)` to keep the moves within a year.

- time:

  `NULL` (the order of the rows) or the time of each record: numbers;
  ties keep the order of the rows.

## Value

An object of class `trnmat`, a list: `count` and `prob`, the matrices of
moves from (rows) and to (columns); `pi`, the stationary distribution
(`NA` when a state is never left or it is not unique); and `summ`
(`moves`, `stay`, `theil`).

## References

Kemeny, J. G. and Snell, J. L. (1960). *Finite Markov Chains*. Van
Nostrand.

Theil, H. (1967). *Economics and Information Theory*. North-Holland.

## Examples

``` r
r <- trnmat(c("N", "N", "E", "N", "W", "W", "E", "E"),
            c(1, 1, 1, 2, 2, 2, 3, 3), time = c(1, 2, 3, 1, 2, 3, 1, 2))
r$prob
#>           E         N         W
#> E 1.0000000 0.0000000 0.0000000
#> N 0.3333333 0.3333333 0.3333333
#> W 0.0000000 0.0000000 1.0000000
r$pi
#>  E  N  W 
#> NA NA NA 
```
