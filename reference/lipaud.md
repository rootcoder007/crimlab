# Individual fairness: a Lipschitz audit

Dwork et al. (2012) ask that similar individuals be treated similarly:
that the gap between two people's scores be at most `lip` times the
distance between them, for a distance that says who is alike for the
task. `lipaud()` checks the given pairs against that rule and reports
the smallest constant the scores satisfy on them.

## Usage

``` r
lipaud(h, d = NULL, pairs = NULL, lip = 1)
```

## Arguments

- h:

  The scores: finite numbers, one per person.

- d:

  The distances: with `pairs` `NULL`, a
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html) object or a
  symmetric matrix between all the people (every pair is audited); with
  `pairs`, `NULL` (all 1) or one number from 0 per pair.

- pairs:

  `NULL`, or the pairs to audit: a matrix of two columns of person
  numbers from 1.

- lip:

  The constant: one number from 0.

## Value

A list: `lreq` (the smallest constant the scores satisfy on these
pairs), `nviol` and `rate` (the pairs breaking the rule at `lip`, and
their share), `ok` (`TRUE` when none does), `maxgap` (the largest gap),
`worst` (the first pair with the largest ratio), and `pairs`, a data
frame with one row per pair: `i`, `j`, `d`, `gap`, `ratio` and `viol`.

## Details

A pair at distance 0 is two people the distance cannot tell apart: any
gap between their scores breaks the rule, and no constant repairs it
(its ratio is `Inf`). The audit is only as good as the distance: it says
whether the scores treat alike those the distance calls alike, not
whether the distance is the right one.

## References

Dwork, C., Hardt, M., Pitassi, T., Reingold, O. and Zemel, R. (2012).
Fairness through awareness. *Proceedings of the 3rd Innovations in
Theoretical Computer Science Conference*, 214-226.
[doi:10.1145/2090236.2090255](https://doi.org/10.1145/2090236.2090255)

## Examples

``` r
x <- c(1, 2, 2.5, 6)
h <- c(0.10, 0.20, 0.60, 0.65)
# scores on a scale of 0 to 1 against distances in years
lipaud(h, dist(x), lip = 0.2)
#> $lreq
#> [1] 0.8
#> 
#> $nviol
#> [1] 2
#> 
#> $rate
#> [1] 0.3333333
#> 
#> $ok
#> [1] FALSE
#> 
#> $maxgap
#> [1] 0.55
#> 
#> $worst
#> [1] 2 3
#> 
#> $pairs
#>   i j   d  gap      ratio  viol
#> 1 1 2 1.0 0.10 0.10000000 FALSE
#> 2 1 3 1.5 0.50 0.33333333  TRUE
#> 3 1 4 5.0 0.55 0.11000000 FALSE
#> 4 2 3 0.5 0.40 0.80000000  TRUE
#> 5 2 4 4.0 0.45 0.11250000 FALSE
#> 6 3 4 3.5 0.05 0.01428571 FALSE
#> 
```
