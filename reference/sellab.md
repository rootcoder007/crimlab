# Evaluating a release rule from the released alone

Failures (a missed court date, a new arrest) are seen only for the
defendants a judge released. A proposed rule that releases only people
the judge released is scored exactly from their outcomes; one that
releases people the judge detained has a failure rate known only up to
the worst and best cases for them, an interval whose width is the share
of the rule's releases the judge detained. Both ends can occur. The
naive score on the rule's observed releases alone always lies inside the
interval, but reports one point of it as if it were known. See
`proof("sellab")`.

## Usage

``` r
sellab(y, rel, rule, w = NULL)
```

## Arguments

- y:

  Failures: 1 or 0 for each released defendant; values for detained
  defendants are ignored and may be `NA`.

- rel:

  Released by the judge: `TRUE` or `FALSE` for each defendant. When `y`,
  `rel`, `rule` (and `w`) are all named with unique names, they are
  matched by name; otherwise by position.

- rule:

  Released by the rule: `TRUE` or `FALSE`, same length.

- w:

  `NULL` for equal weights, or one non-negative weight per defendant.

## Value

A named numeric vector: `obs` (the failure rate among the judge's
releases), `unobs` (the share of the rule's releases the judge
detained), `lower` and `upper` (the bounds on the rule's failure rate),
`width` (`upper - lower`, equal to `unobs`) and `naive` (the failure
rate on the rule's releases that were observed; `NaN` when none were).

## References

Lakkaraju, H., Kleinberg, J., Leskovec, J., Ludwig, J. and Mullainathan,
S. (2017). The selective labels problem. *Proceedings of the 23rd ACM
SIGKDD International Conference on Knowledge Discovery and Data Mining*,
275-284.
[doi:10.1145/3097983.3098066](https://doi.org/10.1145/3097983.3098066)

## Examples

``` r
set.seed(4)
risk <- runif(500)
y <- rbinom(500, 1, risk)
sellab(y, rel = risk < 0.7, rule = risk < 0.8)
#>       obs     unobs     lower     upper     width     naive 
#> 0.3488372 0.1202046 0.3069054 0.4271100 0.1202046 0.3488372 
```
