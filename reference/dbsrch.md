# Database search

What a match found by trawling a database is worth. `np` is the bound of
the second National Research Council report: the match probability times
the number of profiles searched. `any` is the exact chance of at least
one chance match among that many unrelated profiles, never above `np`.
`post` takes the other view (Balding and Donnelly): everyone else in the
database was excluded, so among `pop` possible sources with equal prior
the matching person's probability of being the source goes up, not down.
See `proof("dbsrch")`.

## Usage

``` r
dbsrch(mp, n, pop = NULL)
```

## Arguments

- mp:

  The profile match probability: one number above 0 and at most 1.

- n:

  The number of profiles searched, the match included.

- pop:

  `NULL`, or the number of possible sources, at least `n`.

## Value

A named numeric vector: `np`, `any` and, when `pop` is given, `post`.

## References

Balding, D. J. and Donnelly, P. (1996). Evaluating DNA profile evidence
when the suspect is identified through a database search. *Journal of
Forensic Sciences*, 41(4), 603-607.
[doi:10.1520/JFS13961J](https://doi.org/10.1520/JFS13961J)

## Examples

``` r
dbsrch(mp = 1e-9, n = 5e6, pop = 3e7)
#>          np         any        post 
#> 0.005000000 0.004987521 0.975609756 
```
