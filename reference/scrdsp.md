# A risk score's spread across groups

Each group's number of cases, mean, median and standard deviation of a
score, each mean's gap to a reference group with its interval, and the
one-way analysis of variance of the score on group, with Welch's (1951)
test, which does not assume equal variances.

## Usage

``` r
scrdsp(score, group, ref = NULL, conf = 0.95)
```

## Arguments

- score:

  The scores: finite numbers, one per case.

- group:

  The group of each case: labels.

- ref:

  The reference group: `NULL` for the group with the lowest mean (the
  first on ties), its label, or its position from 1.

- conf:

  The confidence level of the intervals, strictly between 0 and 1.

## Value

A list: `tab`, a data frame with one row per group: `group`, `n`,
`mean`, `median`, `sd` (`NaN` for one case), `gap`, `se`, `lo` and `hi`
(`NaN` for the reference); `anova` and `welch`, named vectors `f`,
`df1`, `df2` and `p`; `spread` (the largest mean minus the smallest);
`ref`.

## Details

The analysis of variance is the standard one
([`stats::anova()`](https://rdrr.io/r/stats/anova.html) of
[`stats::lm()`](https://rdrr.io/r/stats/lm.html)): a group with one case
adds to the between-group sum of squares and nothing to the within-group
one. The gaps' intervals use the pooled within-group variance on its
degrees of freedom. Welch's test needs two or more cases and some spread
in every group, and is `NaN` otherwise.

A gap in mean scores is not by itself a sign of bias: a score built to
predict an outcome differs between groups whose rates of that outcome
differ. Read it with
[`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md),
which compares the scores with what was recorded, or
[`fair()`](https://rootcoder007.github.io/crimlab/reference/fair.md),
which compares them with outcomes case by case.

## References

Welch, B. L. (1951). On the comparison of several mean values: an
alternative approach. *Biometrika*, 38(3-4), 330-336.
[doi:10.1093/biomet/38.3-4.330](https://doi.org/10.1093/biomet/38.3-4.330)

## Examples

``` r
set.seed(7)
g <- rep(c("x", "y", "z"), c(40, 40, 20))
s <- rnorm(100, ifelse(g == "x", 0, ifelse(g == "y", 0.5, 1)))
scrdsp(s, g)
#> $tab
#>   group  n      mean    median        sd       gap        se          lo
#> 1     x 40 0.2708131 0.1684252 1.0461883 0.0000000       NaN         NaN
#> 2     y 40 0.6225100 0.6625734 0.9275965 0.3516969 0.2144616 -0.07395004
#> 3     z 20 0.9068368 0.7723597 0.8266946 0.6360236 0.2626607  0.11471478
#>          hi
#> 1       NaN
#> 2 0.7773438
#> 3 1.1573325
#> 
#> $anova
#>           f         df1         df2           p 
#>  3.18629701  2.00000000 97.00000000  0.04568425 
#> 
#> $welch
#>           f         df1         df2           p 
#>  3.31086371  2.00000000 54.14739082  0.04400091 
#> 
#> $spread
#> [1] 0.6360236
#> 
#> $ref
#> [1] "x"
#> 
```
