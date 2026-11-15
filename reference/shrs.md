# Share of a total, in percent

Expresses each count as a percentage of a total (for example each
division's share of all stops, or each offence's share of all
admissions), with a Wilson confidence interval. A share's denominator is
the total of the same events, so shares over a complete grouping sum to
100. Use a rate (such as
[`grprr()`](https://rootcoder007.github.io/crimlab/reference/grprr.md))
instead when the denominator is a population; "40% of stops" and "40
stops per 1,000 residents" are different claims.

## Usage

``` r
shrs(x, ...)

# S3 method for class 'data.frame'
shrs(x, count = NULL, by = NULL, total = NULL, conflvl = 0.95, ...)

# Default S3 method
shrs(x, total = NULL, conflvl = 0.95, ...)
```

## Arguments

- x:

  A data frame, or a vector of counts.

- ...:

  Unused; anything given is an error.

- count:

  The column of counts (required for a data frame), unquoted or as a
  string: non-negative numbers.

- by:

  `NULL`, or the grouping columns (names).

- total:

  The denominator: `NULL` (the sum of the counts, so the shares sum to
  100), one non-negative number, or (for a data frame) the name of a
  column whose sum is the total, for when some of it is not in the
  table.

- conflvl:

  The level of the interval, in (0, 1).

## Value

A data frame of class `shrs`: the `by` columns (or `label`), `count`,
`total`, `share`, `lower` and `upper` (percent), with an attribute
`share`, a list of `count`, `by`, `total` and `conflvl`.

## Details

For count \\y\\ out of total \\n\\, the share is \\100\\y/n\\ and the
interval is Wilson's score interval for \\\hat p = y/n\\, \$\$\frac{\hat
p + \frac{z^2}{2n} \pm z\sqrt{\frac{\hat p(1-\hat p)}{n} +
\frac{z^2}{4n^2}}}{1 + z^2/n},\$\$ with \\z\\ the standard normal
quantile for `conflvl`, in percent and clipped to \\\[0, 100\]\\
(exactly 0 at a count of 0 and 100 at the whole total). The normal
interval on a proportion is wrong in exactly the cases people reach for
it (small counts, shares near 0 or 100), where it runs past the ends of
the scale; Wilson's stays inside and is accurate there (Wilson, 1927).
It treats each count as a binomial draw from the total, an approximation
when the total is the sum of the counts.

With a data frame, counts are summed within each combination of the `by`
columns (missing values ignored), the groups sorted as
[`stats::aggregate()`](https://rdrr.io/r/stats/aggregate.html) sorts
them. The default method treats each element of `x` as a row; a named
`x` gives a leading `label` column.

Errors: a count column that is not numeric or not found; negative
counts; a total of zero or below; a count larger than the total;
`conflvl` outside (0, 1).

## References

Wilson EB (1927). Probable inference, the law of succession, and
statistical inference. *Journal of the American Statistical Association*
22(158), 209–212.
[doi:10.1080/01621459.1927.10502953](https://doi.org/10.1080/01621459.1927.10502953)

## Examples

``` r
# simulated stops by division
stops <- data.frame(division = c("North", "South", "East"),
                    stops = c(412, 77, 3))
shrs(stops, stops, by = "division")
#> shrs: shares of 492 in percent, 95% Wilson intervals
#>  division count total      share      lower     upper
#>      East     3   492  0.6097561  0.2075843  1.777215
#>     North   412   492 83.7398374 80.2200225 86.736863
#>     South    77   492 15.6504065 12.7074531 19.125597
shrs(stops, "stops", by = "division", total = 10000)
#> shrs: shares of 10000 in percent, 95% Wilson intervals
#>  division count total share      lower      upper
#>      East     3 10000  0.03 0.01020322 0.08817358
#>     North   412 10000  4.12 3.74774662 4.52748907
#>     South    77 10000  0.77 0.61657440 0.96123408
shrs(c(murder = 3, robbery = 7))
#> shrs: shares of 10 in percent, 95% Wilson intervals
#>    label count total share    lower    upper
#>   murder     3    10    30 10.77913 60.32219
#>  robbery     7    10    70 39.67781 89.22087
try(shrs(c(5, 2), total = 4))
#> Error : shrs: a count is larger than the total, so these are not shares of it.
```
