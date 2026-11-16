# Association rules between two categorical variables

The rules the Apriori algorithm (Agrawal et al. 1996) finds when each
observation holds one level of each of two variables: `{} => {b}` and
`{a} => {b}` with `a` a level of one variable and `b` of the other, each
with its support, confidence, coverage and lift. Wongchanapai et al.
(2025) cut the rate of intimate partner violence and the annual PM2.5 of
New Jersey's municipalities into five groups each with
[`cut()`](https://rdrr.io/r/base/cut.html) and kept rules of "at least
10% support ... and at least 20% confidence" (`minsup = 0.1`,
`minconf = 0.2`, the defaults here) in R's arules package; with two
variables these are all its rules, the empty left-hand side included as
arules includes it by default (`empty`).

## Usage

``` r
asrule(x, y, minsup = 0.1, minconf = 0.2, empty = TRUE, nm = c("x", "y"))
```

## Arguments

- x, y:

  The two variables: factors (their levels, in order) or vectors (their
  sorted values), one value per observation, no `NA`.

- minsup, minconf:

  The least support and confidence, from 0 to 1.

- empty:

  `TRUE` to keep the rules with an empty left-hand side.

- nm:

  The names of the two variables in the items.

## Value

A data frame, one row per rule: `lhs`, `rhs` (items `"name=level"`, `""`
for the empty set), `support`, `confidence`, `coverage`, `lift` and
`count`.

## Details

Support is the share of observations holding both sides; confidence, the
support over that of the left-hand side (its coverage); lift, the
confidence over the support of the right-hand side, 1 for a side that
tells nothing of the other. Two levels of one variable never occur
together, so no rule joins them, and a pair of levels never observed
together is no rule.

## References

Agrawal, R., Mannila, H., Srikant, R., Toivonen, H. and Verkamo, A. I.
(1996). Fast discovery of association rules. In *Advances in Knowledge
Discovery and Data Mining*, 307-328. AAAI Press.

Wongchanapai, P., Green, R., Ku, C. S., Peltier, M. R. and Ranjan, S.
(2025). Air pollution and intimate partner violence in New Jersey:
examining the evidence. *Psychology of Violence*, 15(6), 682-688.
[doi:10.1037/vio0000653](https://doi.org/10.1037/vio0000653)

## See also

[`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md)

## Examples

``` r
set.seed(5)
pm <- cut(rnorm(200), 5, labels = c("vl", "l", "m", "mh", "vh"))
ipv <- cut(rnorm(200), 5, labels = c("vl", "l", "m", "mh", "vh"))
asrule(pm, ipv, nm = c("pm25", "ipv"))
#>        lhs     rhs support confidence coverage     lift count
#> 1           pm25=l   0.225  0.2250000    1.000 1.000000    45
#> 2           pm25=m   0.405  0.4050000    1.000 1.000000    81
#> 3          pm25=mh   0.260  0.2600000    1.000 1.000000    52
#> 4            ipv=m   0.375  0.3750000    1.000 1.000000    75
#> 5           ipv=mh   0.365  0.3650000    1.000 1.000000    73
#> 6   pm25=m   ipv=m   0.160  0.3950617    0.405 1.053498    32
#> 7    ipv=m  pm25=m   0.160  0.4266667    0.375 1.053498    32
#> 8   pm25=m  ipv=mh   0.160  0.3950617    0.405 1.082361    32
#> 9   ipv=mh  pm25=m   0.160  0.4383562    0.365 1.082361    32
#> 10 pm25=mh   ipv=m   0.105  0.4038462    0.260 1.076923    21
#> 11   ipv=m pm25=mh   0.105  0.2800000    0.375 1.076923    21
```
