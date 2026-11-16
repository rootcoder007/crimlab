# One call on a published table: check, change, rates and trend

The questions asked of a table an agency publishes are nearly always the
same: are its categories sound, what changed from one period to the next
and how sure can one be, what are the rates where there is a population,
and does a short series trend. `tabrpt()` answers them in one call on a
table in long form (one row per period and group), with the steps a
reviewer needs and that are easy to forget on a deadline: exact
intervals and p-values adjusted over every comparison, and the envelope
that rounding in the release allows.

## Usage

``` r
tabrpt(
  data,
  value,
  prd,
  by = NULL,
  pop = NULL,
  rnd = NULL,
  kind = "nearest",
  per = 1000,
  level = 0.95,
  adj = "BH",
  alpha = 0.05
)

# S3 method for class 'tabrpt'
print(x, ...)
```

## Arguments

- data:

  The table: a data frame in long form.

- value:

  The column of the counts: one name.

- prd:

  The column of the periods (numbers: a year, a month index): one name.

- by:

  `NULL`, or the columns of the groups: names.

- pop:

  `NULL`, or the column of the populations: one name.

- rnd:

  `NULL`, or the rounding base of the published counts.

- kind:

  `"nearest"` or `"random"` rounding.

- per:

  The rates are given per `per` of `pop`: one positive number.

- level:

  The confidence level of every interval.

- adj:

  The adjustment of the p-values for the number of comparisons, as
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) names it.

- alpha:

  The level at which an adjusted p-value is marked.

- x:

  A `tabrpt` result.

- ...:

  Not used.

## Value

An object of class `tabrpt`, a list: `check` (the
[`catchk()`](https://rootcoder007.github.io/crimlab/reference/catchk.md)
result, `NULL` without `by`), `change` (the
[`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
table), `rates` (the
[`grprr()`](https://rootcoder007.github.io/crimlab/reference/grprr.md)
result, `NULL` without `pop`) and `trend` (a list of
[`trndrr()`](https://rootcoder007.github.io/crimlab/reference/trndrr.md)
results named by group, `NULL` for a group too short to test).

## Details

The steps:
[`catchk()`](https://rootcoder007.github.io/crimlab/reference/catchk.md)
on the group columns (with `by`);
[`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)
on the counts by period within each group, with `rnd`, `kind`, `level`,
`adj` and `alpha`; with `pop`,
[`grprr()`](https://rootcoder007.github.io/crimlab/reference/grprr.md)
of the counts over the population pooled by group (by period when there
is no `by`), per `per`; and for each group with three periods or more
and a count above 0,
[`trndrr()`](https://rootcoder007.github.io/crimlab/reference/trndrr.md)
of its totals by period (over its population by period with `pop`). A
group is the combination of its `by` columns, labelled with them joined
by `" / "`. The refusals and warnings of the steps are given under the
name `tabrpt`.

## Examples

``` r
d <- data.frame(year = rep(2019:2023, 2),
                svc = rep(c("A", "B"), each = 5),
                n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33),
                pop = rep(c(50000, 12000), each = 5))
r <- tabrpt(d, "n", "year", by = "svc", pop = "pop", rnd = 5)
r
#> tabrpt: 10 group-periods, 8 changes tested, 0 marked after adjustment
#> groups: no problem found
#> rates per 1000: A 2.38, B 2.62
#> trend A: rate ratio per period 1.01, p 0.839
#> trend B: rate ratio per period 1.01, p 0.901
r$change[, c("group", "prd", "pctchg", "padj")]
#> yoypx: change on the earlier period
#>  group  prd     pctchg      padj
#>      A 2019         NA        NA
#>      A 2020   9.166667 0.7039955
#>      A 2021 -25.190840 0.1369352
#>      A 2022   7.142857 0.7700275
#>      A 2023  33.333333 0.1369352
#>      B 2019         NA        NA
#>      B 2020  -6.666667 0.8956832
#>      B 2021  46.428571 0.2960644
#>      B 2022 -39.024390 0.1707133
#>      B 2023  32.000000 0.5730293
```
