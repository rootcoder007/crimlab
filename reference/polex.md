# A simulated policy panel for the examples

Twenty districts in four regions over 36 months; a reform adopted region
by region in three stages (months 13, 19 and 25), with the fourth region
never adopting it; crime and arrest rates and counts. The numbers are
made up, with a known effect of the reform (a fall of 3 in the crime
rate per 100,000), so the examples and tests can check that each
estimator recovers it.

## Usage

``` r
polex(seed = 20240101)
```

## Arguments

- seed:

  The seed (default 20240101).

## Value

A data frame: `district`, `region`, `month`, `refmonth` (`NA` for
never), `reformed` (0 or 1), `population`, `logpop`, `crimes`,
`crmrate`, `arrrate`.

## Details

The panel is simulated: it describes no real place, district or reform,
and its numbers must not be read as data. The same seed gives the same
panel in every front end of the package (the draws are those of R's
default generator after `set.seed(seed)`), and the policy functions'
examples use it.

## Examples

``` r
# the panel is simulated: made-up numbers with a known effect of -3
p <- polex()
head(p)
#>   district region month refmonth reformed population   logpop crimes crmrate
#> 1      D01     R1     1       13        0     152152 11.93264     80  39.732
#> 2      D01     R1     2       13        0     152152 11.93264     56  42.708
#> 3      D01     R1     3       13        0     152152 11.93264     68  41.228
#> 4      D01     R1     4       13        0     152152 11.93264     68  43.296
#> 5      D01     R1     5       13        0     152152 11.93264     62  42.036
#> 6      D01     R1     6       13        0     152152 11.93264     59  40.197
#>   arrrate
#> 1  15.394
#> 2  16.387
#> 3  16.043
#> 4  17.148
#> 5  15.443
#> 6  15.568
# the same seed gives the same panel
identical(polex(7), polex(7))
#> [1] TRUE
```
