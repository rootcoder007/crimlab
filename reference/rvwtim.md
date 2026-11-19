# Was the review held, and in time?

A stay in a restrictive unit is to be reviewed by a set day. `rvwtim()`
takes the length of each stay and the day of its first review, and gives
how many stays of each length were never reviewed, how many of those
past the day a review was due had none, and how long the reviews that
were held took. Sprott, Doob and Iftene (2021, Tables 2 and 15) asked
this of the reviews by Independent External Decision Makers of stays in
Canada's Structured Intervention Units; the question is the same for any
review a rule sets a day for (a bail review, a segregation review, a
hearing).

## Usage

``` r
rvwtim(days, first, due, by = NULL, cuts = NULL, level = 0.95)
```

## Arguments

- days:

  The length of each stay in days: numbers from 0.

- first:

  The day of the stay on which it was first reviewed, `NA` for a stay
  never reviewed: numbers from 0.

- due:

  The day by which a review is due: one number above 0.

- by:

  `NULL`, or the group of each stay: labels.

- cuts:

  `NULL` (two bands of length, up to `due` and over it) or the lengths
  that part the bands: increasing numbers above 0.

- level:

  The level of the intervals: one number between 0 and 1.

## Value

A list: `cov`, a data frame with one row per group (then `"all"`) and
band of length (`group`, `band`, `n`, `rev` the stays reviewed, `never`,
`share` the share never reviewed, `lo`, `hi`); `sum`, one row per group
(`group`, `n`, `rev`, `past` the stays longer than `due`, `pastnever`
those of them never reviewed, `share`, `lo`, `hi`, `late` the stays
first reviewed after `due`, `med` the median day of the first review);
`due`; and `level`.

## Details

A stay still open is entered with its length to the day the data end.
The shares carry exact (Clopper-Pearson) intervals.

## References

Sprott, J. B., Doob, A. N. and Iftene, A. (2021). *Do Independent
External Decision Makers Ensure that "An Inmate's Confinement in a
Structured Intervention Unit Is to End as Soon as Possible"?* Schulich
School of Law, Dalhousie University, 10 May 2021.

## Examples

``` r
days <- c(12, 40, 75, 90, 130, 150, 64, 20)
first <- c(NA, 31, 48, NA, 61, NA, 33, NA)
rvwtim(days, first, due = 60, cuts = c(60, 120))
#> $cov
#>   group           band n rev never     share          lo        hi
#> 1   all       up to 60 3   1     2 0.6666667 0.094299324 0.9915962
#> 2   all over 60 to 120 3   2     1 0.3333333 0.008403759 0.9057007
#> 3   all       over 120 2   1     1 0.5000000 0.012579117 0.9874209
#> 
#> $sum
#>   group n rev past pastnever share         lo        hi late  med
#> 1   all 8   4    5         2   0.4 0.05274495 0.8533672    1 40.5
#> 
#> $due
#> [1] 60
#> 
#> $level
#> [1] 0.95
#> 
```
