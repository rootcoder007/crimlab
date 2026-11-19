# Do decision makers decide alike?

When cases are shared among several decision makers (reviewers, judges,
panels), the share of their decisions that go one way should differ by
no more than chance if the cases are alike. `dmvar()` gives each
decision maker's rate with its exact interval and the share of all such
decisions that are theirs, the chi-square test that the rates are equal,
the same test without the decision makers who decided few cases, and the
test with each decision maker left out in turn, which shows whether one
of them makes the difference. Sprott, Doob and Iftene (2021, Tables 10
and 11) found the share of reviews that kept a person in a Structured
Intervention Unit to run from 37.5% to over 75% across twelve decision
makers, with and without the one who had decided fewest.

## Usage

``` r
dmvar(dec, dm, out, nmin = 5, level = 0.95)
```

## Arguments

- dec:

  The decision in each case: labels.

- dm:

  The decision maker of each case: labels.

- out:

  The decision whose rate is compared: one of the labels of `dec`.

- nmin:

  The cases below which a decision maker is left out of the second test:
  a whole number from 1.

- level:

  The level of the exact intervals: one number between 0 and 1.

## Value

A list: `tab`, a data frame with one row per decision maker (`dm`, `n`,
`k` the decisions that were `out`, `rate`, `lo`, `hi`, `share` of all
such decisions); `test`, two rows (`set`: `"all"` and `"n >= nmin"`,
with `m` the decision makers, `stat`, `df`, `p` and `low` the cells with
an expected count under 5); `drop`, one row per decision maker left out
(`dm`, `stat`, `df`, `p`); `out`; and `level`.

## Details

The test takes the cases as shared without regard to who decides; where
they are not, a difference may be one of the cases and not of the
decision makers
([`jdgiv()`](https://rootcoder007.github.io/crimlab/reference/jdgiv.md)
and
[`jmono()`](https://rootcoder007.github.io/crimlab/reference/jmono.md)
treat that question).

## References

Sprott, J. B., Doob, A. N. and Iftene, A. (2021). *Do Independent
External Decision Makers Ensure that "An Inmate's Confinement in a
Structured Intervention Unit Is to End as Soon as Possible"?* Schulich
School of Law, Dalhousie University, 10 May 2021.

## Examples

``` r
dm <- rep(c("a", "b", "c", "d"), c(30, 28, 25, 4))
dec <- rep(c("remain", "remove", "remain", "remove", "remain", "remove",
             "remain"), c(24, 6, 10, 18, 15, 10, 4))
dmvar(dec, dm, out = "remain")
#> $tab
#>   dm  n  k      rate        lo        hi     share
#> 1  a 30 24 0.8000000 0.6143335 0.9228645 0.4528302
#> 2  b 28 10 0.3571429 0.1864067 0.5593497 0.1886792
#> 3  c 25 15 0.6000000 0.3866535 0.7887452 0.2830189
#> 4  d  4  4 1.0000000 0.3976354 1.0000000 0.0754717
#> 
#> $test
#>         set m     stat df           p low
#> 1       all 4 14.63426  3 0.002157435   2
#> 2 n >= nmin 3 11.75893  2 0.002796286   0
#> 
#> $drop
#>   dm      stat df            p
#> 1  a  7.270408  2 0.0263785505
#> 2  b  4.356395  2 0.1132454520
#> 3  c 14.672556  2 0.0006514707
#> 4  d 11.758926  2 0.0027962859
#> 
#> $out
#> [1] "remain"
#> 
#> $level
#> [1] 0.95
#> 
```
