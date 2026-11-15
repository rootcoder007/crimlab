# A Community Well-Being score from the census profile

Indigenous Services Canada's Community Well-Being (CWB) index
([`cwbscr`](https://rootcoder007.github.io/crimlab/reference/cwbscr.md))
computed from the public 2021 Census Indigenous Population Profile, for
each geography and group.

## Usage

``` r
cnscwb(
  x,
  dwell = NULL,
  groups = c("First Nations", "Métis", "Inuit", "Non-Indigenous"),
  bounds = c(2650, 75000),
  resid = "Total - Residence on or off reserve",
  minpop = 250
)
```

## Arguments

- x:

  A
  [`cnsprof`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md)
  of a population table.

- dwell:

  Optional
  [`cnsprof`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md)
  of the matching household table.

- groups:

  The groups to score (default First Nations, Metis, Inuit and
  non-Indigenous people).

- bounds:

  Income bounds of
  [`cwbscr`](https://rootcoder007.github.io/crimlab/reference/cwbscr.md).

- resid:

  On or off reserve, in the tables that split it (default the total).

- minpop:

  The smallest population scored (default 250; the Department uses 65).

## Value

A data frame of class `cnscwb`: `geo`, `dguid`, `group`; the components
`income`, `education`, `housing`, `labour` and the score `cwb` (0 to
100, `NA` where not scored); and the inputs `incpercap` (dollars),
`highsch25p`, `uni25p`, `population`, `participation`, `employment`,
`nomajrep`, `notcro` (percentages).

## Details

The Department computes the CWB from census microdata; the public
profile gives close but not the same inputs, and this function says
where it departs:

- income:

  per capita: the number of total income recipients aged 15 and over
  times their average total income, over the population in private
  households (all ages), as the CWB does.

- education:

  the shares with at least a high school diploma (all but those with no
  certificate, diploma or degree) and with a bachelor's degree or
  higher, both among those aged 25 and over: the profile has no age
  group of 20 and over, which the CWB uses for the first.

- labour force:

  the participation rate and the employment rate of participants (100
  less the unemployment rate) of those aged 25 to 64: the CWB uses 20 to
  64.

- housing:

  the shares of households with one person or fewer per room and needing
  no major repairs, from the household profile (`dwell`), which splits
  Indigenous from non-Indigenous households only, so First Nations,
  Metis and Inuit people all get the pooled Indigenous households'
  housing: the CWB counts persons, not households. Without `dwell` the
  housing component and the score are `NA`.

A group of fewer than `minpop` people in a geography gets no score. The
Department scores communities of 65 people or more; the default here is
250, because the public profile's counts are rounded to multiples of 5
and its long-form counts rest on a few dozen households at that size.
Geographies are matched by DGUID, so places that share a name stay
apart. The inputs are those of the 2021 profile (98-510-X2021001 to
-016): other editions name their characteristics differently. Read the
result as a CWB-style description of socio-economic conditions, close to
but not the Department's index, and never as a measure of the wellness
of First Nations, Metis or Inuit people as they define it (see
[`cwbscr`](https://rootcoder007.github.io/crimlab/reference/cwbscr.md)).

## References

Indigenous Services Canada (2019). Report on trends in First Nations
communities, 1981 to 2016 (Community Well-Being index).  
sac-isc.gc.ca/eng/1345816651029/1557323327644

## See also

[`cwbscr`](https://rootcoder007.github.io/crimlab/reference/cwbscr.md).

## Examples

``` r
# made-up numbers in the layout of a population profile table, for one
# place and two groups (the published tables hold many more rows)
nm <- c(
  "Total - Age groups of the population - 100% data",
  paste("Number of total income recipients aged 15 years and over in",
        "private households in 2020 - 25% sample data"),
  "Average total income in 2020 among recipients ($)",
  paste("Total - Highest certificate, diploma or degree for the population",
        "aged 25 years and over in private households - 25% sample data"),
  "  No certificate, diploma or degree",
  "  Bachelor's degree or higher",
  "Participation rate",
  "Unemployment rate"
)
x <- data.frame(
  GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
  AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3),
                rep("25 to 64 years", 2)),
  RESIDENCE = "Total - Residence on or off reserve",
  CHARACTERISTIC_ID = 1:8, CHARACTERISTIC_NAME = nm,
  `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14),
  SYMBOL_IND_IDENT = "",
  NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6),
  SYMBOL_NON_IND = "",
  check.names = FALSE
)
p <- cnsprof(x)
cnscwb(p, groups = c("Indigenous identity", "Non-Indigenous"), minpop = 0)
#> Community Well-Being index, from the census
#> -------------------------------------------
#> Rows:       2
#> Not scored: 2
#> Scale:      0 to 100 for each component and the score
#> 
#>      geo                group  income  education  housing  labour  cwb
#> Lakeside  Indigenous identity   71.06      58.10       NA      74   NA
#> Lakeside       Non-Indigenous   82.37      67.62       NA      80   NA
#> 
#> No group has a score to compare.
#> 
#> cwb is the mean of income, education, housing and labour force
#> activity (Indigenous Services Canada's Community Well-Being index).
#> CWB-style, close to but not the Department's index: the public
#> profile's inputs differ (ages 25 and over; housing from households,
#> pooled for First Nations, Métis and Inuit people).
#> Never a measure of wellness as First Nations, Métis or Inuit
#> people define it. A group below minpop people gets no score.
#> Housing needs the household table (dwell).
```
