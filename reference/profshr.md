# The share of a population with a characteristic, from the census profile

The share of a group with a characteristic of the 2021 Census Indigenous
Population Profile (foster children among children, speakers of an
Indigenous language, households needing major repairs, ...), with an
approximate confidence interval that carries the long-form sample and
the random rounding of the counts.

## Usage

``` r
profshr(
  x,
  num,
  den = NULL,
  smpfrac = 0.25,
  deff = 1,
  conflvl = 0.95,
  roubase = 5
)
```

## Arguments

- x:

  A
  [`cnsprof`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md).

- num:

  The characteristic counted: its `chaid`, or its name exactly.

- den:

  The characteristic it is a share of (id or name); default its parent.

- smpfrac:

  The long-form sampling fraction (default 0.25); full-count
  characteristics use 1 whatever is given.

- deff:

  Design effect (default 1).

- conflvl:

  Confidence level of the share's limits (default 0.95).

- roubase:

  The base of the random rounding (default 5, Statistics Canada's); 1
  for counts that are not rounded.

## Value

A data frame of class `profshr` with `geo`, `dguid`, `residence`,
`gender`, `age`, `group`, `numerator`, `denominator`, `estimate`
(percent), `lower`, `upper`, `quality`; ready for
[`indgap`](https://rootcoder007.github.io/crimlab/reference/indgap.md).

## Details

The share is the numerator's count over the denominator's (by default
the numerator's parent in the profile's hierarchy), for each geography,
gender, age group and group. Statistics Canada does not publish standard
errors for the profile, so the interval is an approximation built from
what is known about how the counts were made:

- a long-form ("25% sample") count is estimated from a sample of about
  one household in four, so the share has variance \\d (1 - f) p (1 - p)
  / (f N)\\ for a sampling fraction \\f\\ (`smpfrac`, 0.25) of a
  population of \\N\\ (the denominator) and a design effect \\d\\
  (`deff`, 1 by default: the clustering of people in households raises
  it, and a value of 1.5 to 2 is prudent for characteristics shared
  within a household). On reserves (rows whose `residence` is
  `"On reserve"`) every household got the long form, so \\f = 1\\ there
  and the sampling term vanishes;

- every count is randomly rounded, up or down without bias, to a
  multiple of \\b\\ (`roubase`, 5). Averaged over the remainders, that
  adds variance \\c = (b^2 - 1)/6\\ to each count (4 for \\b = 5\\), so
  \\c (1 + p^2) / N^2\\ to the share (the only error of a full-count
  share).

In all, \$\$V = \frac{d (1 - f) p (1 - p)}{f N} + \frac{c (1 +
p^2)}{N^2}.\$\$ The interval is Wilson's at the effective sample size
\\n_e = \max(p(1-p), 1/N) / V\\ (the floor \\1/N\\ keeps a share of 0 or
100% from an infinite \\n_e\\), so it stays inside 0 and 100. It does
not carry non-response (see `tnrlf`) or the edit and imputation error,
so it is narrower than the full error; treat it as a lower bound on the
uncertainty.

## References

Statistics Canada (2022). Guide to the Census of Population, 2021,
chapter 9 (data quality) and chapter 10 (confidentiality and random
rounding). Catalogue no. 98-304-X.

Wilson EB (1927). Probable inference, the law of succession, and
statistical inference. Journal of the American Statistical Association
22(158):209-212.
[doi:10.1080/01621459.1927.10502953](https://doi.org/10.1080/01621459.1927.10502953)

## See also

[`cnsprof`](https://rootcoder007.github.io/crimlab/reference/cnsprof.md).

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
# the share with no certificate, of the population aged 25 and over
profshr(p, num = "No certificate, diploma or degree")
#> Share of a population with a characteristic
#> -------------------------------------------
#> Numerator:   No certificate, diploma or degree
#> Denominator: Total - Highest certificate, diploma or degree for the population aged 25 years and over in private households - 25% sample data
#> Rows:        2
#> Estimate:    percent
#> 
#>      geo  dguid                            residence          gender                age                group  numerator  denominator  estimate      interval  quality
#> Lakeside  X0001  Total - Residence on or off reserve  Total - Gender  25 years and over  Indigenous identity        150          700     21.43  [16.6, 27.2]       ok
#> Lakeside  X0001  Total - Residence on or off reserve  Total - Gender  25 years and over       Non-Indigenous        400         3500     11.43  [9.73, 13.4]       ok
#> 
#> Highest share: Lakeside, 25 years and over, Indigenous identity (21.43); lowest: Lakeside, 25 years and over, Non-Indigenous (11.43).
#> Their intervals do not overlap: the gap is more than sampling error.
#> 
#> The interval carries the long-form sample and the random rounding,
#> not non-response: read it as a lower bound on the uncertainty.
```
