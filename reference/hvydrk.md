# Heavy drinking, and drinking above the low-risk guideline

Two flags from a number of standard drinks and a sex, as the Canadian
surveys of alcohol and drug use derive them (Health Canada 2022).
`hvydrk()`: heavy drinking, four or more drinks for women and five or
more for men on one occasion, the measure of Wechsler et al. (1995).
`lrdg()`: more drinks in a day than the daily limit of Canada's Low-Risk
Alcohol Drinking Guidelines for chronic harm, two for women and three
for men (Stockwell et al. 2012).

## Usage

``` r
hvydrk(drinks, sex, f = 4, m = 5)

lrdg(drinks, sex, f = 2, m = 3)
```

## Arguments

- drinks:

  The number of standard drinks: on one occasion (`hvydrk()`), in a day
  (`lrdg()`).

- sex:

  `"male"` or `"female"`, one or one per value of `drinks`.

- f, m:

  The thresholds for women and for men.

## Value

A logical vector: `TRUE` where the drinking is heavy (`hvydrk()`) or
over the limit (`lrdg()`); `NA` where `drinks` is missing.

## Details

The first counts the threshold itself as heavy (`>=`); the second is
over the limit only above it (`>`). Both thresholds are arguments, for
another definition or another country's guideline. The guidelines of
2011 that `lrdg()` follows were replaced in 2023 by Canada's Guidance on
Alcohol and Health, which gives a continuum of risk by drinks a week and
no daily limit by sex; the survey variables, and so the defaults here,
are those of 2011.

Errors: `drinks` not numeric or negative; `sex` not `"male"` or
`"female"`; `f` or `m` not one number from 0; `drinks` and `sex` of
different lengths (a length of 1 is repeated).

A `clpre` list from
[`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
may be given as `drinks`: its elements fill the arguments.

## References

Wechsler, H., Dowdall, G. W., Davenport, A. and Rimm, E. B. (1995). A
gender-specific measure of binge drinking among college students.
*American Journal of Public Health*, 85(7), 982-985.
[doi:10.2105/AJPH.85.7.982](https://doi.org/10.2105/AJPH.85.7.982)

Stockwell, T., Butt, P., Beirness, D., Gliksman, L. and Paradis, C.
(2012). The basis for Canada's new low-risk drinking guidelines: a
relative risk approach to estimating hazardous levels and patterns of
alcohol use. *Drug and Alcohol Review*, 31(2), 126-134.
[doi:10.1111/j.1465-3362.2011.00342.x](https://doi.org/10.1111/j.1465-3362.2011.00342.x)

Health Canada (2022). *Canadian Postsecondary Education Alcohol and Drug
Use Survey, 2021-22: Public Use Microdata File user guide*. Table 3,
derived variables `alc12_30d_prev`, `u_d_lrdg` and `alc13_lrdg`.

## See also

[`ebac()`](https://rootcoder007.github.io/crimlab/reference/ebac.md) and
[`stddrk()`](https://rootcoder007.github.io/crimlab/reference/stddrk.md);
[`auditc()`](https://rootcoder007.github.io/crimlab/reference/auditc.md)
for the screening questionnaires.

## Examples

``` r
d <- c(2, 3, 4, 5, NA)
hvydrk(d, "female")
#> [1] FALSE FALSE  TRUE  TRUE    NA
hvydrk(d, "male")
#> [1] FALSE FALSE FALSE  TRUE    NA
lrdg(d, c("female", "male", "male", "female", "male"))
#> [1] FALSE FALSE  TRUE  TRUE    NA
# another definition: six drinks for everyone
hvydrk(d, "female", f = 6, m = 6)
#> [1] FALSE FALSE FALSE FALSE    NA
try(hvydrk(-1, "male"))
#> Error : hvydrk: `drinks` must be numbers from 0.
# a survey's public file, by its preset: the codes 997 to 999 are not
# answers, and the sex is coded 1 and 2
d <- data.frame(alc13a = c(3, 5, 997, 2), dvdemq01 = c(1, 2, 1, 2))
lrdg(clprep(d, "cpads-lrdg"))
#> [1]  TRUE  TRUE FALSE
```
