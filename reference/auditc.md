# Scores of the AUDIT, the DAST-10, the CUDIT-R and the ASSIST

Scores four screening questionnaires from their item answers, one row
per respondent and one column per item, in the order of the
questionnaire. An answer outside the scoring of its item is refused, not
scored: a miscoded column must not pass for a low score.

## Usage

``` r
auditc(x, sex = NULL)

dast10(x)

cudit(x)

canast(x)
```

## Arguments

- x:

  The answers: a numeric matrix or data frame, a row per respondent and
  a column per item.

- sex:

  `auditc()`: `"male"` or `"female"`, one or one per row, for the
  AUDIT-C screen; `NULL` (the default) leaves it out.

## Value

A data frame with a row per respondent. `auditc()`: `total` (`NA` with
three items), `auditc`, `zone` (`"low"`, `"hazardous"`, `"harmful"`,
`"dependence"`; `NA` with three items) and `cpos`, the AUDIT-C screen
(`NA` without `sex`). `dast10()`: `total` and `level` (`"none"`,
`"low"`, `"moderate"`, `"substantial"`, `"severe"`). `cudit()`: `total`
and `level` (`"low"`, `"hazardous"`, `"disorder"`). `canast()`: `total`
and `level` (`"low"`, `"moderate"`, `"high"`).

## Details

`auditc()`: the Alcohol Use Disorders Identification Test (Saunders et
al. 1993). Ten items scored 0 to 4 (items 9 and 10: 0, 2 or 4), total 0
to 40, and the zones of the World Health Organization's manual: 0 to 7
low risk, 8 to 15 hazardous, 16 to 19 harmful, 20 and over possible
dependence. The first three items are the AUDIT-C (Bush et al. 1998), 0
to 12; with `sex` its screen is positive from 4 for men and from 3 for
women (Bradley et al. 2003). Three columns are scored as the AUDIT-C
alone.

`dast10()`: the Drug Abuse Screening Test in its ten-item form (Skinner
1982). Ten yes or no items, 1 for yes; item 3 ("are you always able to
stop using drugs when you want to?") counts when the answer is no. Total
0 to 10: 0 none, 1 to 2 low, 3 to 5 moderate, 6 to 8 substantial, 9 to
10 severe.

`cudit()`: the Cannabis Use Disorders Identification Test, revised
(Adamson et al. 2010). Eight items scored 0 to 4 (item 8: 0, 2 or 4),
total 0 to 32: 8 and over hazardous use, 12 and over a possible cannabis
use disorder.

`canast()`: the cannabis score of the Alcohol, Smoking and Substance
Involvement Screening Test of the World Health Organization (WHO ASSIST
Working Group 2002), as the Canadian surveys of alcohol and drug use
derive it (Health Canada 2022, `can_assist`). Six items. The first four
are answered 0 never, 1 once or twice, 2 monthly, 3 weekly, 4 daily or
almost daily, over the past three months: use, a strong desire to use,
problems from use, and failing what was expected. The last two are
answered 0 never, 1 yes but not in the past three months, 2 yes in the
past three months: the concern of others, and a failed attempt to cut
down. The answers are weighted as the test weights them (0, 2, 3, 4, 6;
0, 3, 4, 5, 6; 0, 4, 5, 6, 7; 0, 5, 6, 7, 8; and 0, 3, 6 for each of the
last two), total 0 to 39: 0 to 3 low risk, 4 to 26 moderate, 27 and over
high (Humeniuk et al. 2008).

A respondent with a missing item has missing scores. These are screens:
a score says who to assess, not who has a disorder.

Errors: `x` not a matrix or data frame of numbers with the columns of
the questionnaire (`canast()`: six); an answer outside its item's
scoring; `sex` not `"male"` or `"female"`, one or one per row.

`canast()` takes a `clpre` list from
[`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
as `x` too: its elements fill the arguments.

## References

Saunders, J. B., Aasland, O. G., Babor, T. F., de la Fuente, J. R. and
Grant, M. (1993). Development of the Alcohol Use Disorders
Identification Test (AUDIT): WHO collaborative project on early
detection of persons with harmful alcohol consumption-II. *Addiction*,
88(6), 791-804.
[doi:10.1111/j.1360-0443.1993.tb02093.x](https://doi.org/10.1111/j.1360-0443.1993.tb02093.x)

Bush, K., Kivlahan, D. R., McDonell, M. B., Fihn, S. D. and Bradley, K.
A. (1998). The AUDIT alcohol consumption questions (AUDIT-C): an
effective brief screening test for problem drinking. *Archives of
Internal Medicine*, 158(16), 1789-1795.
[doi:10.1001/archinte.158.16.1789](https://doi.org/10.1001/archinte.158.16.1789)

Bradley, K. A., Bush, K. R., Epler, A. J., Dobie, D. J., Davis, T. M.,
Sporleder, J. L., Maynard, C., Burman, M. L. and Kivlahan, D. R. (2003).
Two brief alcohol-screening tests from the Alcohol Use Disorders
Identification Test (AUDIT). *Archives of Internal Medicine*, 163(7),
821-829.
[doi:10.1001/archinte.163.7.821](https://doi.org/10.1001/archinte.163.7.821)

Skinner, H. A. (1982). The drug abuse screening test. *Addictive
Behaviors*, 7(4), 363-371.
[doi:10.1016/0306-4603(82)90005-3](https://doi.org/10.1016/0306-4603%2882%2990005-3)

Adamson, S. J., Kay-Lambkin, F. J., Baker, A. L., Lewin, T. J.,
Thornton, L., Kelly, B. J. and Sellman, J. D. (2010). An improved brief
measure of cannabis misuse: the Cannabis Use Disorders Identification
Test-Revised (CUDIT-R). *Drug and Alcohol Dependence*, 110(1-2),
137-143.
[doi:10.1016/j.drugalcdep.2010.02.017](https://doi.org/10.1016/j.drugalcdep.2010.02.017)

WHO ASSIST Working Group (2002). The Alcohol, Smoking and Substance
Involvement Screening Test (ASSIST): development, reliability and
feasibility. *Addiction*, 97(9), 1183-1194.
[doi:10.1046/j.1360-0443.2002.00185.x](https://doi.org/10.1046/j.1360-0443.2002.00185.x)

Humeniuk, R., Ali, R., Babor, T. F., Farrell, M., Formigoni, M. L.,
Jittiwutikarn, J., de Lacerda, R. B., Ling, W., Marsden, J., Monteiro,
M., Nhiwatiwa, S., Pal, H., Poznyak, V. and Simon, S. (2008). Validation
of the Alcohol, Smoking and Substance Involvement Screening Test
(ASSIST). *Addiction*, 103(6), 1039-1047.
[doi:10.1111/j.1360-0443.2007.02114.x](https://doi.org/10.1111/j.1360-0443.2007.02114.x)

Health Canada (2022). *Canadian Postsecondary Education Alcohol and Drug
Use Survey, 2021-22: Public Use Microdata File user guide*. Table 3,
derived variable `can_assist`.

## See also

[`ebac()`](https://rootcoder007.github.io/crimlab/reference/ebac.md),
[`stddrk()`](https://rootcoder007.github.io/crimlab/reference/stddrk.md),
[`hvydrk()`](https://rootcoder007.github.io/crimlab/reference/hvydrk.md)
and
[`lrdg()`](https://rootcoder007.github.io/crimlab/reference/hvydrk.md).

## Examples

``` r
# made-up answers of three respondents
a <- rbind(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
           c(3, 2, 2, 1, 1, 0, 1, 1, 0, 2),
           c(4, 4, 4, 3, 3, 2, 3, 3, 4, 4))
auditc(a, sex = c("female", "male", "male"))
#>   total auditc       zone  cpos
#> 1     1      1        low FALSE
#> 2    13      7  hazardous  TRUE
#> 3    34     12 dependence  TRUE
# the first three items alone are the AUDIT-C
auditc(a[, 1:3], sex = "female")
#>   total auditc zone  cpos
#> 1    NA      1 <NA> FALSE
#> 2    NA      7 <NA>  TRUE
#> 3    NA     12 <NA>  TRUE
dast10(rbind(c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
             c(1, 1, 0, 1, 1, 0, 1, 1, 0, 1)))
#>   total       level
#> 1     1         low
#> 2     8 substantial
cudit(rbind(c(1, 0, 0, 0, 0, 0, 0, 0),
            c(4, 3, 2, 1, 2, 1, 2, 4)))
#>   total    level
#> 1     1      low
#> 2    19 disorder
# weekly use, a monthly desire, no problems: moderate risk
canast(rbind(c(1, 0, 0, 0, 0, 0), c(3, 2, 0, 0, 0, 0),
             c(4, 4, 3, 3, 2, 2)))
#>   total    level
#> 1     2      low
#> 2     8 moderate
#> 3    37     high
try(auditc(cbind(a[, 1:9], 3)))
#> Error : auditc: item 10 must be scored 0, 2 or 4.
```
