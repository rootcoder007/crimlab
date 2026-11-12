# Bounds on a sentence's effect without assumptions

Each person receives one of two sentences and shows the outcome (say,
reconviction) under that one only. With no assumption about who receives
which, the chance of the outcome under a sentence is known only up to an
interval whose width is the share who received the other sentence, and
the effect of one sentence against the other lies in an interval of
width exactly 1 that always contains 0: observational data alone cannot
even sign it. See `proof("sntbnd")`.

## Usage

``` r
sntbnd(y, z, treat = NULL, w = NULL)
```

## Arguments

- y:

  Outcomes: 0 or 1 for each person.

- z:

  Sentences: a vector with exactly two distinct values, same length:
  character, numbers, logical or a factor (its labels; level order and
  whether it is ordered are ignored), not dates nor a labelled vector of
  codes. When `y`, `z` (and `w`) are all named with unique names, they
  are matched by name.

- treat:

  The value of `z` taken as the treatment, matched exactly. `NULL` is
  allowed only where the treatment is explicit in `z`: `TRUE` for a
  logical `z`, `1` for a `z` of 0 and 1; otherwise it must be given,
  since a treatment chosen by sorting the values would be arbitrary.

- w:

  `NULL` for equal weights, or one non-negative weight per person with a
  positive total.

## Value

A data frame with three rows, the comparison sentence, the treatment and
`effect` (treatment minus comparison): `term`, `share` (the weighted
share sentenced; `NA` for the effect), `lower`, `upper` and `observed`
(the observed rate under each sentence and their difference).

## References

Manski, C. F. and Nagin, D. S. (1998). Bounding disagreements about
treatment effects: a case study of sentencing and recidivism.
*Sociological Methodology*, 28(1), 99-137.
[doi:10.1111/0081-1750.00043](https://doi.org/10.1111/0081-1750.00043)

## Examples

``` r
set.seed(1)
z <- sample(c("community", "custody"), 500, TRUE, c(0.6, 0.4))
y <- rbinom(500, 1, ifelse(z == "custody", 0.55, 0.35))
sntbnd(y, z, treat = "custody")
#>        term share  lower upper  observed
#> 1 community 0.624  0.208 0.584 0.3333333
#> 2   custody 0.376  0.192 0.816 0.5106383
#> 3    effect    NA -0.392 0.608 0.1773050
```
