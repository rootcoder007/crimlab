# Collider bias of arrest data

Two causes of arrest with prevalences `a` and `e`, independent in the
population, and a background cause of prevalence `bg`; a person is
arrested when any is present. Among arrestees the two causes are no
longer independent: their odds ratio is exactly `bg`, below 1, and the
negative association grows as arrests for other reasons become rare.
Studying arrestees alone manufactures it (Hernan and Robins, Fine Point
8.2).

## Usage

``` r
colarr(a, e, bg)
```

## Arguments

- a, e:

  The prevalences of the two causes: numbers strictly between 0 and 1,
  one per scenario or one for all.

- bg:

  The prevalence of the background cause: numbers from 0 to 1, one per
  scenario or one for all.

## Value

A data frame with one row per scenario: `a`, `e`, `bg`, `or` the odds
ratio among arrestees, `lor` its log, and the shares of arrestees with
both causes (`c11`), `a` only (`c10`), `e` only (`c01`) and neither
(`c00`).

## References

Hernan, M. A. and Robins, J. M. (2020). *Causal Inference: What If*.
Chapman and Hall/CRC. <https://miguelhernan.org/whatifbook>

## Examples

``` r
colarr(a = 0.3, e = 0.2, bg = c(0.5, 0.1, 0.01))
#>     a   e   bg   or        lor        c11       c10       c01        c00
#> 1 0.3 0.2 0.50 0.50 -0.6931472 0.08333333 0.3333333 0.1944444 0.38888889
#> 2 0.3 0.2 0.10 0.10 -2.3025851 0.12096774 0.4838710 0.2822581 0.11290323
#> 3 0.3 0.2 0.01 0.01 -4.6051702 0.13464991 0.5385996 0.3141831 0.01256732
```
