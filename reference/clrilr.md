# Log-ratio transforms of compositions

Element shares from portable X-ray fluorescence of bone or soil, or any
parts of a whole, live on the simplex: they sum to one, and ordinary
statistics on the raw shares mislead. The centred log-ratio (`"clr"`)
takes each part's log over the geometric mean of its row; the isometric
log-ratio (`"ilr"`, Egozcue et al. 2003) gives `d - 1` coordinates in an
orthonormal basis (the default pivot basis: coordinate `i` contrasts the
first `i` parts with part `i + 1`), on which the usual methods apply.
Each row is closed to sum 1 after `ps` is added to every part, so that a
part below the detection limit (0) has a log.

## Usage

``` r
clrilr(x, type = "clr", ps = 1e-06)
```

## Arguments

- x:

  The compositions: a numeric matrix or data frame of numbers from 0,
  one row per sample and one column per part, at least two parts.

- type:

  `"clr"` or `"ilr"`.

- ps:

  The pseudocount added to every part: one number from 0.

## Value

A matrix with one row per sample: `d` columns named as the parts
(`"clr"`), or `d - 1` columns `ilr1`, `ilr2`, ... (`"ilr"`).

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.
[doi:10.1007/978-94-009-4109-0](https://doi.org/10.1007/978-94-009-4109-0)

Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and
Barcelo-Vidal, C. (2003). Isometric logratio transformations for
compositional data analysis. *Mathematical Geology*, 35(3), 279-300.
[doi:10.1023/A:1023818214614](https://doi.org/10.1023/A%3A1023818214614)

## Examples

``` r
x <- rbind(c(Ca = 60, P = 25, Fe = 10, Zn = 5),
           c(Ca = 80, P = 12, Fe = 6, Zn = 2))
clrilr(x)
#>            Ca         P         Fe        Zn
#> [1,] 1.288034 0.4125649 -0.5037257 -1.196873
#> [2,] 2.044066 0.1469466 -0.5462005 -1.644812
clrilr(x, "ilr")
#>           ilr1     ilr2     ilr3
#> [1,] 0.6190499 1.105557 1.382030
#> [2,] 1.3414664 1.340448 1.899266
```
