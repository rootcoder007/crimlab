# Machine-checked proofs behind the package

Names the 'Lean' theorem each function rests on. The theorems ship with
the package: `file` is the path below
`system.file(package = "crimlab")`, and each file compiles against
'Mathlib'.

## Usage

``` r
proof(fn = NULL)
```

## Arguments

- fn:

  `NULL` for every function, or function names.

## Value

A data frame: `fn`, `theorem` and `file`.

## Examples

``` r
proof("vodtst")
#>       fn                                 theorem           file
#> 1 vodtst Crimlab.Veil.odds_ratio_identifies_bias lean/Veil.lean
#> 2 vodtst                 Crimlab.Veil.share_sign lean/Veil.lean
#> 3 vodtst   Crimlab.Veil.exposure_shift_confounds lean/Veil.lean
```
