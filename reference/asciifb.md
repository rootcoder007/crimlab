# Use text as it is, falling back to ASCII when it is not valid UTF-8

Returns `x` unchanged when it is valid UTF-8 (so an accented name is
kept), and transliterates to plain ASCII with
[`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md)
only the elements that are not valid UTF-8 (an encoding error), or every
element when `force = TRUE` (for an ASCII-only destination). Names keep
their accents wherever UTF-8 is supported and degrade gracefully where
it is not.

## Usage

``` r
asciifb(x, force = FALSE)
```

## Arguments

- x:

  A character vector (other types are converted with
  [`as.character()`](https://rdrr.io/r/base/character.html)).

- force:

  `TRUE` transliterates every element; `FALSE` (the default) only those
  that are not valid UTF-8. Anything but `TRUE` counts as `FALSE`.

## Value

A character vector the length of `x`.

## Details

Each element is tested on its bytes with
[`validUTF8()`](https://rdrr.io/r/base/validUTF8.html), so the outcome
does not depend on the session locale. A missing value is kept as `NA`
unless `force = TRUE`, when it becomes `""`.

## See also

[`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md).

## Examples

``` r
asciifb("\u00c1ngela")
#> [1] "Ángela"
asciifb("\u00c1ngela", force = TRUE)
#> [1] "Angela"
bad <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
validUTF8(bad)
#> [1] FALSE
asciifb(bad)
#> [1] "caf"
```
