# Transliterate text to plain ASCII

Converts text to plain 7-bit ASCII: the common Latin accents and
ligatures are folded to their nearest ASCII letters (an accented capital
A becomes `"A"`, sharp s `"ss"`, the ae ligature `"ae"`, thorn `"th"`),
and every character with no folding is dropped. Use it for destinations
that accept only ASCII (a file name, a fixed-width export, a legacy
system), so that a name written with an accented letter arrives with the
plain letter rather than as an encoding error or a question mark.

## Usage

``` r
toascii(x)
```

## Arguments

- x:

  A character vector (other types are converted with
  [`as.character()`](https://rdrr.io/r/base/character.html)).

## Value

A character vector the length of `x`, holding only printable ASCII
(bytes 32 to 126).

## Details

Bytes that are not valid UTF-8 at all (a mislabelled file) are dropped
first. The folding table covers the Latin-1 letters and the Latin
Extended-A letters of Central European names (U+0100 to U+017E); other
scripts are dropped rather than romanised. The result is then stripped
of anything outside the printable ASCII range (space to tilde), so
control characters such as line breaks and tabs are removed too. The
result never depends on the session's locale. A missing value becomes
`""`; other types are converted with
[`as.character()`](https://rdrr.io/r/base/character.html) first.

## See also

[`asciifb()`](https://rootcoder007.github.io/crimlab/reference/asciifb.md)
to transliterate only when needed,
[`wrtxtfb()`](https://rootcoder007.github.io/crimlab/reference/wrtxtfb.md)
to write text.

## Examples

``` r
toascii("Prof. \u00c1ngela Zorro Medina")
#> [1] "Prof. Angela Zorro Medina"
toascii(c("Se\u00e1n", "Zo\u00eb", "na\u00efve", "\u0141\u00f3d\u017a"))
#> [1] "Sean"  "Zoe"   "naive" "Lodz" 
all(charToRaw(toascii("caf\u00e9")) < 128)
#> [1] TRUE
toascii(c(NA, "two\nlines"))
#> [1] ""         "twolines"
```
