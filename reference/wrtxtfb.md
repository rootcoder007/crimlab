# Write text as UTF-8, falling back to ASCII where it is not UTF-8

Writes `text` to `path` as UTF-8, one element per line, for generated
text files (a README, a citation, a script) whose content may hold names
with accents. An element that is not valid UTF-8 is written as its ASCII
transliteration
([`asciifb()`](https://rootcoder007.github.io/crimlab/reference/asciifb.md)),
so the file is always valid UTF-8.

## Usage

``` r
wrtxtfb(text, path)
```

## Arguments

- text:

  The lines; [`as.character()`](https://rdrr.io/r/base/character.html)
  converts any other type.

- path:

  The file, one non-empty string; its directory must exist. A name with
  no suffix is given `.txt`.

## Value

The path written, invisibly.

## Details

The lines are joined with `"\n"`, a final newline is added, and the
bytes are written through a binary connection, so the file is UTF-8
whatever the session locale and its line endings are `"\n"` on every
platform. An existing file is overwritten. A missing element is written
as the two letters `NA`, and an empty `text` writes a file holding a
single newline. A `path` that is not one non-empty string is refused
before anything is written; a directory that does not exist is an error
of the write.

## See also

[`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md).

## Examples

``` r
p <- wrtxtfb(c("line one", "line two"), tempfile(fileext = ".txt"))
readLines(p)
#> [1] "line one" "line two"
q <- wrtxtfb("Author: \u00c1ngela Zorro Medina", tempfile())
readBin(q, "raw", 20)
#>  [1] 41 75 74 68 6f 72 3a 20 c3 81 6e 67 65 6c 61 20 5a 6f 72 72
try(wrtxtfb("x", c("a.txt", "b.txt")))
#> Error : wrtxtfb: `path` must be one non-empty file path.
```
