#' Transliterate text to plain ASCII
#'
#' Converts text to plain 7-bit ASCII: the common Latin accents and
#' ligatures are folded to their nearest ASCII letters (an accented capital
#' A becomes `"A"`, sharp s `"ss"`, the ae ligature `"ae"`, thorn `"th"`),
#' and every character with no folding is dropped. Use it for destinations
#' that accept only ASCII (a file name, a fixed-width export, a legacy
#' system), so that a name written with an accented letter arrives with the
#' plain letter rather than as an encoding error or a question mark.
#'
#' Bytes that are not valid UTF-8 at all (a mislabelled file) are dropped
#' first. The folding table covers the Latin-1 letters and the Latin
#' Extended-A letters of Central European names (U+0100 to U+017E); other
#' scripts are dropped rather than romanised. The result is then stripped
#' of anything outside the printable ASCII range (space to tilde), so
#' control characters such as line breaks and tabs are removed too. The
#' result never depends on the session's locale. A missing value becomes
#' `""`; other types are converted with [as.character()] first.
#'
#' @param x A character vector (other types are converted with
#'   [as.character()]).
#' @return A character vector the length of `x`, holding only printable
#'   ASCII (bytes 32 to 126).
#' @seealso [asciifb()] to transliterate only when needed, [wrtxtfb()] to
#'   write text.
#' @examples
#' toascii("Prof. \u00c1ngela Zorro Medina")
#' toascii(c("Se\u00e1n", "Zo\u00eb", "na\u00efve", "\u0141\u00f3d\u017a"))
#' all(charToRaw(toascii("caf\u00e9")) < 128)
#' toascii(c(NA, "two\nlines"))
#' @export
toascii <- function(x) {
  .Call(C_clab_r_ascii, .txtin(x, "toascii"))
}

#' Use text as it is, falling back to ASCII when it is not valid UTF-8
#'
#' Returns `x` unchanged when it is valid UTF-8 (so an accented name is
#' kept), and transliterates to plain ASCII with [toascii()] only the
#' elements that are not valid UTF-8 (an encoding error), or every element
#' when `force = TRUE` (for an ASCII-only destination). Names keep their
#' accents wherever UTF-8 is supported and degrade gracefully where it is
#' not.
#'
#' Each element is tested on its bytes with [validUTF8()], so the outcome
#' does not depend on the session locale. A missing value is kept as `NA`
#' unless `force = TRUE`, when it becomes `""`.
#'
#' @inheritParams toascii
#' @param force `TRUE` transliterates every element; `FALSE` (the default)
#'   only those that are not valid UTF-8. Anything but `TRUE` counts as
#'   `FALSE`.
#' @return A character vector the length of `x`.
#' @seealso [toascii()].
#' @examples
#' asciifb("\u00c1ngela")
#' asciifb("\u00c1ngela", force = TRUE)
#' bad <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
#' validUTF8(bad)
#' asciifb(bad)
#' @export
asciifb <- function(x, force = FALSE) {
  x <- .txtin(x, "asciifb")
  if (isTRUE(force)) {
    return(toascii(x))
  }
  bad <- !is.na(x) & !validUTF8(x)
  x[bad] <- toascii(x[bad])
  x
}

#' Write text as UTF-8, falling back to ASCII where it is not UTF-8
#'
#' Writes `text` to `path` as UTF-8, one element per line, for generated
#' text files (a README, a citation, a script) whose content may hold names
#' with accents. An element that is not valid UTF-8 is written as its ASCII
#' transliteration ([asciifb()]), so the file is always valid UTF-8.
#'
#' The lines are joined with `"\n"`, a final newline is added, and the
#' bytes are written through a binary connection, so the file is UTF-8
#' whatever the session locale and its line endings are `"\n"` on every
#' platform. An existing file is overwritten. A missing element is written
#' as the two letters `NA`, and an empty `text` writes a file holding a
#' single newline. A `path` that is not one non-empty string is refused
#' before anything is written; a directory that does not exist is an error
#' of the write.
#'
#' @param text The lines; [as.character()] converts any other type.
#' @param path The file, one non-empty string; its directory must exist.
#' @return `path`, invisibly.
#' @seealso [toascii()].
#' @examples
#' p <- wrtxtfb(c("line one", "line two"), tempfile(fileext = ".txt"))
#' readLines(p)
#' q <- wrtxtfb("Author: \u00c1ngela Zorro Medina", tempfile())
#' readBin(q, "raw", 20)
#' try(wrtxtfb("x", c("a.txt", "b.txt")))
#' @export
wrtxtfb <- function(text, path) {
  text <- asciifb(.txtin(text, "wrtxtfb", "text"))
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) {
    stop("wrtxtfb: `path` must be one non-empty file path.", call. = FALSE)
  }
  b <- charToRaw(enc2utf8(paste0(paste(text, collapse = "\n"), "\n")))
  writeBin(b, path)
  invisible(path)
}

# as.character() of an atomic vector; a function, list or other object is
# refused
.txtin <- function(x, fn, nm = "x") {
  if (!is.null(x) && !is.atomic(x)) {
    stop(fn, ": `", nm, "` must be text (or numbers, TRUE and FALSE).",
         call. = FALSE)
  }
  as.character(x)
}
