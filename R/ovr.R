#' Text of a police oversight report
#'
#' Turns the page of a report published by a police oversight body (a
#' civilian investigation unit, an ombudsman, an inspector) into the plain
#' text its other tools read: markup removed, each paragraph, heading,
#' list item and table row on a line of its own, so that section headings
#' stand alone, and character entities decoded. Use it to read or search
#' saved reports without a browser, or before the rules that code them.
#'
#' Scripts, styles and comments are dropped. The start and end of a block
#' element (a paragraph, a heading, a division, a list item, a table row,
#' a line break and the like) end a line; every other tag is a space, and
#' the page's own line breaks are spaces, as a browser shows them. Named
#' entities common on report pages and every numeric entity are decoded
#' (curly quotes and dashes to their plain forms). Spaces are collapsed,
#' lines trimmed, and at most one empty line is kept between paragraphs.
#' Plain text (with `html = FALSE`, or with no markup when `html = NA`) is
#' only tidied that way, its line breaks kept. Nothing is fetched and
#' nothing is written; the pass is linear in the length of the page.
#'
#' @param x The pages or texts: a character vector. `NA` stays `NA`.
#' @param html `NA` (the default) to read each element as HTML when it
#'   holds markup (a `<` followed by a letter, `/` or `!`), `TRUE` or
#'   `FALSE` to say.
#' @return A character vector the length of `x`, in UTF-8.
#' @examples
#' ovrtxt("<h2>The Team</h2><p>Number of investigators assigned: 3</p>")
#' ovrtxt("<p>WO\n#1 Interviewed</p><p>WO&nbsp;#2 Interviewed</p>")
#' ovrtxt("a  b\r\nc", html = FALSE)
#' @export
ovrtxt <- function(x, html = NA) {
  f <- "ovrtxt"
  x <- .ovrin(x, "x", f)
  .Call(C_clab_r_ovrtxt, x, .ovrhtml(html, f))
}

# Text: a character vector (a vector of NA only is text that is missing),
# valid in its encoding, as UTF-8.
.ovrin <- function(x, what, f) {
  if (is.logical(x) && length(x) && all(is.na(x))) {
    x <- as.character(x)
  }
  if (!is.character(x) || !is.null(dim(x))) {
    .ckstop(f, "`%s` must be text (a character vector), not %s", what,
            class(x)[1L])
  }
  enc2utf8(.cktxt(x, what, f))
}

# NA, TRUE or FALSE.
.ovrhtml <- function(html, f) {
  if (!is.logical(html) || length(html) != 1L) {
    .ckstop(f, "`html` must be NA, TRUE or FALSE")
  }
  html
}

#' Dates written out in reports, as ISO dates
#'
#' Finds the first date written in each text, in English or French, and
#' gives it as `YYYY-MM-DD`, the form that sorts and that [as.Date()]
#' reads. `"January 5, 2023"`, `"Jan. 5th 2023"`, `"5 January 2023"`,
#' `"5th of January, 2023"`, `"5 janvier 2023"` and `"2023-01-05"` all
#' give `"2023-01-05"`.
#'
#' Month names are English or French, whole or abbreviated (with or
#' without a full stop), in any case, French accents optional; the day
#' may carry an ordinal (`"1st"`, `"1er"`). A day that does not exist
#' (`"February 30, 2023"`) is passed over, and the search goes on to the
#' next date in the text. The reading never depends on the session's
#' locale.
#'
#' @param x The texts: a character vector.
#' @return A character vector the length of `x`: the first date of each
#'   as `YYYY-MM-DD`, `NA` where there is none (or `x` is `NA`).
#' @examples
#' ovrdte(c("January 5, 2023", "le 1er février 2021", "no date"))
#' ovrdte("On February 30, 2023 or rather March 2, 2023, the SIU ...")
#' as.Date(ovrdte("Date: April 28, 2023"))
#' @export
ovrdte <- function(x) {
  .Call(C_clab_r_ovrdte, .ovrin(x, "x", "ovrdte"))
}

#' Field schemas for coding police oversight reports
#'
#' The fields a preset reads from each report, and the rule, cues and
#' section of each: a data frame that [ovrpar()] takes as `schema`. Copy
#' it and change rows, or write one of the same columns, to code the
#' reports of any oversight body; nothing in the rules is tied to one
#' agency.
#'
#' A field is read by one rule. `cue` holds the phrases the rule looks for,
#' alternatives separated by `"|"`, matched as written (case counts) and
#' as whole words at their ends. With `sect` (headings separated by
#' `"|"`), the rule reads only the section under the last line that is
#' one of those headings (a page that opens with a table of contents
#' repeats each heading), up to the next line that looks like a heading
#' (a short line starting with a capital, with no digit, not ending in a
#' stop or a comma); a report without the heading leaves the field
#' missing. The rules:
#'
#' * `after`: the text after a cue (past spaces and a colon), up to a
#'   comma, a semicolon, a full stop or the end of the line.
#' * `num`: the first number (digits, or a number word in English or
#'   French up to twenty) within 40 characters after a cue on its line.
#' * `date`: the first date ([ovrdte()]) in the sentence of a cue, or
#'   with no cue the first date in the section.
#' * `tag`: the highest number of a numbered tag (`"WO #3"`, `"WO 3"`,
#'   `"AT no 3"`): the count of the people the tags number.
#' * `sent`: the sentence that holds the first cue.
#' * `prev`: the line before the last line that is a cue (a signature:
#'   the name above a title).
#' * `lines`: every line that starts with a cue, distinct, joined by
#'   `"; "`.
#' * `name`: a name of up to four capitalised words that ends in a cue
#'   (`"Barrie Police Service"` for the cue `"Police Service"`), without
#'   a leading "The".
#' * `flag`: `FALSE` when a phrase of `alt` occurs, else `TRUE` when a
#'   cue does.
#' * `before`: the number just before a cue (`"34"` in `"34-year-old"`).
#' * `word`: the word just after a cue (`"man"` in `"34-year-old man"`).
#'
#' The preset `"siu"` reads the director's reports of Ontario's Special
#' Investigations Unit as published in English since 2020: the case
#' number, the police service, the dates of the incident, of the
#' notification and of the decision, the investigators assigned, the
#' witness officials, civilian witnesses and subject officials, the age
#' and sex of the affected person, whether the director found grounds
#' for charges, the director's name, the municipality, the injuries and
#' the provisions considered.
#'
#' @param preset The name of a preset: `"siu"`.
#' @return A data frame with one row per field: `name`, `rule`, `cue`,
#'   `alt`, `sect` and `desc` (what the field records).
#' @examples
#' s <- ovrsch()
#' s[, c("name", "rule")]
#' # a field of one's own: the date a complaint came in
#' rbind(s, data.frame(name = "date_received", rule = "date",
#'                     cue = "received the complaint", alt = "", sect = "",
#'                     desc = "the date the complaint came in"))
#' @export
ovrsch <- function(preset = "siu") {
  f <- "ovrsch"
  preset <- .ckarg(preset, "siu", "preset", f)
  as.data.frame(.Call(C_clab_r_ovrsch, preset), stringsAsFactors = FALSE)
}
