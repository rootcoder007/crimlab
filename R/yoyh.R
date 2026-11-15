# A change table rendered to a self-contained HTML page or to a PDF
# written by the package itself (no device, no external program), and
# the writer that picks the format from the file name. Colour encodes the
# DIRECTION of an improvement, not the sign of a change: with `better`
# set, a rise in days in segregation or in reported assaults is drawn as
# worse and a rise in a completion rate as better.

# the palettes are .YOY_PALETTES (R/plotx.R), shared with the plots

.ypick <- function(x, opt, nm, fn) {
  if (!is.character(x) || length(x) != 1L || !x %in% opt) {
    stop(fn, ": `", nm, "` must be one of ", toString(dQuote(opt, FALSE)),
         ".", call. = FALSE)
  }
  x
}

# the class of each row: bad, good or flat
.ycls <- function(chg, better) {
  out <- rep("flat", length(chg))
  if (better == "neutral") {
    return(out)
  }
  up <- !is.na(chg) & chg > 0
  dn <- !is.na(chg) & chg < 0
  out[up] <- if (better == "higher") "good" else "bad"
  out[dn] <- if (better == "higher") "bad" else "good"
  out
}

.yesc <- function(s) {
  s <- as.character(s)
  s[is.na(s)] <- ""
  s <- .yctl(s, ws = TRUE)
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  s <- gsub("\"", "&quot;", s, fixed = TRUE)
  gsub("'", "&#39;", s, fixed = TRUE)
}

.yfirst <- function(d, lag) {
  if (!("flag" %in% names(d)) ||
        !any(!is.na(d$flag) & d$flag == "no earlier period")) {
    return(NULL)
  }
  sprintf(
    "Note: the first %s%s %s nothing to compare with, so %s no change.",
    if (lag == 1L) "period" else sprintf("%d periods", lag),
    if ("group" %in% names(d)) " of each group" else "",
    if (lag == 1L) "has" else "have",
    if (lag == 1L) "it shows" else "they show"
  )
}

.ycss <- function(p) {
  paste0(
    ":root{--bad:", p$bad, ";--good:", p$good, ";--flat:", p$flat,
    ";--badbg:", p$badbg, ";--goodbg:", p$goodbg,
    ";--ink:#1a1a1a;--dim:#5c5c5c;--rule:#e0ddd8;--ground:#faf9f7;",
    "--panel:#ffffff}",
    "@media (prefers-color-scheme:dark){:root:not([data-theme=light]){",
    "--ink:#eceae6;--dim:#a6a29b;--rule:#332f2b;--ground:#171513;",
    "--panel:#201d1a;--badbg:#331f1a;--goodbg:#142a26}}",
    "*{box-sizing:border-box}",
    "body{margin:0;background:var(--ground);color:var(--ink);",
    "font:14px/1.5 ui-sans-serif,-apple-system,'Segoe UI',Roboto,",
    "'Helvetica Neue',Arial,sans-serif}",
    "main{max-width:1000px;margin:0 auto;padding:2.5rem 1.25rem 4rem}",
    "h1{font-size:1.6rem;margin:0 0 .25rem;letter-spacing:-.01em}",
    "h2{font-size:1rem;text-transform:uppercase;letter-spacing:.08em;",
    "color:var(--dim);margin:2.5rem 0 .75rem;font-weight:600}",
    "p.sub{color:var(--dim);margin:0 0 1.75rem}",
    ".scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:6px;",
    "background:var(--panel)}",
    "table{border-collapse:collapse;width:100%;min-width:520px}",
    "th{text-align:left;font-size:.72rem;text-transform:uppercase;",
    "letter-spacing:.06em;color:var(--dim);font-weight:600;",
    "padding:.65rem .75rem;border-bottom:1px solid var(--rule);",
    "white-space:nowrap}",
    "td{padding:.55rem .75rem;border-bottom:1px solid var(--rule);",
    "vertical-align:middle}",
    "tbody tr:last-child td{border-bottom:0}",
    "td.num,td.ci{text-align:right;font-variant-numeric:tabular-nums}",
    "td.per{font-variant-numeric:tabular-nums;white-space:nowrap}",
    "td.grp{font-weight:600}",
    "td.ci{color:var(--dim);font-size:.85rem;white-space:nowrap}",
    "td.pct{position:relative;text-align:right;white-space:nowrap;",
    "font-variant-numeric:tabular-nums}",
    "td.pct .arw{margin-right:.3rem;font-size:.8rem}",
    "td.pct.bad{color:var(--bad)}td.pct.good{color:var(--good)}",
    "td.pct.flat{color:var(--flat)}",
    ".bar{position:absolute;left:0;bottom:0;height:2px;opacity:.55}",
    ".bar.bad{background:var(--bad)}.bar.good{background:var(--good)}",
    ".bar.flat{background:var(--flat)}",
    "tr.bad td.pct{background:var(--badbg)}",
    "tr.good td.pct{background:var(--goodbg)}",
    "tr.note td{border-bottom:1px solid var(--rule);color:var(--dim);",
    "font-size:.8rem;padding:.2rem .75rem .5rem}",
    ".cards{display:grid;gap:.75rem;",
    "grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}",
    ".card{border:1px solid var(--rule);border-radius:6px;padding:.9rem;",
    "background:var(--panel)}",
    ".card .lab{font-size:.72rem;text-transform:uppercase;",
    "letter-spacing:.06em;color:var(--dim);font-weight:600}",
    ".card .big{font-size:1.5rem;margin:.2rem 0 .3rem;",
    "font-variant-numeric:tabular-nums}",
    ".card .cap{color:var(--dim);font-size:.8rem}",
    "ul.notes{color:var(--dim);font-size:.85rem;padding-left:1.1rem}",
    "p.gen{color:var(--dim);font-size:.75rem;margin-top:2.5rem}"
  )
}

.yargs <- function(title, sub, notes, fn) {
  .ytxt1(title, "title", fn)
  .ytxt1(sub, "sub", fn)
  if (!is.null(notes) && (!is.character(notes) || anyNA(notes))) {
    stop(fn, ": `notes` must be text (a character vector), or NULL.",
         call. = FALSE)
  }
}

#' Render a change table to HTML or PDF
#'
#' `yoyhtml()` writes one self-contained HTML page (every style inlined,
#' nothing fetched, so the file renders the same years later): the table,
#' a card per series with its change over the span ([yoysmry()]) and the
#' notes. `yoypdf()` writes a PDF of the table, paginated when it does not
#' fit, with the standard Helvetica font and no other program: crimlab
#' writes the file itself.
#'
#' Colour encodes the direction of an improvement, not the sign of the
#' change: with `better = "lower"` (fewer assaults, fewer days in
#' segregation) a rise is drawn in the colour for worse, with `"higher"`
#' (a completion rate) in the colour for better, and with `"neutral"` (the
#' default) every change in the neutral colour. The palette `"safe"`
#' (blue and orange) stays apart in the common forms of colour blindness;
#' `"mono"` is for print. The first period of each series, with nothing to
#' compare with, is said once under the table rather than on every first
#' row. The PDF holds ASCII only: dashes and arrows become `-` and `+`,
#' other characters their transliteration ([toascii()]).
#'
#' @param x A change table ([yoypx()]).
#' @param file The file path.
#' @param title The heading (and the page title).
#' @param sub The line under the heading.
#' @param pal The palette: `"diverging"`, `"safe"` or `"mono"`.
#' @param better Which direction is better: `"neutral"`, `"higher"` or
#'   `"lower"`.
#' @param digits The decimals of the change shown.
#' @param bars `TRUE` draws a bar under each change, proportional to its
#'   size.
#' @param intv `TRUE` shows the interval column of a count table.
#' @param notes `NULL`, or footnotes.
#' @param lag The periods per step of the table, for the summary and the
#'   note on the first periods.
#' @param width,height The PDF page in inches, each from 2 to 200.
#' @return `file`, invisibly.
#' @seealso [yoywr()], [yoycsv()].
#' @examples
#' # simulated police-reported counts by offence and year, any agency
#' d <- data.frame(year = rep(2021:2023, 2),
#'                 offence = rep(c("Assault", "Robbery"), each = 3),
#'                 n = c(812, 845, 790, 120, 131, 98))
#' y <- yoypx(d$n, d$year, d$offence)
#' f <- tempfile(fileext = ".html")
#' yoyhtml(y, f, title = "Reported offences (simulated)", better = "lower")
#' g <- tempfile(fileext = ".pdf")
#' yoypdf(y, g, title = "Reported offences (simulated)")
#' # simulated days in segregation by fiscal year
#' s <- yoylbl(yoypx(c(4100, 3900, 4400), 2021:2023), fylbl)
#' yoyhtml(s, f, title = "Days in segregation (simulated)", better = "lower")
#' unlink(c(f, g))
#' @name yoyhtml
#' @export
yoyhtml <- function(x, file, title = "Year-over-year change", sub = "",
                    pal = "diverging", better = "neutral", digits = 1L,
                    bars = TRUE, intv = TRUE, notes = NULL, lag = 1) {
  d <- .ytab(x, "yoyhtml")
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop("yoyhtml: `file` must be one file path.", call. = FALSE)
  }
  .yargs(title, sub, notes, "yoyhtml")
  p <- .YOY_PALETTES[[.ypick(pal, names(.YOY_PALETTES), "pal", "yoyhtml")]]
  better <- .ypick(better, c("neutral", "higher", "lower"), "better",
                   "yoyhtml")
  digits <- .ydig(digits, "yoyhtml")
  bars <- .yflag(bars, "bars", "yoyhtml")
  intv <- .yflag(intv, "intv", "yoyhtml")
  s <- .yshow(d, digits, intv, md = FALSE)
  chg <- d[[.ychg(d)]]
  cls <- .ycls(chg, better)
  sc <- suppressWarnings(max(abs(chg), na.rm = TRUE))
  if (!is.finite(sc) || sc == 0) sc <- 1
  hc <- names(s$cols)
  pc <- which(hc %in% c("change %", "points"))
  rows <- vapply(seq_len(nrow(d)), function(r) {
    cells <- vapply(seq_along(hc), function(j) {
      v <- s$cols[[j]][r]
      if (j == pc) {
        a <- if (is.na(chg[r])) "" else if (chg[r] > 0) "\u25b2" else
          if (chg[r] < 0) "\u25bc" else "\u2013"
        bar <- if (bars && !is.na(chg[r]) && is.finite(chg[r])) {
          sprintf("<span class=\"bar %s\" style=\"width:%.1f%%\"></span>",
                  cls[r], min(100, 100 * abs(chg[r]) / sc))
        } else {
          ""
        }
        return(paste0("<td class=\"pct ", cls[r], "\"><span class=\"arw\">",
                      a, "</span><span class=\"val\">", v, "</span>", bar,
                      "</td>"))
      }
      k <- switch(hc[j], group = "grp", prd = "per", interval = "ci", "num")
      paste0("<td class=\"", k, "\">", if (k %in% c("grp", "per")) .yesc(v)
             else v, "</td>")
    }, "")
    f <- s$flag[r]
    note <- if (is.na(f) || f == "no earlier period") "" else
      sprintf("<tr class=\"note\"><td colspan=\"%d\">note: %s</td></tr>",
              length(hc), .yesc(f))
    paste0("<tr class=\"", cls[r], "\">", paste(cells, collapse = ""),
           "</tr>", note)
  }, "")
  sm <- yoysmry(d, lag)
  cards <- vapply(seq_len(nrow(sm)), function(r) {
    paste0(
      "<div class=\"card\"><div class=\"lab\">",
      .yesc(if ("group" %in% names(sm)) sm$group[r] else "overall"),
      "</div><div class=\"big\">", .ypct(sm$totalpct[r], 1L, FALSE),
      "</div><div class=\"cap\">",
      .yesc(sprintf("%s to %s", sm$from[r], sm$to[r])),
      "</div><div class=\"cap\">",
      .yesc(sprintf("compound %s per period",
                    .ypct(sm$cagrpct[r], 1L, FALSE))),
      "</div><div class=\"cap\">",
      .yesc(sprintf("%d up, %d down, %d unchanged", sm$up[r], sm$down[r],
                    sm$flat[r])),
      "</div></div>")
  }, "")
  nt <- c(.yfirst(d, lag), notes)
  html <- c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    paste0("<title>", .yesc(title), "</title>"),
    paste0("<style>", .ycss(p), "</style>"),
    "</head><body><main>",
    paste0("<h1>", .yesc(title), "</h1>"),
    paste0("<p class=\"sub\">", .yesc(sub), "</p>"),
    "<div class=\"scroll\"><table><thead><tr>",
    paste(paste0("<th>", .yesc(hc), "</th>"), collapse = ""),
    "</tr></thead><tbody>",
    rows,
    "</tbody></table></div>",
    paste0("<h2>Across the span</h2><div class=\"cards\">",
           paste(cards, collapse = ""), "</div>"),
    if (length(nt)) {
      c("<ul class=\"notes\">", paste0("<li>", .yesc(nt), "</li>"), "</ul>")
    },
    paste0("<p class=\"gen\">",
           .yesc(sprintf("crimlab %s", utils::packageVersion("crimlab"))),
           "</p>"),
    "</main></body></html>"
  )
  writeBin(charToRaw(enc2utf8(paste0(paste(html, collapse = "\n"), "\n"))),
           file)
  invisible(file)
}

# the PDF's text: dashes and arrows as ASCII, then toascii()
.yasc <- function(x) {
  out <- as.character(x)
  for (p in list(c("\u2014", "-"), c("\u2013", "-"), c("\u25b2", "+"),
                 c("\u25bc", "-"), c("\u00b7", "|"), c("\u2500", "-"),
                 c("\u21b3", ">"))) {
    out <- gsub(p[1L], p[2L], out, fixed = TRUE)
  }
  toascii(out)
}

# Helvetica widths (1/1000 em) of the printable ASCII characters
.yhw <- c(278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278,
          333, 278, 278, rep(556, 10), 278, 278, 584, 584, 584, 556, 1015,
          667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833,
          722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611,
          278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556,
          556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278,
          556, 500, 722, 500, 500, 500, 334, 260, 334, 584)

.ywid <- function(s, size) {
  vapply(s, function(t) {
    k <- utf8ToInt(t) - 31L
    sum(.yhw[k[k >= 1L & k <= 95L]]) * size / 1000
  }, 0, USE.NAMES = FALSE)
}

.ypdfs <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("(", "\\(", s, fixed = TRUE)
  gsub(")", "\\)", s, fixed = TRUE)
}

.yrgb <- function(hex) {
  v <- strtoi(substring(hex, c(2L, 4L, 6L), c(3L, 5L, 7L)), 16L) / 255
  paste(formatC(v, format = "f", digits = 3), collapse = " ")
}

#' @rdname yoyhtml
#' @export
yoypdf <- function(x, file, title = "Year-over-year change", sub = "",
                   pal = "diverging", better = "neutral", digits = 1L,
                   intv = TRUE, notes = NULL, width = 11, height = 8.5,
                   lag = 1) {
  d <- .ytab(x, "yoypdf")
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop("yoypdf: `file` must be one file path.", call. = FALSE)
  }
  .yargs(title, sub, notes, "yoypdf")
  p <- .YOY_PALETTES[[.ypick(pal, names(.YOY_PALETTES), "pal", "yoypdf")]]
  better <- .ypick(better, c("neutral", "higher", "lower"), "better",
                   "yoypdf")
  digits <- .ydig(digits, "yoypdf")
  intv <- .yflag(intv, "intv", "yoypdf")
  for (nm in c("width", "height")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 2 ||
          v > 200) {
      stop("yoypdf: `", nm, "` must be one number from 2 to 200.",
           call. = FALSE)
    }
  }
  s <- .yshow(d, digits, intv, md = FALSE)
  hdr <- .yasc(names(s$cols))
  body <- lapply(s$cols, .yasc)
  cls <- .ycls(d[[.ychg(d)]], better)
  nt <- .yasc(c(.yfirst(d, lag), notes))
  W <- 72 * width
  H <- 72 * height
  m <- 0.6 * 72 / 6
  nc <- length(hdr)
  numfrom <- as.integer("group" %in% names(d)) + 2L
  w <- rep(1, nc)
  w[seq_len(numfrom - 1L)] <- 1.3
  xs <- m + (W - 2 * m) * cumsum(c(0, w / sum(w)))
  rowh <- 0.22 * 72
  per <- max(5L, as.integer(floor((0.82 * height - 0.22 * 1.1) / 0.22)))
  pages <- split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / per))
  if (!length(pages)) pages <- list(integer())
  np <- length(pages)
  txt <- function(x0, y0, s, size, font = "F1", rgb = "0 0 0", right = FALSE) {
    if (right) x0 <- x0 - .ywid(s, size)
    sprintf("BT /%s %g Tf %s rg %.2f %.2f Td (%s) Tj ET", font, size, rgb,
            x0, y0, .ypdfs(s))
  }
  dim <- .yrgb(p$flat)
  streams <- vapply(seq_len(np), function(k) {
    top <- H - m - 0.12 * H
    out <- c(txt(m, H - m - 18, .yasc(title), 18, "F2"),
             txt(m, H - m - 32, .yasc(sub), 9, rgb = dim),
             sprintf("%s RG 0.7 w %.2f %.2f m %.2f %.2f l S", dim, m,
                     top + 10, W - m, top + 10))
    for (j in seq_len(nc)) {
      r <- j >= numfrom && hdr[j] != "note"
      out <- c(out, txt(if (r) xs[j + 1L] - 3 else xs[j] + 3, top, hdr[j], 8,
                        "F2", dim, r))
    }
    out <- c(out, sprintf("%s RG 0.7 w %.2f %.2f m %.2f %.2f l S", dim, m,
                          top - 5, W - m, top - 5))
    y <- top - rowh * 1.1
    pc <- which(hdr %in% c("change %", "points"))
    for (i in pages[[k]]) {
      for (j in seq_len(nc)) {
        r <- j >= numfrom
        col <- if (j == pc) .yrgb(p[[cls[i]]]) else "0 0 0"
        out <- c(out, txt(if (r) xs[j + 1L] - 3 else xs[j] + 3, y,
                          body[[j]][i], 9, rgb = col, right = r))
      }
      y <- y - rowh
    }
    if (np > 1L) {
      out <- c(out, txt(W - m, m, sprintf("page %d of %d", k, np), 7,
                        rgb = dim, right = TRUE))
    }
    if (k == np && length(nt)) {
      yy <- m + 10 + 10 * length(nt)
      for (q in nt) {
        out <- c(out, txt(m, yy, q, 7.5, rgb = dim))
        yy <- yy - 10
      }
    }
    paste(out, collapse = "\n")
  }, "")
  .ypdfw(streams, W, H, .yasc(title), file)
}

# the PDF file: catalog, pages, two standard fonts, a page and a content
# stream per page, the cross-reference table
.ypdfw <- function(streams, W, H, title, file) {
  np <- length(streams)
  kids <- 5L + 2L * (seq_len(np) - 1L)
  obj <- c(
    "<< /Type /Catalog /Pages 2 0 R >>",
    sprintf("<< /Type /Pages /Kids [%s] /Count %d >>",
            paste(kids, "0 R", collapse = " "), np),
    paste("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica",
          "/Encoding /WinAnsiEncoding >>"),
    paste("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold",
          "/Encoding /WinAnsiEncoding >>")
  )
  for (k in seq_len(np)) {
    obj <- c(obj, sprintf(paste0(
      "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] ",
      "/Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> /Contents %d 0 R >>"
    ), W, H, kids[k] + 1L),
    sprintf("<< /Length %d >>\nstream\n%s\nendstream",
            nchar(streams[k], type = "bytes") + 0L, streams[k]))
  }
  obj <- c(obj, sprintf("<< /Title (%s) /Producer (crimlab) >>",
                        .ypdfs(title)))
  head <- "%PDF-1.4\n"
  parts <- sprintf("%d 0 obj\n%s\nendobj\n", seq_along(obj), obj)
  off <- nchar(head, type = "bytes") +
    cumsum(c(0, nchar(parts, type = "bytes")))[seq_along(parts)]
  xref <- paste0("xref\n0 ", length(obj) + 1L, "\n0000000000 65535 f \n",
                 paste(sprintf("%010d 00000 n \n", off), collapse = ""))
  at <- nchar(head, type = "bytes") + sum(nchar(parts, type = "bytes"))
  tr <- sprintf(paste0("trailer\n<< /Size %d /Root 1 0 R /Info %d 0 R >>\n",
                       "startxref\n%d\n%%%%EOF\n"),
                length(obj) + 1L, length(obj), at)
  writeBin(charToRaw(paste0(head, paste(parts, collapse = ""), xref, tr)),
           file)
  invisible(file)
}

#' Write a change table in the format its file name says
#'
#' Picks the writer from the extension of `file` (`.html` or `.htm`,
#' `.pdf`, `.csv`, `.tsv` or `.tab`, `.json`, `.md` or `.markdown`), or
#' from `format`, and passes the other arguments to it.
#'
#' @param x A change table ([yoypx()]).
#' @param file The file path.
#' @param format `"auto"` (from the extension), or one of `"html"`, `"pdf"`,
#'   `"csv"`, `"tsv"`, `"json"`, `"markdown"`.
#' @param ... The arguments of the chosen writer.
#' @return `file`, invisibly.
#' @seealso [yoycsv()], [yoyhtml()].
#' @examples
#' # simulated police-reported counts by offence and year
#' y <- yoypx(c(812, 845, 790), 2021:2023)
#' f <- tempfile(fileext = ".csv")
#' yoywr(y, f)
#' readLines(f)
#' unlink(f)
#' @export
yoywr <- function(x, file, format = "auto", ...) {
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop("yoywr: `file` must be one file path.", call. = FALSE)
  }
  format <- tolower(.ytxt1(format, "format", "yoywr"))
  if (format == "auto") {
    map <- c(html = "html", htm = "html", pdf = "pdf", csv = "csv",
             tsv = "tsv", tab = "tsv", json = "json", md = "markdown",
             markdown = "markdown")
    ext <- tolower(sub(".*\\.", "", basename(file)))
    if (!ext %in% names(map)) {
      stop("yoywr: cannot tell the format from \"", basename(file),
           "\"; pass `format`.", call. = FALSE)
    }
    format <- map[[ext]]
  }
  f <- switch(format, html = yoyhtml, pdf = yoypdf, csv = yoycsv,
              tsv = yoytsv, json = yoyjs, markdown = yoymd,
              stop("yoywr: unknown format: ", format, ".", call. = FALSE))
  f(x, file, ...)
  invisible(file)
}
