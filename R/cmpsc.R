#' Congruent matching profile segments (CMPS) of two bullet signatures
#'
#' How many pieces of one land engraved area's signature line up, at one
#' common shift, with another's: the CMPS score of Chen et al. (2019), the
#' count of congruent matching profile segments, with its share of the
#' segments. A high score points to the same barrel, a low one to
#' different barrels; no threshold is implied, since the scores of same-
#' and different-source pairs depend on the data and the settings (Ju and
#' Hofmann 2022).
#'
#' The comparison signature `a`, from its first value to its last, is cut
#' into consecutive basis segments of `len` points (a shorter remainder at
#' the end is dropped). For each, the
#' cross-correlation with the reference `b` is taken at every lag that
#' keeps the segment inside `b` (Pearson correlation over the pairs with
#' neither value missing), and the lags of its highest local maxima are
#' found: `npeaks[1]` of them at the segment's own length, then, for each
#' further element of `npeaks`, that many at twice the previous length
#' about the segment's centre (a length longer than `a` is left out). A
#' segment's candidate is the best peak of its longest length that lies
#' within `tx` points of a peak at every shorter length (the consistent
#' correlation peak; with one element of `npeaks`, all its peaks are
#' candidates). The congruent lag is the candidate that most segments have
#' a candidate within `tx` of (the middle one among ties); the segments
#' that do are the congruent matching profile segments, and their number
#' is the score. Signatures are taken as they come (detrended profiles,
#' as `bulletxtrctr` extracts them; Hare et al. 2017), in points of one
#' spacing; [chmb()] compares the raw profiles of tool marks instead.
#'
#' @param a The comparison signature: numbers, `NA` for a missing point.
#' @param b The reference signature, the same way.
#' @param len The length of a basis segment in points: one whole number
#'   from 2, at most the points of `a` from its first value to its last.
#' @param npeaks The number of peaks at each segment length: one to eight
#'   whole numbers from 1 (`c(5, 3, 1)`, Chen et al.'s choice).
#' @param tx The tolerance in points: one whole number from 0.
#' @return An object of class `cmpsc`, a list: `score`, the number of
#'   congruent matching profile segments; `nseg`, the number of basis
#'   segments (the highest score possible); `scaled`, `score / nseg`;
#'   `pos`, the congruent lag (`b[t + pos]` lines up with `a[t]`; `NA`
#'   when no segment has a candidate); and `segs`, a data frame with one
#'   row per basis segment: `seg`, `start` (its first point in `a`),
#'   `ccp` (its candidate lag nearest `pos`, `NA` for none) and `cmps`.
#' @references
#'   Chen, Z., Chu, W., Soons, J. A., Thompson, R. M., Song, J. and Zhao,
#'   X. (2019). Fired bullet signature correlation using the Congruent
#'   Matching Profile Segments (CMPS) method. *Forensic Science
#'   International*, 305, 109964. \doi{10.1016/j.forsciint.2019.109964}
#'
#'   Ju, W. and Hofmann, H. (2022). An open-source implementation of the
#'   CMPS algorithm for assessing similarity of bullets. *The R Journal*,
#'   14(2), 268-287. \doi{10.32614/RJ-2022-035}
#'
#'   Hare, E., Hofmann, H. and Carriquiry, A. (2017). Automatic matching of
#'   bullet land impressions. *The Annals of Applied Statistics*, 11(4),
#'   2332-2356. \doi{10.1214/17-AOAS1080}
#' @examples
#' set.seed(1)
#' a <- cumsum(rnorm(600, sd = 0.5)) + sin(seq_len(600) / 9)
#' b <- a[31:600] + rnorm(570, sd = 0.1)
#' cmpsc(a, b)
#' cmpsc(a, cumsum(rnorm(570, sd = 0.5)))$score
#' @export
cmpsc <- function(a, b, len = 50, npeaks = c(5, 3, 1), tx = 25) {
  .Call(C_clab_r_cmpsc, a, b, len, npeaks, tx)
}

#' @export
print.cmpsc <- function(x, ...) {
  cat(
    "cmpsc: ",
    x[["score"]],
    " of ",
    x[["nseg"]],
    " segments congruent (",
    format(x[["scaled"]]),
    ") at lag ",
    x[["pos"]],
    "\n",
    sep = ""
  )
  invisible(x)
}
