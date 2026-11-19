#' Credit for pre-sentence custody
#'
#' The credit a Canadian court gives for time held in custody before
#' sentence: each day held counts for at most one day of the sentence
#' (Criminal Code s. 719(3)) or, where the circumstances justify it, for
#' at most one and a half days (s. 719(3.1)). With the sentence imposed,
#' the time left to serve after the credit is given too.
#'
#' The credit is `days * ratio`, in days, unrounded: courts state the
#' credit they grant (s. 719(3.3)) and round it in their own way. A ratio
#' of 1.5 is the enhanced credit that R v. Summers (2014) holds is
#' generally available to an offender denied bail, because time held
#' before sentence earns no remission and no parole eligibility; a ratio
#' below 1 is a court's choice to give less than the maximum. Any further
#' mitigation for harsh conditions of custody is a matter of the sentence,
#' not of this arithmetic. The time left is the sentence less the credit,
#' and nothing when the credit is the larger.
#'
#' @param days The days held in custody before sentence: numbers from 0.
#' @param ratio The credit for each day: one number, or one per value of
#'   `days`, from 0 to 1.5.
#' @param sentence `NULL`, or the sentence imposed, in days: one number, or
#'   one per value of `days`, from 0.
#' @return A data frame with one row per value of `days`: `days`, `ratio`,
#'   `credit` and, with `sentence`, `sentence` and `remain` (the days left
#'   to serve).
#' @references
#'   Criminal Code, R.S.C. 1985, c. C-46, s. 719(3), (3.1) and (3.3).
#'
#'   R. v. Summers, 2014 SCC 26, \[2014\] 1 S.C.R. 575.
#' @examples
#' pscred(c(10, 30, 45))
#' pscred(c(10, 30, 45), ratio = c(1, 1.5, 1.5), sentence = 60)
#' @export
pscred <- function(days, ratio = 1.5, sentence = NULL) {
  f <- "pscred"
  days <- as.vector(.cknum(days, "days", f, scalar = FALSE, min = 0))
  ratio <- as.vector(.cknum(ratio, "ratio", f, scalar = FALSE, min = 0,
                            max = 1.5))
  n <- length(days)
  .pscone(ratio, n, "ratio", f)
  if (!is.null(sentence)) {
    sentence <- as.vector(.cknum(sentence, "sentence", f, scalar = FALSE,
                                 min = 0))
    .pscone(sentence, n, "sentence", f)
  }
  k <- .Call(C_clab_r_pscred, days, ratio, sentence)
  out <- data.frame(days = days, ratio = rep_len(ratio, n),
                    credit = k$credit)
  if (!is.null(sentence)) {
    out$sentence <- rep_len(sentence, n)
    out$remain <- k$remain
  }
  out
}

#' One value, or one per value of `days`.
#' @noRd
.pscone <- function(x, n, what, f) {
  if (length(x) != 1L && length(x) != n) {
    .ckstop(f, "`%s` must have one value or one per value of `days`", what)
  }
  invisible(x)
}
