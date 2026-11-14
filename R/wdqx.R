#' Weighted displacement quotient: did crime move next door?
#'
#' Asks whether an intervention in an area pushed crime into the area
#' around it (displacement) or reduced it there too (a diffusion of
#' benefits). Counts before and after the intervention are needed for the
#' treated area `a`, a buffer `b` around it and a control area `c` that
#' neither can reach (see [bufas()] for drawing them). Following Bowers and
#' Johnson (2003), with subscripts 0 before and 1 after:
#'
#' * the success measure `sm = a1 / c1 - a0 / c0`, the change in the
#'   treated area relative to the control; negative when the intervention
#'   worked;
#' * the buffer displacement measure `bdm = b1 / c1 - b0 / c0`; positive
#'   when crime rose in the buffer (displacement), negative when it fell
#'   (diffusion);
#' * the quotient `wdq = bdm / sm`. With `sm` negative, `wdq` above 0 is a
#'   diffusion of benefits, below 0 displacement, and below -1 displacement
#'   larger than the gain in the treated area;
#' * the counts expected after in the treated area and the buffer had they
#'   followed the control (`ea = a0 c1 / c0`, `eb`), and `net`, the crimes
#'   prevented in both together, `(ea - a1) + (eb - b1)`.
#'
#' The quotient means little when `sm` is close to 0 (`NA` at exactly 0)
#' or the counts are small; the control must be like the other two areas
#' in trend, and far enough not to be reached by either.
#'
#' @param a,b,c The counts before and after in the treated area, the
#'   buffer and the control: two numbers from 0 each, before first; those of
#'   `c` positive.
#' @return A named numeric vector: `sm`, `bdm`, `wdq`, `ea`, `eb` and
#'   `net`.
#' @references Bowers, K. J. and Johnson, S. D. (2003). Measuring the
#'   geographical displacement and diffusion of benefit effects of crime
#'   prevention activity. *Journal of Quantitative Criminology*, 19(3),
#'   275-301. \doi{10.1023/A:1024909009240}
#' @examples
#' # burglaries: treated 100 -> 60, buffer 50 -> 40, control 200 -> 180
#' wdqx(a = c(100, 60), b = c(50, 40), c = c(200, 180))
#' @export
wdqx <- function(a, b, c) {
  .Call(C_clab_r_wdqx, a, b, c)
}
