#' Weighted displacement quotient: did crime move next door?
#'
#' Asks whether an intervention in an area pushed crime into the area
#' around it (displacement) or reduced it there too (a diffusion of
#' benefits). Counts before and after the intervention are needed for the
#' treated area `a`, a buffer `b` around it and a control area `c` that
#' neither can reach (see [bufas()] for drawing them). With subscripts 0
#' before and 1 after, the measures of Bowers and Johnson (2003), as
#' Guerette (2009, Appendix B, equations 1 to 4) states them:
#'
#' * the success measure `sm = a1 / c1 - a0 / c0`, the change in the
#'   treated area relative to the control; negative when the intervention
#'   worked;
#' * the buffer displacement measure `bdm = b1 / c1 - b0 / c0`; positive
#'   when crime rose in the buffer (displacement), negative when it fell
#'   (diffusion);
#' * the weighted displacement quotient `wdq = bdm / sm`. With `sm`
#'   negative, `wdq` above 0 is a diffusion of benefits, below 0
#'   displacement, and below -1 displacement larger than the gain in the
#'   treated area;
#' * the gross effect `ge = a0 - a1`, the fall in the treated area, and the
#'   net effect `ne = a0 / c0 - a1 / c1` (that is `-sm`), the fall relative
#'   to the control;
#' * the total net effect `tne = (a0 c1 / c0 - a1) + (b0 c1 / c0 - b1)`,
#'   the crimes prevented in the treated area and the buffer together
#'   (positive when crime was prevented), with `ea = a0 c1 / c0` and `eb =
#'   b0 c1 / c0` the counts expected after had both areas followed the
#'   control. `ge` and `tne` count crimes, so they need counts, not rates.
#'
#' The quotient means little when `sm` is close to 0 (`NA` at exactly 0)
#' or the counts are small; the control must be like the other two areas
#' in trend, and far enough not to be reached by either. Guerette's worked
#' example (Appendix D: the treated area 207 before and 110 after, the
#' buffer 178 and 157, the control 308 and 318) gives `ge` 97, `ne` 0.326
#' and `tne` 130.5, as printed there; the guide prints the quotient as
#' 0.411, which its own formula and counts do not give (they give 0.258).
#'
#' @param a,b,c The counts before and after in the treated area, the
#'   buffer and the control: two numbers from 0 each, before first; those of
#'   `c` positive.
#' @return A named numeric vector: `sm`, `bdm`, `wdq`, `ge`, `ne`, `ea`,
#'   `eb` and `tne`.
#' @references
#'   Bowers, K. J. and Johnson, S. D. (2003). Measuring the geographical
#'   displacement and diffusion of benefit effects of crime prevention
#'   activity. *Journal of Quantitative Criminology*, 19(3), 275-301.
#'   \doi{10.1023/A:1024909009240}
#'
#'   Guerette, R. T. (2009). *Analyzing Crime Displacement and Diffusion*.
#'   Tool Guide No. 10. Center for Problem-Oriented Policing, Arizona State
#'   University.
#' @examples
#' # burglaries: treated 100 -> 60, buffer 50 -> 40, control 200 -> 180
#' wdqx(a = c(100, 60), b = c(50, 40), c = c(200, 180))
#' # Guerette (2009), Appendix D
#' wdqx(a = c(207, 110), b = c(178, 157), c = c(308, 318))
#' @export
wdqx <- function(a, b, c) {
  .Call(C_clab_r_wdqx, a, b, c)
}
