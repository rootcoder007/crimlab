#' Hypercube queue of patrol units
#'
#' Larson's hypercube model of patrol cars: calls arrive from each atom
#' (beat, reporting area) at rate `lam` and go to the first free unit in
#' the atom's order of preference; each unit serves at rate `mu`. The
#' state is the set of busy units, \eqn{2^n} states for `n` units, and
#' the stationary distribution is found exactly. From it come each unit's
#' workload (the share of time busy), the share of calls each unit answers
#' from each atom, the probability that every unit is busy (a lost call,
#' or with `queue = "wait"` a call that waits; Larson's infinite-line
#' model, the queue lumped into the all-busy state), the mean queue and
#' wait, and with travel times the mean response travel time. Whatever the
#' orders, the number of busy units follows Erlang's loss or delay system.
#' The noisy-OR chance that patrols detect an event is [nordet()].
#'
#' @param lam The call rate of each atom: positive numbers.
#' @param pref The dispatch orders: a matrix with one row per atom, each a
#'   permutation of the unit numbers 1 to `n` (at most 10 units), or a
#'   list of such vectors.
#' @param mu The service rate of a unit (1 / mean service time): one
#'   positive number.
#' @param queue `"loss"` (a call finding every unit busy is lost) or
#'   `"wait"` (it waits; the calls must be fewer than the units serve).
#' @param tt `NULL` or the travel times: a matrix with one row per unit
#'   and one column per atom.
#' @return An object of class `hypq`, a list: `state`, each set of busy
#'   units (`busy`, one character per unit, `"1"` busy) with its
#'   probability `p`; `wl`, the workloads; `disp`, the share of all calls
#'   answered at once by each unit (rows) from each atom (columns);
#'   `queue`; `pfull`, the probability that every unit is busy; `lq` and
#'   `wq`, the mean queue and wait (0 for `"loss"`); `rho`, the calls over
#'   what the units can serve; `travel`, `NULL` or the mean travel time
#'   overall (`mean`) and by atom (`atom`).
#' @references
#'   Larson, R. C. (1974). A hypercube queuing model for facility location
#'   and redistricting in urban emergency services. *Computers and
#'   Operations Research*, 1(1), 67-95.
#'   \doi{10.1016/0305-0548(74)90076-8}
#' @examples
#' pref <- rbind(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1), c(3, 1, 2))
#' tt <- matrix(c(2, 5, 8, 6, 3, 7, 9, 6, 2, 7, 8, 3), 3)
#' hypq(c(0.4, 0.7, 0.3, 0.6), pref, mu = 1.5, tt = tt)
#' @export
hypq <- function(lam, pref, mu = 1, queue = "loss", tt = NULL) {
  .Call(C_clab_r_hypq, lam, pref, mu, queue, tt)
}

#' @export
print.hypq <- function(x, ...) {
  cat("hypq: ", length(x[["wl"]]), " units, rho ", format(x[["rho"]]),
      "; P(all busy) ", format(x[["pfull"]]),
      if (x[["queue"]] == "wait") {
        paste0(", mean queue ", format(x[["lq"]]))
      },
      "\nworkloads: ", paste(format(x[["wl"]]), collapse = " "), "\n",
      sep = "")
  if (!is.null(x[["travel"]])) {
    cat("mean travel time ", format(x[["travel"]][["mean"]]), "\n", sep = "")
  }
  invisible(x)
}

#' Square-root law of response distance
#'
#' Kolesar and Blum's square-root law: with available units spread like a
#' random (Poisson) pattern over a region of area `area`, the mean
#' distance from a call to the nearest available unit is \eqn{c
#' \sqrt{A / (u (1 - b))}} for `units` u of which a share `busy` b is busy;
#' \eqn{c = 1/2} in straight lines and \eqn{\sqrt{2 \pi} / 4 \approx 0.627}
#' on a rectangular street grid. Quadrupling the units halves the distance
#' and scaling area and units together leaves it unchanged
#' (`proof("sqrtlw")`). The busy share can come from [hypq()]'s workloads.
#'
#' @param area The area served: one positive number (squared distance
#'   units).
#' @param units The units deployed: one positive number.
#' @param busy The share of units busy: one number from 0, below 1.
#' @param metric `"euclid"` (straight lines) or `"grid"` (right-angle
#'   streets).
#' @param speed `NULL` or the travel speed (distance per unit time): one
#'   positive number.
#' @return A list: `dist`, the mean distance; `const`, the constant `c`;
#'   `avail`, the units available; `time`, the mean travel time (`NA`
#'   without `speed`).
#' @references
#'   Kolesar, P. and Blum, E. H. (1973). Square root laws for fire engine
#'   response distances. *Management Science*, 19(12), 1368-1378.
#'   \doi{10.1287/mnsc.19.12.1368}
#' @examples
#' sqrtlw(100, 4)
#' sqrtlw(100, 8, busy = 0.3, metric = "grid", speed = 0.5)
#' @export
sqrtlw <- function(area, units, busy = 0, metric = "euclid", speed = NULL) {
  .Call(C_clab_r_sqrtlw, area, units, busy, metric, speed)
}
