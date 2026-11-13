#' Hot spots and how long they last
#'
#' Counts incidents in square grid cells of side `cell` (from the corner
#' `org`) for each period and marks as hot, in each period, the `k =
#' ceiling(top * m)` busiest of the `m` cells that ever hold an incident,
#' with every cell tied with the k-th; when a tie would make more than `2
#' k` cells hot the counts do not single them out and the tied cells are
#' left out (counted in `tiedout`). Each period gets its concentration (the
#' fewest cells holding half its incidents), the share of its incidents in
#' hot cells, and the Jaccard overlap of its hot set with the period
#' before (`NA` when either is empty); each cell gets the number of
#' periods it was hot and its longest run of successive hot periods, the
#' persistence that tells chronic hot spots from flare-ups (Weisburd,
#' Bushway, Lum and Yang 2004).
#'
#' @inheritParams cgstat
#' @param prd The period of each incident: labels (a year, a month, a
#'   quarter), taken in sorted order (a factor in the order of its levels),
#'   with no missing values.
#' @param cell The side of a cell: positive, in the unit of the coordinates
#'   (metres when `lonlat = TRUE`).
#' @param top The share of the active cells counted hot, above 0 and at
#'   most 1.
#' @param org `NULL` for the lowest `x` and `y`, or the grid's corner
#'   `c(x, y)` (projected, as the coordinates are, in metres when `lonlat =
#'   TRUE`).
#' @return A list of class `hsper`: `periods`, a data frame with one row per
#'   period: `period`, `n` (incidents), `active` (cells with any), `nhalf`
#'   and `shalf` (the fewest cells holding half the incidents, and that
#'   number over the active cells of all periods), `hot`, `tiedin` (hot
#'   cells beyond `k` from a tie), `tiedout`, `hotshr` (the share of
#'   incidents in hot cells) and `jacc`; `cells`, a data frame with one row
#'   per active cell: its centre `x`, `y` (in degrees when `lonlat =
#'   TRUE`), `n`, `nhot` and `run`; `hot`, a logical matrix of cells by
#'   periods; `cell`, each incident's cell (a row of `cells`); and `k`.
#' @references Weisburd, D., Bushway, S., Lum, C. and Yang, S.-M. (2004).
#'   Trajectories of crime at places: a longitudinal study of street
#'   segments in the city of Seattle. *Criminology*, 42(2), 283-322.
#'   \doi{10.1111/j.1745-9125.2004.tb00521.x}
#' @examples
#' set.seed(1)
#' n <- 600
#' p <- cbind(x = c(rnorm(300, 2, 0.3), runif(300, 0, 10)),
#'            y = c(rnorm(300, 2, 0.3), runif(300, 0, 10)))
#' h <- hsper(p, prd = rep(2021:2024, length.out = n), cell = 1, top = 0.1)
#' h
#' subset(h$cells, nhot == 4)
#' @export
hsper <- function(x, y = NULL, prd, cell, top = 0.05, org = NULL,
                  lonlat = FALSE) {
  .Call(C_clab_r_hsper, x, y, prd, cell, top, org, lonlat)
}

#' @export
print.hsper <- function(x, ...) {
  cat("hsper: ", nrow(x[["cells"]]), " active cells, the ", x[["k"]],
      " busiest hot in each of ", nrow(x[["periods"]]), " periods\n",
      sep = "")
  print(x[["periods"]], row.names = FALSE)
  invisible(x)
}
