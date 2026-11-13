#' DBSCAN clusters of incidents
#'
#' Density-based clusters (Ester, Kriegel, Sander and Xu 1996): a core
#' point has at least `minpts` points within distance `eps`, itself
#' included; core points within `eps` of each other share a cluster, and a
#' point that is not core but lies within `eps` of a core point (a border
#' point) joins the first cluster that reaches it, clusters being grown
#' from the core points in input order. The rest is noise. With `border =
#' FALSE` border points are noise too (DBSCAN*, Campello, Moulavi and
#' Sander 2013), which makes the result independent of the order. These
#' are the clusters of `dbscan::dbscan()`, numbered the same way.
#'
#' @inheritParams cgstat
#' @param eps The neighbourhood radius: positive, in the unit of the
#'   coordinates (metres when `lonlat = TRUE`).
#' @param minpts The fewest points (the point itself included) within `eps`
#'   of a core point: a whole number from 1.
#' @param border `TRUE` to put border points in clusters.
#' @return A list of class `dbscn`: `cl` (each point's cluster from 1, `NA`
#'   for noise), `core` (whether each point is a core point) and `sizes`
#'   (the points in each cluster).
#' @references Ester, M., Kriegel, H.-P., Sander, J. and Xu, X. (1996). A
#'   density-based algorithm for discovering clusters in large spatial
#'   databases with noise. *Proceedings of the Second International
#'   Conference on Knowledge Discovery and Data Mining*, 226-231.
#'
#'   Campello, R. J. G. B., Moulavi, D. and Sander, J. (2013). Density-based
#'   clustering based on hierarchical density estimates. *Advances in
#'   Knowledge Discovery and Data Mining*, 160-172.
#'   \doi{10.1007/978-3-642-37456-2_14}
#' @examples
#' set.seed(1)
#' p <- rbind(cbind(x = rnorm(30, 0, 0.2), y = rnorm(30, 0, 0.2)),
#'            cbind(x = rnorm(30, 3, 0.2), y = rnorm(30, 1, 0.2)),
#'            cbind(x = runif(10, -2, 5), y = runif(10, -2, 3)))
#' d <- dbscn(p, eps = 0.3, minpts = 5)
#' d
#' table(d$cl, useNA = "ifany")
#' @export
dbscn <- function(x, y = NULL, eps, minpts = 5, border = TRUE,
                  lonlat = FALSE) {
  .Call(C_clab_r_dbscn, x, y, eps, minpts, border, lonlat)
}

#' @export
print.dbscn <- function(x, ...) {
  cat("dbscn: ", length(x[["sizes"]]), " clusters of ",
      length(x[["cl"]]), " points, ", sum(is.na(x[["cl"]])), " noise\n",
      sep = "")
  if (length(x[["sizes"]])) {
    print(data.frame(cluster = seq_along(x[["sizes"]]),
                     points = x[["sizes"]]), row.names = FALSE)
  }
  invisible(x)
}
