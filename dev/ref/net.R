# Reference values of tests/testthat/test-net.R, "measures agree with
# igraph on new random graphs": igraph's strength, betweenness, closeness
# and harmonic centrality on the graphs the test draws. From the package
# root, with igraph in the library:
#   Rscript dev/ref/net.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- "igraph"
source(file.path("dev", "ref", "put.R"))

set.seed(11)
out <- list()
for (k in 1:5) {
  n <- 30 + 10 * k
  m <- 3 * n
  f <- sample(n, m, replace = TRUE)
  t <- sample(n, m, replace = TRUE)
  w <- stats::runif(m, 0.5, 2)
  lab <- paste0("p", seq_len(n + 2))
  for (dir in c(FALSE, TRUE)) {
    g <- igraph::graph_from_data_frame(
      data.frame(f = lab[f], t = lab[t], weight = w), directed = dir,
      vertices = lab
    )
    g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"))
    out[[sprintf("%d%s", k, if (dir) "d" else "u")]] <- list(
      inp = sum(f, t, w),
      wdeg = unname(igraph::strength(g, mode = "all")),
      btw = unname(igraph::betweenness(g, weights = NA)),
      cls = unname(igraph::closeness(g, weights = NA)),
      hrm = unname(igraph::harmonic_centrality(g, weights = NA))
    )
  }
}
put(list(rand = out), "net", need)
