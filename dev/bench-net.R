# Speed of the network functions against igraph and sna, on fixed seeded
# graphs of 1e3, 1e4 and 1e5 nodes with heavy-tailed degrees (Chung-Lu:
# 3n edges, ends drawn with probability proportional to Pareto weights of
# exponent 2.5). Writes inst/bench/net.csv; nothing in it is typed by hand.
#   Rscript dev/bench-net.R [LIB]     (from the package root; LIB holds the
#                                      installed crimlab)
# Columns: seconds for crimlab (from labels: each call builds the graph),
# for the reference on a graph it has already built (ref_s) and for the
# reference including building its graph from the same labels (ref_build_s),
# each the median of three runs (one run when a run takes over 20 s; a
# fast call is repeated until a run lasts 0.2 s and divided);
# speed = ref_s / crimlab_s (above 1: crimlab is faster); faster is "yes"
# or "no"; note gives the reason a measure is not faster or was not run.
a <- commandArgs(TRUE)
if (length(a)) .libPaths(c(a[1], .libPaths()))
library(crimlab)
library(igraph)
has_sna <- requireNamespace("sna", quietly = TRUE)

tm <- function(expr) {
  e <- substitute(expr)
  p <- parent.frame()
  # repeated until a run lasts 0.2 s, so that fast calls are not lost in
  # the clock's resolution
  f <- function() {
    k <- 1
    repeat {
      t <- system.time(for (i in seq_len(k)) eval(e, p))[["elapsed"]]
      if (t >= 0.2) return(t / k)
      k <- k * 4
    }
  }
  t1 <- f()
  if (t1 > 20) return(t1)
  stats::median(c(t1, f(), f()))
}

graph <- function(n) {
  set.seed(n)
  w <- (1 - stats::runif(n))^(-1 / 1.5)
  m <- 3 * n
  f <- sample.int(n, m, replace = TRUE, prob = w)
  t <- sample.int(n, m, replace = TRUE, prob = w)
  lab <- paste0("n", seq_len(n))
  list(f = lab[f], t = lab[t], w = stats::runif(m, 0.5, 2), nodes = lab)
}

rows <- list()
add <- function(meas, ref, n, m, cs, rs, rbs, note = "") {
  sp <- if (is.na(rs) || is.na(cs)) NA else rs / cs
  fast <- if (is.na(sp)) NA else if (sp > 1) "yes" else "no"
  if (identical(fast, "no") && note == "") {
    note <- if (!is.na(rbs) && cs < rbs) {
      paste("slower than the measure alone; faster than building the",
            "graph from labels and measuring")
    } else {
      "not investigated"
    }
  }
  rows[[length(rows) + 1L]] <<- data.frame(
    measure = meas, reference = ref, nodes = n, edges = m,
    crimlab_s = cs, ref_s = rs, ref_build_s = rbs, speed = sp,
    faster = fast, note = note
  )
}

for (n in c(1e3, 1e4, 1e5)) {
  g <- graph(n)
  m <- length(g$f)
  bld <- function(dir = FALSE) {
    x <- graph_from_data_frame(data.frame(g$f, g$t, weight = g$w),
                               directed = dir, vertices = g$nodes)
    simplify(x, edge.attr.comb = list(weight = "sum"))
  }
  tb <- tm(bld())
  ig <- bld()
  igd <- bld(TRUE)
  tbd <- tm(bld(TRUE))
  cat("n =", n, "built\n")

  add("ntdeg", "igraph::degree + strength", n, m,
      tm(ntdeg(g$f, g$t, w = g$w, nodes = g$nodes)),
      tm({degree(ig); strength(ig)}), tb + tm({degree(ig); strength(ig)}))
  if (n <= 1e4) {
    r <- tm(betweenness(ig, weights = NA))
    add("ntbtw", "igraph::betweenness", n, m,
        tm(ntbtw(g$f, g$t, nodes = g$nodes)), r, tb + r)
    r <- tm(closeness(ig, weights = NA))
    add("ntcls", "igraph::closeness", n, m,
        tm(ntcls(g$f, g$t, nodes = g$nodes)), r, tb + r)
  } else {
    for (x in c("ntbtw", "ntcls")) {
      add(x, "igraph", n, m, NA, NA, NA,
          "not run: O(nodes x edges), about 3e10 steps for either package")
    }
  }
  r <- tm(eigen_centrality(ig))
  add("ntevc", "igraph::eigen_centrality", n, m,
      tm(ntevc(g$f, g$t, w = g$w, nodes = g$nodes)), r, tb + r)
  r <- tm(coreness(ig))
  add("ntcore", "igraph::coreness", n, m,
      tm(ntcore(g$f, g$t, nodes = g$nodes)), r, tb + r)
  r <- tm({components(igd, "weak"); components(igd, "strong")})
  add("ntcmp", "igraph::components weak + strong", n, m,
      tm(ntcmp(g$f, g$t, nodes = g$nodes, dir = TRUE)), r, tbd + r)
  r <- tm({set.seed(1); cluster_louvain(ig)})
  add("ntlvn", "igraph::cluster_louvain", n, m,
      tm(ntlvn(g$f, g$t, w = g$w, nodes = g$nodes, seed = 1)), r, tb + r)
  memb <- rep_len(1:20, n)
  r <- tm(modularity(ig, memb))
  add("ntmod", "igraph::modularity", n, m,
      tm(ntmod(g$f, g$t, memb, w = g$w, nodes = g$nodes)), r, tb + r)
  grp <- rep_len(c("a", "b", "c"), n)
  if (has_sna && n <= 1e3) {
    adj <- as.matrix(as_adjacency_matrix(igd, sparse = FALSE))
    adj[adj != 0] <- 1
    r <- tm(sna::brokerage(adj, grp))
    add("ntbrk", "sna::brokerage", n, m,
        tm(ntbrk(g$f, g$t, grp, nodes = g$nodes, dir = TRUE)), r, NA)
  } else {
    add("ntbrk", "sna::brokerage", n, m,
        tm(ntbrk(g$f, g$t, grp, nodes = g$nodes, dir = TRUE)), NA, NA,
        if (has_sna) "sna not run: it needs a dense nodes x nodes matrix"
        else "sna not installed")
  }
  # the largest component, for the Cheeger bounds
  big <- components(ig)
  keep <- names(big$membership)[big$membership == which.max(big$csize)]
  e <- g$f %in% keep & g$t %in% keep
  if (n <= 1e4) {
    cv <- ntchg(g$f[e], g$t[e], w = g$w[e], maxit = 2e4)$conv
    add("ntchg", "none", n, sum(e),
        tm(ntchg(g$f[e], g$t[e], w = g$w[e], maxit = 2e4)), NA, NA,
        paste("no counterpart in igraph or sna; converged within 2e4",
              "steps:", cv))
  }
  # co-offending: incidents of 1 to 4 offenders drawn by weight
  set.seed(n + 1)
  k <- sample(1:4, n, replace = TRUE)
  inc <- rep(seq_len(n), k)
  off <- sample(g$nodes, sum(k), replace = TRUE)
  bip <- function() {
    x <- graph_from_data_frame(data.frame(paste0("i", inc), off),
                               directed = FALSE)
    V(x)$type <- startsWith(V(x)$name, "n")
    bipartite_projection(x, which = "true", multiplicity = TRUE)
  }
  r <- tm(bip())
  add("coofg", "igraph::bipartite_projection (with building)", n, length(inc),
      tm(coofg(inc, off)), r, r)
  cat("n =", n, "done\n")
}

out <- do.call(rbind, rows)
out$platform <- paste(R.version$platform, "R", getRversion(), "igraph",
                      packageVersion("igraph"))
dir.create("inst/bench", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "inst/bench/net.csv", row.names = FALSE)
print(out[c("measure", "nodes", "crimlab_s", "ref_s", "ref_build_s", "speed",
            "faster")])
