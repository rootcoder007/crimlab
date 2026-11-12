# Reference values for the network functions, from igraph, sna and base R
# (eigen() of the dense normalised Laplacian), on fixed graphs. Writes
# tests/testthat/netref.rds; test-net.R compares crimlab with it.
#   Rscript dev/netref.R        (from the package root; needs igraph, sna)
# Nothing here calls crimlab: the references are made without it.
library(igraph)
set.seed(20261112)

# A graph by the rules of crimlab: loops dropped, repeated edges merged
# with summed weights, nodes in the order given.
mk <- function(f, t, w, nodes, dir) {
  g <- graph_from_data_frame(data.frame(f = f, t = t, weight = w),
                             directed = dir, vertices = nodes)
  simplify(g, edge.attr.comb = list(weight = "sum"))
}

one <- function(f, t, w = NULL, nodes, dir = FALSE, grp = NULL) {
  ww <- if (is.null(w)) rep(1, length(f)) else w
  g <- mk(f, t, ww, nodes, dir)
  v <- V(g)$name
  nmd <- function(x) stats::setNames(as.numeric(x), v)
  r <- list(edges = data.frame(from = f, to = t), nodes = nodes, dir = dir)
  r$edges$w <- w
  if (dir) {
    r$indeg <- nmd(degree(g, mode = "in"))
    r$outdeg <- nmd(degree(g, mode = "out"))
    r$instr <- nmd(strength(g, mode = "in"))
    r$outstr <- nmd(strength(g, mode = "out"))
  } else {
    r$deg <- nmd(degree(g))
    r$str <- nmd(strength(g))
  }
  r$btw <- nmd(betweenness(g, weights = NA))
  r$nbtw <- nmd(betweenness(g, weights = NA, normalized = TRUE))
  r$cls <- nmd(closeness(g, weights = NA))
  r$ncls <- nmd(closeness(g, weights = NA, normalized = TRUE))
  r$hrm <- nmd(harmonic_centrality(g, weights = NA))
  r$weak <- nmd(components(g, "weak")$membership)
  r$strong <- nmd(components(g, "strong")$membership)
  if (!is.null(grp)) {
    a <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
    a[a != 0] <- 1
    k <- sna::brokerage(a, grp)$raw.nli
    rownames(k) <- v
    r$grp <- stats::setNames(grp, v)
    r$brk <- k
  }
  if (!dir) {
    e <- eigen_centrality(g)
    r$evc <- nmd(e$vector)
    r$evcval <- e$value
    r$core <- nmd(coreness(g))
    memb <- sample(4, length(v), replace = TRUE)
    r$memb <- stats::setNames(memb, v)
    r$mod <- modularity(g, memb)
    r$mod5 <- modularity(g, memb, resolution = 0.5)
    lv <- lapply(1:20, function(s) {
      set.seed(s)
      cluster_louvain(g)
    })
    r$louv <- vapply(lv, modularity, 0)
    best <- lv[[which.max(r$louv)]]
    r$louvmemb <- nmd(membership(best))
    if (is_connected(g) && length(v) <= 400) {
      a <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
      s <- 1 / sqrt(rowSums(a))
      l <- diag(length(v)) - a * outer(s, s)
      r$lam2 <- sort(eigen(l, symmetric = TRUE, only.values = TRUE)$values)[2]
    }
  }
  r
}

ref <- list()

# Zachary's karate club, labels "1" .. "34"
z <- as_edgelist(make_graph("Zachary"))
ref$karate <- one(as.character(z[, 1]), as.character(z[, 2]),
                  nodes = as.character(1:34),
                  grp = rep(c("a", "b", "c"), length.out = 34))

# weighted, undirected, with repeats, loops and isolated nodes, nodes given
# in a shuffled order
n <- 60
lab <- paste0("v", 1:n)
m <- 200
f <- sample(n - 3, m, replace = TRUE)
t <- sample(n - 3, m, replace = TRUE)
f[1:5] <- t[1:5]
f <- c(f, t[6:15])
t <- c(t, f[6:15])
ref$wrand <- one(lab[f], lab[t], w = round(stats::runif(length(f), 0.5, 3), 2),
                 nodes = sample(lab), grp = sample(c("x", "y", "z"), n, TRUE))

# directed, weighted, with repeats and loops
n <- 50
lab <- paste0("d", 1:n)
f <- sample(n, 180, replace = TRUE)
t <- sample(n, 180, replace = TRUE)
f[1:4] <- t[1:4]
f <- c(f, f[10:20])
t <- c(t, t[10:20])
ref$drand <- one(lab[f], lab[t], w = stats::runif(length(f), 1, 2),
                 nodes = lab, dir = TRUE,
                 grp = sample(c("x", "y", "z"), n, TRUE))

# bipartite: two roles, ties only across them, connected
b1 <- paste0("o", 1:20)
b2 <- paste0("s", 1:15)
f <- c(b1, sample(b1, 40, replace = TRUE))
t <- c(sample(b2, 20, replace = TRUE), sample(b2, 40, replace = TRUE))
f <- c(f, rep("o1", 15))
t <- c(t, b2)
ref$bip <- one(f, t, nodes = unique(c(rbind(f, t))))

# two K6 joined by one edge: one best partition
k <- t(utils::combn(6, 2))
f <- c(k[, 1], k[, 1] + 6, 1)
t <- c(k[, 2], k[, 2] + 6, 7)
ref$barbell <- one(as.character(f), as.character(t),
                   nodes = unique(as.character(c(rbind(f, t)))))

# preferential attachment, a realistic degree distribution
g <- sample_pa(300, m = 2, directed = FALSE)
e <- as_edgelist(g)
ref$pa <- one(paste0("n", e[, 1]), paste0("n", e[, 2]),
              nodes = paste0("n", 1:300))

attr(ref, "made") <- c(igraph = as.character(packageVersion("igraph")),
                       sna = as.character(packageVersion("sna")),
                       R = as.character(getRversion()))
saveRDS(ref, "tests/testthat/netref.rds", version = 2)
