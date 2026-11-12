# Values of the network kernels are checked in C++ (src/t_net.cpp) against
# closed forms and brute force. Here: the rules of the R interface (labels,
# node order, repeated edges, loops, isolated nodes, missing values), the
# messages, and agreement with igraph and sna: stored references written
# by dev/netref.R, and live comparisons when the packages are installed.

ref <- readRDS(test_path("netref.rds"))

# a column of a result, as a vector named by node
byn <- function(d, col) stats::setNames(d[[col]], d[["node"]])

# the same partition, whatever the numbering
samep <- function(a, b) {
  length(a) == length(b) &&
    nrow(unique(data.frame(a, b))) == length(unique(a)) &&
    length(unique(a)) == length(unique(b))
}

test_that("measures agree with igraph and sna (stored references)", {
  for (nm in names(ref)) {
    r <- ref[[nm]]
    e <- r$edges
    nd <- r$nodes
    v <- nd
    d <- ntdeg(e$from, e$to, w = e$w, nodes = nd, dir = r$dir)
    expect_identical(d$node, v, info = nm)
    if (r$dir) {
      expect_equal(byn(d, "in"), r$indeg[v], tolerance = 0, info = nm)
      expect_equal(byn(d, "out"), r$outdeg[v], tolerance = 0, info = nm)
      expect_equal(byn(d, "win"), r$instr[v], tolerance = 1e-12, info = nm)
      expect_equal(byn(d, "wout"), r$outstr[v], tolerance = 1e-12, info = nm)
    } else {
      expect_equal(byn(d, "deg"), r$deg[v], tolerance = 0, info = nm)
      expect_equal(byn(d, "wdeg"), r$str[v], tolerance = 1e-12, info = nm)
    }
    b <- ntbtw(e$from, e$to, nodes = nd, dir = r$dir)
    expect_equal(byn(b, "btw"), r$btw[v], tolerance = 1e-10, info = nm)
    b <- ntbtw(e$from, e$to, nodes = nd, dir = r$dir, norm = TRUE)
    expect_equal(byn(b, "btw"), r$nbtw[v], tolerance = 1e-10, info = nm)
    cl <- ntcls(e$from, e$to, nodes = nd, dir = r$dir)
    expect_equal(byn(cl, "cls"), r$cls[v], tolerance = 1e-10, info = nm)
    expect_equal(byn(cl, "ncls"), r$ncls[v], tolerance = 1e-10, info = nm)
    expect_equal(byn(cl, "hrm"), r$hrm[v], tolerance = 1e-10, info = nm)
    cm <- ntcmp(e$from, e$to, nodes = nd, dir = r$dir)
    expect_true(samep(cm$weak, r$weak[v]), info = nm)
    expect_true(samep(cm$strong, r$strong[v]), info = nm)
    if (!is.null(r$grp)) {
      k <- ntbrk(e$from, e$to, grp = r$grp, nodes = nd, dir = r$dir)
      expect_equal(unname(as.matrix(k[3:8])), unname(r$brk[v, ]),
                   tolerance = 0, info = nm)
    }
    if (r$dir) next
    ev <- ntevc(e$from, e$to, w = e$w, nodes = nd)
    expect_true(ev$conv, info = nm)
    expect_equal(byn(ev$cent, "evc"), r$evc[v], tolerance = 1e-9, info = nm)
    expect_equal(ev$value, r$evcval, tolerance = 1e-10, info = nm)
    expect_equal(byn(ntcore(e$from, e$to, nodes = nd), "core"), r$core[v],
                 tolerance = 0, info = nm)
    expect_equal(ntmod(e$from, e$to, r$memb, w = e$w, nodes = nd), r$mod,
                 tolerance = 1e-12, info = nm)
    expect_equal(ntmod(e$from, e$to, r$memb, w = e$w, nodes = nd, res = 0.5),
                 r$mod5, tolerance = 1e-12, info = nm)
    # Louvain: the best of several seeds against igraph's best, and never
    # below igraph's worst by more than a little
    q <- vapply(0:9, function(s) {
      ntlvn(e$from, e$to, w = e$w, nodes = nd, seed = s)$q
    }, 0)
    expect_gte(max(q), max(r$louv) - 0.01)
    if (!is.null(r$lam2)) {
      ch <- ntchg(e$from, e$to, w = e$w, nodes = nd)
      expect_equal(ch$bound[["lambda2"]], r$lam2, tolerance = 1e-10, info = nm)
    }
  }
})

test_that("a graph with one best partition gives igraph's partition", {
  r <- ref[["barbell"]]
  l <- ntlvn(r$edges$from, r$edges$to)
  expect_true(samep(l$memb$comm, r$louvmemb[l$memb$node]))
  expect_equal(l$q, max(r$louv), tolerance = 1e-12)
})

test_that("measures agree with igraph on new random graphs", {
  skip_if_not_installed("igraph")
  set.seed(11)
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
      d <- ntdeg(lab[f], lab[t], w = w, nodes = lab, dir = dir)
      expect_equal(d$wdeg, unname(igraph::strength(g, mode = "all")),
                   tolerance = 1e-12)
      expect_equal(ntbtw(lab[f], lab[t], nodes = lab, dir = dir)$btw,
                   unname(igraph::betweenness(g, weights = NA)),
                   tolerance = 1e-10)
      cl <- ntcls(lab[f], lab[t], nodes = lab, dir = dir)
      expect_equal(cl$cls, unname(igraph::closeness(g, weights = NA)),
                   tolerance = 1e-10)
      expect_equal(cl$hrm,
                   unname(igraph::harmonic_centrality(g, weights = NA)),
                   tolerance = 1e-10)
    }
  }
})

test_that("results are by label, in order of first appearance", {
  d <- ntdeg(c("z", "a", "m"), c("a", "y", "z"))
  expect_identical(d$node, c("z", "a", "y", "m"))
  # a factor is read by its labels, not its codes; unused levels are not
  # nodes
  f <- factor(c("z", "a"), levels = c("a", "q", "z"))
  expect_identical(ntdeg(f, c("a", "b"))$node, c("z", "a", "b"))
  # numbers are labels too
  expect_identical(ntdeg(c(10, 2), c(2, 7))$node, c("10", "2", "7"))
  # nodes sets the order and adds isolated nodes
  d <- ntdeg("a", "b", nodes = c("q", "b", "a"))
  expect_identical(d$node, c("q", "b", "a"))
  expect_identical(d$deg, c(0L, 1L, 1L))
  # nodes may carry names; they are dropped
  expect_null(names(ntdeg("a", "b", nodes = c(x = "a", y = "b"))$node))
  # every node isolated
  expect_identical(ntcore(character(), character(), nodes = 1:2)$core,
                   c(0L, 0L))
})

test_that("repeated edges merge, self-loops drop", {
  d <- ntdeg(c("a", "b", "a", "a"), c("b", "a", "a", "c"), w = c(1, 2, 5, 1))
  expect_identical(d$node, c("a", "b", "c"))
  expect_identical(d$deg, c(2L, 1L, 1L))
  expect_identical(d$wdeg, c(4, 3, 1))
  # directed, a -> b and b -> a stay apart
  d <- ntdeg(c("a", "b", "a"), c("b", "a", "b"), dir = TRUE)
  expect_identical(d$out, c(1L, 1L))
  expect_identical(d$wout, c(2, 1))
  # a loop's node stays
  expect_identical(ntdeg("a", "a")$deg, 0L)
  expect_identical(ntbtw(c("a", "a", "b"), c("b", "b", "c"))$btw, c(0, 1, 0))
})

test_that("the co-offending projection counts shared incidents", {
  g <- coofg(c("i1", "i1", "i1", "i2", "i2", "i2", "i3", "i4"),
             c("b", "a", "c", "a", "b", "a", "d", "c"))
  expect_identical(g$nodes, c("b", "a", "c", "d"))
  expect_identical(g$edges$from, c("b", "b", "a"))
  expect_identical(g$edges$to, c("a", "c", "c"))
  expect_identical(g$edges$w, c(2, 1, 1))
  # factors by their labels
  g2 <- coofg(factor(c(1, 1)), factor(c("x", "y"), levels = c("y", "x")))
  expect_identical(g2$edges$from, "x")
  # nobody offends together: no edges, every offender a node
  g3 <- coofg(1:3, c("a", "b", "a"))
  expect_identical(nrow(g3$edges), 0L)
  expect_identical(g3$nodes, c("a", "b"))
})

test_that("per-node values match by name or by position", {
  e <- c("a", "b", "c")
  f <- c("b", "c", "d")
  m1 <- ntmod(e, f, c(d = 2, c = 2, b = 1, a = 1, zz = 9))
  m2 <- ntmod(e, f, c(1, 1, 2, 2))
  expect_identical(m1, m2)
  k <- ntbrk(e, f, grp = factor(c("u", "u", "v", "v")))
  expect_identical(k$grp, c("u", "u", "v", "v"))
})

#' @srrstats {G5.4c} src/t_net.cpp checks published values on Zachary's
#'   karate club: the betweenness of member 1 and the optimum modularity
#'   (Brandes et al. 2008); netref.rds holds igraph's and sna's values.
#' @srrstats {G5.9b} ntlvn() with the same seed gives the same partition;
#'   other seeds may stop at other local optima, never below the
#'   modularity of the singletons (src/t_net.cpp).
test_that("communities and seeds", {
  e <- c(1, 1, 2, 4, 4, 5, 3)
  f <- c(2, 3, 3, 5, 6, 6, 4)
  a <- ntlvn(e, f, seed = 3)
  expect_identical(a, ntlvn(e, f, seed = 3))
  expect_identical(a$memb$comm, rep(1:2, each = 3))
  expect_equal(a$q, ntmod(e, f, a$memb$comm), tolerance = 1e-15)
  expect_true(is.nan(ntlvn("a", "a")$q))
  expect_true(is.nan(ntmod("a", "a", 1)))
})

test_that("the Cheeger bounds bracket the sweep and the conductance", {
  e <- c(1, 1, 2, 4, 4, 5, 3)
  f <- c(2, 3, 3, 5, 6, 6, 4)
  ch <- ntchg(e, f)
  b <- ch$bound
  expect_lte(b[["lower"]], b[["sweep"]])
  expect_lte(b[["sweep"]], b[["upper"]])
  expect_identical(sort(ch$set), c("4", "5", "6"))
  # the cut of one edge over the volume of a triangle with its bridge
  expect_equal(b[["sweep"]], 1 / 7, tolerance = 1e-15)
  # lambda2 against the dense eigenvalues of the normalised Laplacian
  a <- matrix(0, 6, 6)
  a[cbind(e, f)] <- 1
  a <- a + t(a)
  dg <- 1 / sqrt(rowSums(a))
  l <- diag(6) - a * outer(dg, dg)
  expect_equal(b[["lambda2"]], sort(eigen(l, symmetric = TRUE)$values)[2],
               tolerance = 1e-10)
  expect_false(ntchg(e, f, maxit = 2)$conv)
})

test_that("eigenvector centrality on bipartite and disconnected graphs", {
  s <- ntevc(rep("c", 4), c("a", "b", "d", "e"))
  expect_equal(s$cent$evc, c(1, rep(0.5, 4)), tolerance = 1e-12)
  expect_equal(s$value, 2, tolerance = 1e-12)
  d <- ntevc(c("a", "a", "b", "x"), c("b", "c", "c", "y"), nodes = c("a", "b",
    "c", "x", "y", "z"))
  expect_identical(d$cent$evc, c(1, 1, 1, 0, 0, 0))
  expect_identical(ntevc("a", "a")$cent$evc, 1)
  expect_false(ntevc(1:3, 2:4, maxit = 1)$conv)
})

test_that("network functions stop on bad input, naming the argument", {
  errs("ntdeg", list(
    list("`from`", c("a", NA), c("b", "c")),
    list("`to`", "a", NA),
    list("`from`", list("a"), "b"),
    list("`from`", Sys.Date(), "b"),
    list("`from`", matrix("a", 2, 2), rep("b", 4)),
    list("length of `from`", "a", c("b", "c")),
    list("`to` names a node", "a", "b", nodes = "a"),
    list("`from` names a node", "a", "b", nodes = "b"),
    list("must not repeat", "a", "b", nodes = c("a", "b", "a")),
    list("`nodes`", "a", "b", nodes = c("a", NA)),
    list("no nodes", character(), character()),
    list("`w`", "a", "b", w = 0),
    list("`w`", "a", "b", w = c(1, 2)),
    list("`dir`", "a", "b", dir = NA),
    list("`dir`", "a", "b", dir = c(TRUE, FALSE))
  ))
  errs("ntbtw", list(list("`norm`", "a", "b", norm = 1)))
  errs("ntcls", list(list("`dir`", "a", "b", dir = "yes")))
  errs("ntevc", list(
    list("`tol`", "a", "b", tol = 0),
    list("`maxit`", "a", "b", maxit = 0.5)
  ))
  errs("ntcore", list(list("`from`", NULL, "b")))
  errs("ntcmp", list(list("`to`", "a", NULL)))
  errs("ntlvn", list(
    list("`res`", "a", "b", res = 0),
    list("`seed`", "a", "b", seed = -1)
  ))
  errs("ntmod", list(
    list("one value per node", "a", "b", 1),
    list("no value for node b", "a", "b", c(a = 1)),
    list("must not repeat a name", "a", "b", c(a = 1, b = 1, a = 2)),
    list("`memb`", "a", "b", c(1, NA)),
    list("`res`", "a", "b", c(1, 1), res = -1)
  ))
  errs("ntbrk", list(list("`grp`", "a", "b", NULL)))
  errs("ntchg", list(
    list("connected", c("a", "c"), c("b", "d")),
    list("two nodes", character(), character(), nodes = "a"),
    list("`tol`", "a", "b", tol = 1),
    list("`maxit`", "a", "b", maxit = -1)
  ))
  errs("coofg", list(
    list("length of `inc`", 1, 1:2),
    list("at least one incident", character(), character()),
    list("`off`", 1, NA)
  ))
})

test_that("network results have the documented shape", {
  e <- c("a", "b")
  f <- c("b", "c")
  shp <- function(x) paste(names(x), vapply(x, typeof, ""), collapse = " ")
  expect_identical(shp(ntdeg(e, f)), "node character deg integer wdeg double")
  expect_identical(shp(ntdeg(e, f, dir = TRUE)), paste(
    "node character in integer out integer deg integer win double",
    "wout double wdeg double"
  ))
  expect_identical(shp(ntbtw(e, f)), "node character btw double")
  expect_identical(shp(ntcls(e, f)), paste(
    "node character reach integer cls double ncls double",
    "hrm double"
  ))
  expect_identical(shp(ntevc(e, f)),
                   "cent list value double iter integer conv logical")
  expect_identical(shp(ntcore(e, f)), "node character core integer")
  expect_identical(shp(ntcmp(e, f)),
                   "node character weak integer strong integer")
  expect_identical(shp(ntlvn(e, f)), "memb list q double levels integer")
  expect_identical(shp(ntbrk(e, f, 1:3)), paste(
    "node character grp character coord double itin double repr double",
    "gate double liais double total double"
  ))
  expect_identical(shp(ntchg(e, f)),
                   "bound double set character iter integer conv logical")
  expect_identical(shp(coofg(1, 1)), "edges list nodes character")
  expect_s3_class(ntdeg(e, f), "data.frame")
})
