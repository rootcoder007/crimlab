#' Co-offending graph from an incident table
#'
#' The one-mode projection of offenders onto incidents: two offenders are
#' joined when they are recorded in the same incident, and the weight of
#' the tie is the number of incidents they share. An offender recorded
#' twice in one incident counts once there. Offenders who never offend
#' with anyone appear in `nodes` but in no edge; give `nodes` to the
#' network functions to keep them as isolated nodes.
#'
#' Ties are listed in the order each pair is first met: incidents in the
#' order they first appear, and within an incident the pairs of its
#' offenders in the order of the rows; `from` is the offender who appears
#' first in `off`. Labels are compared as text (`as.character()`; a factor
#' by its labels).
#'
#' @param inc Incident identifiers, one per row: a vector of labels with no
#'   missing values.
#' @param off Offender identifiers, one per row: labels, same length.
#' @return A list: `edges`, a data frame with `from`, `to` (offender
#'   labels) and `w` (incidents shared), and `nodes`, every offender in
#'   order of first appearance.
#' @references Reiss, A. J. (1988). Co-offending and criminal careers.
#'   *Crime and Justice*, 10, 117-170. \doi{10.1086/449145}
#' @examples
#' g <- coofg(inc = c(1, 1, 1, 2, 2, 3), off = c("a", "b", "c", "a", "b", "d"))
#' g$edges
#' ntdeg(g$edges$from, g$edges$to, w = g$edges$w, nodes = g$nodes)
#' @export
coofg <- function(inc, off) {
  .Call(C_clab_r_coofg, inc, off)
}

#' Degree and strength of each node
#'
#' A graph is an edge list: `from[k]` is tied to `to[k]`. The rules below
#' hold for every network function in crimlab.
#'
#' * Nodes are identified by their labels, compared as text
#'   (`as.character()`; a factor by its labels, unused levels ignored).
#'   Results have one row per node, in the order nodes first appear: each
#'   edge's `from`, then its `to`, edge by edge. With `nodes`, the order is
#'   that of `nodes`, which must hold every label of the edges once; labels
#'   in `nodes` that no edge names are isolated nodes.
#' * A self-loop (`from[k] == to[k]`) is dropped; its node stays.
#' * Repeated edges (the same pair, in either direction when undirected)
#'   are merged into one edge whose weight is the sum of their weights (1
#'   each when `w` is `NULL`). So `deg` counts distinct neighbours and
#'   `wdeg` counts edges or their weight.
#' * A missing value in any argument is an error.
#'
#' @param from,to Labels of the two ends of each edge: vectors of the same
#'   length (zero when `nodes` gives every node) with no missing values.
#' @param w `NULL`, or one positive weight per edge (a tie's strength, such
#'   as the number of co-offences).
#' @param nodes `NULL`, or every node's label once: the node order of the
#'   results and the way to include isolated nodes.
#' @param dir `TRUE` for a directed graph (`from` sends a tie to `to`).
#' @return A data frame with one row per node: `node`, `deg` (distinct
#'   neighbours) and `wdeg` (summed weight); directed, `node`, `in`, `out`,
#'   `deg` (`in + out`), `win`, `wout` and `wdeg`. The degrees sum to twice
#'   the number of distinct edges: see `proof("ntdeg")`.
#' @references Wasserman, S. and Faust, K. (1994). *Social Network Analysis:
#'   Methods and Applications*. Cambridge University Press.
#'   \doi{10.1017/CBO9780511815478}
#' @examples
#' ntdeg(c("a", "a", "b", "a"), c("b", "c", "c", "b"))
#' ntdeg(c("a", "a"), c("b", "c"), nodes = c("c", "b", "a", "z"), dir = TRUE)
#' @export
ntdeg <- function(from, to, w = NULL, nodes = NULL, dir = FALSE) {
  .Call(C_clab_r_ntdeg, from, to, w, nodes, dir)
}

#' Betweenness of each node
#'
#' The sum, over pairs of other nodes, of the share of their shortest paths
#' that pass through the node (Freeman 1977), by Brandes' algorithm in
#' O(nodes x edges) time. Paths are counted in steps: weights are tie
#' strengths, not lengths, and are not taken. An undirected pair is counted
#' once. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntdeg
#' @param norm `TRUE` to divide by the number of pairs of other nodes,
#'   `(n - 1)(n - 2)` directed or `(n - 1)(n - 2) / 2` undirected, as
#'   `igraph::betweenness(normalized = TRUE)` does; `NaN` for fewer than 3
#'   nodes.
#' @return A data frame: `node`, `btw`.
#' @references Brandes, U. (2001). A faster algorithm for betweenness
#'   centrality. *Journal of Mathematical Sociology*, 25(2), 163-177.
#'   \doi{10.1080/0022250X.2001.9990249}
#'
#'   Freeman, L. C. (1977). A set of measures of centrality based on
#'   betweenness. *Sociometry*, 40(1), 35-41. \doi{10.2307/3033543}
#' @examples
#' ntbtw(c(1, 2, 3), c(2, 3, 4))
#' @export
ntbtw <- function(from, to, nodes = NULL, dir = FALSE, norm = FALSE) {
  .Call(C_clab_r_ntbtw, from, to, nodes, dir, norm)
}

#' Closeness of each node
#'
#' Distances are counted in steps (weights are not taken), along edges'
#' directions when `dir = TRUE` (out-distances). A disconnected graph is
#' handled in the two standard ways, both returned: `cls` and `ncls` count
#' only the nodes a node reaches (as `igraph::closeness()` does), so a node
#' close to the few it reaches scores high; `hrm`, the harmonic centrality,
#' sums `1 / distance` over every other node, an unreachable one adding 0,
#' so reaching more nodes always scores more. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntdeg
#' @return A data frame: `node`, `reach` (nodes reached), `cls` (1 / the
#'   summed distance to them), `ncls` (`reach` / that sum, the inverse mean
#'   distance; `cls` and `ncls` are `NaN` when nothing is reached) and `hrm`
#'   (the harmonic centrality).
#' @references Freeman, L. C. (1978). Centrality in social networks:
#'   conceptual clarification. *Social Networks*, 1(3), 215-239.
#'   \doi{10.1016/0378-8733(78)90021-7}
#'
#'   Marchiori, M. and Latora, V. (2000). Harmony in the small-world.
#'   *Physica A*, 285(3-4), 539-546. \doi{10.1016/S0378-4371(00)00311-3}
#' @examples
#' ntcls(c(1, 2, 3), c(2, 3, 4), nodes = 1:5)
#' @export
ntcls <- function(from, to, nodes = NULL, dir = FALSE) {
  .Call(C_clab_r_ntcls, from, to, nodes, dir)
}

#' Eigenvector centrality of each node
#'
#' The leading eigenvector of the (weighted) adjacency matrix of the
#' undirected graph, scaled to a maximum of 1 (Bonacich 1972). It is found
#' by power iteration on `A + c I`, `c` the largest weight, which has the
#' same eigenvectors and whose leading eigenvalue is strictly dominant even
#' on a bipartite graph (a co-offending network between two roles, a star),
#' where plain power iteration oscillates. The iteration starts from 1 for
#' every node and stops when no value changes by `tol` or more between two
#' steps, or after `maxit` steps (then `conv` is `FALSE`).
#'
#' On a disconnected graph the leading eigenvector lives on the components
#' with the largest leading eigenvalue; every other component gets 0 (each
#' component is iterated and scaled on its own, and a component whose
#' eigenvalue falls short of the largest by more than 1e-10 of it is set
#' to 0). When several components share the largest eigenvalue each is
#' scaled to a maximum of 1: any mix of them is an eigenvector too, so
#' compare nodes within a component. A graph without edges gives 1 for
#' every node. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntdeg
#' @param tol The stopping tolerance: one number between 0 and 1.
#' @param maxit The most steps: one whole number, at least 1.
#' @return A list: `cent`, a data frame with `node` and `evc`; `value`, the
#'   leading eigenvalue; `iter`, the steps taken; `conv`, whether the
#'   tolerance was met.
#' @references Bonacich, P. (1972). Factoring and weighting approaches to
#'   status scores and clique identification. *Journal of Mathematical
#'   Sociology*, 2(1), 113-120. \doi{10.1080/0022250X.1972.9989806}
#' @examples
#' ntevc(c("hub", "hub", "hub", "a"), c("a", "b", "c", "b"))
#' @export
ntevc <- function(from, to, w = NULL, nodes = NULL, tol = 1e-13,
                  maxit = 100000) {
  .Call(C_clab_r_ntevc, from, to, w, nodes, tol, maxit)
}

#' k-core of each node
#'
#' The coreness of a node is the largest k such that it belongs to a
#' subgraph in which every node has at least k neighbours (Seidman 1983),
#' computed in O(nodes + edges) by the algorithm of Batagelj and
#' Zaversnik. The cores are nested: the (k + 1)-core lies inside the
#' k-core (`proof("ntcore")`). Edges are taken as undirected, without
#' weights. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntdeg
#' @return A data frame: `node`, `core`.
#' @references Seidman, S. B. (1983). Network structure and minimum
#'   degree. *Social Networks*, 5(3), 269-287.
#'   \doi{10.1016/0378-8733(83)90028-X}
#'
#'   Batagelj, V. and Zaversnik, M. (2003). An O(m) algorithm for cores
#'   decomposition of networks. \doi{10.48550/arXiv.cs/0310049}
#' @examples
#' ntcore(c(1, 1, 1, 2, 2, 3, 4), c(2, 3, 4, 3, 4, 4, 5))
#' @export
ntcore <- function(from, to, nodes = NULL) {
  .Call(C_clab_r_ntcore, from, to, nodes)
}

#' Connected components
#'
#' Weak components join nodes linked by edges in either direction; strong
#' components (Tarjan's algorithm, without recursion) join nodes that
#' reach each other along the edges' directions. Undirected, the two
#' coincide. Components are numbered 1, 2, ... in the order of their first
#' node. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntdeg
#' @return A data frame: `node`, `weak`, `strong`.
#' @references Tarjan, R. (1972). Depth-first search and linear graph
#'   algorithms. *SIAM Journal on Computing*, 1(2), 146-160.
#'   \doi{10.1137/0201010}
#' @examples
#' ntcmp(c("a", "b", "c", "x"), c("b", "c", "a", "y"), dir = TRUE)
#' @export
ntcmp <- function(from, to, nodes = NULL, dir = FALSE) {
  .Call(C_clab_r_ntcmp, from, to, nodes, dir)
}

#' Communities by the Louvain method
#'
#' Groups nodes to raise the modularity (see [ntmod()]) of the undirected,
#' weighted graph (Blondel et al. 2008). Each level visits the nodes in
#' node order when `seed` is `NULL`, or in a random order drawn from `seed`
#' (one stream for the whole run), and moves each to the neighbouring
#' community with the largest gain in modularity, the first met among
#' equal gains, when that gain beats staying by more than 1e-12 of the
#' node's strength; passes repeat until none moves a node. The communities
#' then become the nodes of the next level, and a level that moves nothing
#' ends the run. The result is the same on every platform for the same
#' input and seed.
#'
#' The method is a heuristic: different orders reach different local
#' optima (on Zachary's karate club, modularity from about 0.39 to the
#' optimum 0.4198), so run several seeds and keep the best. Graph rules as
#' in [ntdeg()]; isolated nodes stay alone.
#'
#' @inheritParams ntdeg
#' @param res The resolution: one positive number; above 1 favours smaller
#'   communities.
#' @param seed `NULL` for node order, or one whole number from 0.
#' @return A list: `memb`, a data frame with `node` and `comm`
#'   (communities numbered 1, 2, ... in the order of their first node);
#'   `q`, its modularity; `levels`, the levels that moved a node.
#' @references Blondel, V. D., Guillaume, J.-L., Lambiotte, R. and
#'   Lefebvre, E. (2008). Fast unfolding of communities in large networks.
#'   *Journal of Statistical Mechanics*, P10008.
#'   \doi{10.1088/1742-5468/2008/10/P10008}
#' @examples
#' e <- c(1, 1, 2, 4, 4, 5, 3)
#' f <- c(2, 3, 3, 5, 6, 6, 4)
#' ntlvn(e, f)
#' @export
ntlvn <- function(from, to, w = NULL, nodes = NULL, res = 1, seed = NULL) {
  .Call(C_clab_r_ntlvn, from, to, w, nodes, res, seed)
}

#' Modularity of a partition
#'
#' `Q = sum over communities of in / 2m - res (tot / 2m)^2`, where `in` is
#' the weight of edges inside the community (counted from both ends),
#' `tot` the summed strength of its nodes and `2m` the total strength
#' (Newman and Girvan 2004), on the undirected, weighted graph. At `res =
#' 1` it lies between -1/2 and 1 (`proof("ntmod")`). `NaN` for a graph
#' without edges. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntdeg
#' @param memb The community of each node: labels, either named by node
#'   (every node named once; other names are ignored) or unnamed, one per
#'   node in the order the network functions report nodes.
#' @param res The resolution: one number, at least 0.
#' @return The modularity, one number.
#' @references Newman, M. E. J. and Girvan, M. (2004). Finding and
#'   evaluating community structure in networks. *Physical Review E*,
#'   69(2), 026113. \doi{10.1103/PhysRevE.69.026113}
#' @examples
#' ntmod(c(1, 1, 2, 4, 4, 5, 3), c(2, 3, 3, 5, 6, 6, 4),
#'       memb = c(1, 1, 1, 2, 2, 2))
#' @export
ntmod <- function(from, to, memb, w = NULL, nodes = NULL, res = 1) {
  .Call(C_clab_r_ntmod, from, to, memb, w, nodes, res)
}

#' Brokerage roles
#'
#' Counts, for each node v, the two-paths `i -> v -> j` with `i != j` and
#' no tie `i -> j` (v brokers between them), by the roles of Gould and
#' Fernandez (1989), given each node's group: `coord` (coordinator: i, v
#' and j in one group), `itin` (itinerant: i and j in one group, v in
#' another), `repr` (representative: v in i's group, j outside), `gate`
#' (gatekeeper: v in j's group, i outside) and `liais` (liaison: three
#' groups). Undirected, each tie runs both ways, so each unordered pair of
#' neighbours counts twice and `repr` equals `gate`.
#' Weights are not taken. Graph rules as in [ntdeg()]. Raw counts only:
#' the expected counts and z-scores of Gould and Fernandez are not given.
#'
#' @inheritParams ntdeg
#' @param grp The group of each node: labels, named by node or one per
#'   node in node order, as `memb` in [ntmod()].
#' @return A data frame: `node`, `grp` (as text), `coord`, `itin`, `repr`,
#'   `gate`, `liais` and `total`.
#' @references Gould, R. V. and Fernandez, R. M. (1989). Structures of
#'   mediation: a formal approach to brokerage in transaction networks.
#'   *Sociological Methodology*, 19, 89-126. \doi{10.2307/270949}
#' @examples
#' ntbrk(c("a", "b", "b"), c("b", "c", "d"),
#'       grp = c(a = "x", b = "y", c = "y", d = "z"), dir = TRUE)
#' @export
ntbrk <- function(from, to, grp, nodes = NULL, dir = FALSE) {
  .Call(C_clab_r_ntbrk, from, to, grp, nodes, dir)
}

#' Cheeger bounds on a network's bottleneck
#'
#' The conductance of a set of nodes S is the weight of the edges leaving
#' it over the smaller of its volume and that of the rest (volume: summed
#' strength); the graph's conductance h is the least over all sets, a
#' measure of how cleanly the network splits in two. Finding h is hard,
#' but the second-smallest eigenvalue `lambda2` of the normalised Laplacian
#' `I - D^-1/2 A D^-1/2` brackets it: `lambda2 / 2 <= h <= sqrt(2
#' lambda2)` (Cheeger's inequality; the lower half is `proof("ntchg")`).
#' The sweep over the eigenvector finds a set whose conductance (`sweep`)
#' is at most `sqrt(2 lambda2)`; h is at most `sweep`.
#'
#' `lambda2` is found by power iteration on `(I + D^-1/2 A D^-1/2) / 2`
#' kept orthogonal to `D^1/2 1`, from a fixed pseudo-random start, until
#' the residual is below `tol` or after `maxit` steps (`conv` `FALSE`).
#' Convergence is slow when `lambda2` and the next eigenvalue are close
#' (long paths, large sparse graphs). The graph is undirected and must be
#' connected; weights are tie strengths. Graph rules as in [ntdeg()].
#'
#' @inheritParams ntevc
#' @return A list: `bound`, a named vector `lambda2`, `lower`
#'   (`lambda2 / 2`), `upper` (`sqrt(2 lambda2)`) and `sweep`; `set`, the
#'   labels of the sweep set; `iter`; `conv`.
#' @references Chung, F. R. K. (1997). *Spectral Graph Theory*. CBMS
#'   Regional Conference Series in Mathematics 92. American Mathematical
#'   Society. \doi{10.1090/cbms/092}
#' @examples
#' # two triangles joined by one edge
#' ntchg(c(1, 1, 2, 4, 4, 5, 3), c(2, 3, 3, 5, 6, 6, 4))
#' @export
ntchg <- function(from, to, w = NULL, nodes = NULL, tol = 1e-10,
                  maxit = 100000) {
  .Call(C_clab_r_ntchg, from, to, w, nodes, tol, maxit)
}
