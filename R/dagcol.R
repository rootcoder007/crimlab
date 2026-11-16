#' Colliders of a directed acyclic graph
#'
#' `dagcol()` lists the colliders of a causal graph: every pair of causes
#' `a -> c <- b` of a common effect `c`. Conditioning on a collider (or on
#' a consequence of it), as a study of arrestees or of convicted people
#' does, opens a path between its causes and makes them associated when
#' they are not (Hernan and Robins 2020, ch. 8; see [colarr()]). A
#' collider whose causes are not adjacent is unshielded, a v-structure:
#' the v-structures and the skeleton of a graph are what the data can
#' tell apart (Verma and Pearl 1991).
#'
#' @details The graph is given by its edges, as the network functions take
#'   them ([ntdeg()]): `from[k] -> to[k]`. Repeated edges count once. The
#'   nodes are numbered in order of first appearance, each edge's `from`
#'   and then its `to`, and the colliders listed by `c` in that order and,
#'   within it, by the pair of causes in that order. The compiled kernel
#'   refuses a graph with a cycle (Kahn's algorithm), a self-loop
#'   included.
#'
#' @param from,to The labels of the cause and the effect of each edge:
#'   vectors of the same length (at least one edge) with no missing values.
#' @param triple `NULL` or three labels of nodes, `a`, `c` and `b`: whether
#'   `a -> c <- b` is a collider of the graph, and whether it is shielded.
#' @param x For `print()`, a `dagcol` result.
#' @inheritParams cemtch
#' @return An object of class `dagcol`, a list: `col`, a data frame with one
#'   row per collider (`a`, `c`, `b` and `shielded`, whether `a` and `b`
#'   are adjacent); `nodes`, the labels in order; `nedge`, the edges once
#'   repeats are merged; and `query`, `NULL` or for `triple` the logical
#'   values `collider` and `shielded`.
#' @references Pearl, J. (2009). *Causality: Models, Reasoning, and
#'   Inference* (2nd ed.). Cambridge University Press.
#'   \doi{10.1017/CBO9780511803161}
#'
#'   Verma, T. and Pearl, J. (1991). Equivalence and synthesis of causal
#'   models. *Proceedings of the Sixth Conference on Uncertainty in
#'   Artificial Intelligence*, 255-268.
#'
#'   Hernan, M. A. and Robins, J. M. (2020). *Causal Inference: What If*.
#'   Chapman and Hall/CRC.
#'   <https://miguelhernan.org/whatifbook>
#' @examples
#' # offending and mental illness both raise the chance of arrest; a study
#' # of arrestees conditions on the collider "arrest"
#' g <- dagcol(c("offend", "mental", "policing", "arrest"),
#'             c("arrest", "arrest", "arrest", "record"))
#' g
#' dagcol(c("a", "b", "a"), c("c", "c", "b"), triple = c("a", "c", "b"))$query
#' @export
dagcol <- function(from, to, triple = NULL) {
  f <- "dagcol"
  ok <- is.atomic(from) && is.atomic(to) && length(from) == length(to) &&
    length(from) > 0L
  if (!ok) {
    .ckstop(f, paste("`from` and `to` must be labels of the same length, one",
                     "edge each, at least one edge"))
  }
  if (anyNA(from) || anyNA(to)) {
    .ckstop(f, "`from` and `to` must have no missing values")
  }
  from <- as.character(from)
  to <- as.character(to)
  nodes <- unique(as.vector(rbind(from, to)))
  r <- .Call(C_clab_r_dagcol, match(from, nodes) - 1L, match(to, nodes) - 1L,
             length(nodes))
  col <- data.frame(a = nodes[r$a + 1L], c = nodes[r$c + 1L],
                    b = nodes[r$b + 1L], shielded = r$sh,
                    stringsAsFactors = FALSE)
  q <- NULL
  if (!is.null(triple)) {
    t <- as.character(triple)
    if (!is.atomic(triple) || length(t) != 3L || anyNA(triple) ||
          !all(t %in% nodes)) {
      .ckstop(f, "`triple` must be three labels of nodes: a, c, b")
    }
    hit <- col$c == t[2L] & ((col$a == t[1L] & col$b == t[3L]) |
                               (col$a == t[3L] & col$b == t[1L]))
    q <- c(collider = any(hit), shielded = any(col$shielded[hit]))
  }
  structure(list(col = col, nodes = nodes, nedge = r$nedge, query = q),
            class = "dagcol")
}

#' @rdname dagcol
#' @export
print.dagcol <- function(x, ...) {
  d <- x[["col"]]
  cat("dagcol: ", nrow(d), " colliders among ", length(x[["nodes"]]),
      " nodes, ", sum(!d$shielded), " unshielded (v-structures)\n", sep = "")
  if (nrow(d)) print(d, row.names = FALSE)
  q <- x[["query"]]
  if (!is.null(q)) {
    cat("triple: collider ", q[["collider"]], ", shielded ", q[["shielded"]],
        "\n", sep = "")
  }
  invisible(x)
}
