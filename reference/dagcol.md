# Colliders of a directed acyclic graph

`dagcol()` lists the colliders of a causal graph: every pair of causes
`a -> c <- b` of a common effect `c`. Conditioning on a collider (or on
a consequence of it), as a study of arrestees or of convicted people
does, opens a path between its causes and makes them associated when
they are not (Hernan and Robins 2020, ch. 8; see
[`colarr()`](https://rootcoder007.github.io/crimlab/reference/colarr.md)).
A collider whose causes are not adjacent is unshielded, a v-structure:
the v-structures and the skeleton of a graph are what the data can tell
apart (Verma and Pearl 1991).

## Usage

``` r
dagcol(from, to, triple = NULL)

# S3 method for class 'dagcol'
print(x, ...)
```

## Arguments

- from, to:

  The labels of the cause and the effect of each edge: vectors of the
  same length (at least one edge) with no missing values.

- triple:

  `NULL` or three labels of nodes, `a`, `c` and `b`: whether
  `a -> c <- b` is a collider of the graph, and whether it is shielded.

- x:

  For [`print()`](https://rdrr.io/r/base/print.html), a `dagcol` result.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `dagcol`, a list: `col`, a data frame with one row
per collider (`a`, `c`, `b` and `shielded`, whether `a` and `b` are
adjacent); `nodes`, the labels in order; `nedge`, the edges once repeats
are merged; and `query`, `NULL` or for `triple` the logical values
`collider` and `shielded`.

## Details

The graph is given by its edges, as the network functions take them
([`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md)):
`from[k] -> to[k]`. Repeated edges count once. The nodes are numbered in
order of first appearance, each edge's `from` and then its `to`, and the
colliders listed by `c` in that order and, within it, by the pair of
causes in that order. The compiled kernel refuses a graph with a cycle
(Kahn's algorithm), a self-loop included.

## References

Pearl, J. (2009). *Causality: Models, Reasoning, and Inference* (2nd
ed.). Cambridge University Press.
[doi:10.1017/CBO9780511803161](https://doi.org/10.1017/CBO9780511803161)

Verma, T. and Pearl, J. (1991). Equivalence and synthesis of causal
models. *Proceedings of the Sixth Conference on Uncertainty in
Artificial Intelligence*, 255-268.

Hernan, M. A. and Robins, J. M. (2020). *Causal Inference: What If*.
Chapman and Hall/CRC. <https://miguelhernan.org/whatifbook>

## Examples

``` r
# offending and mental illness both raise the chance of arrest; a study
# of arrestees conditions on the collider "arrest"
g <- dagcol(c("offend", "mental", "policing", "arrest"),
            c("arrest", "arrest", "arrest", "record"))
g
#> dagcol: 3 colliders among 5 nodes, 3 unshielded (v-structures)
#>       a      c        b shielded
#>  offend arrest   mental    FALSE
#>  offend arrest policing    FALSE
#>  mental arrest policing    FALSE
dagcol(c("a", "b", "a"), c("c", "c", "b"), triple = c("a", "c", "b"))$query
#> collider shielded 
#>     TRUE     TRUE 
```
