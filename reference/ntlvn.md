# Communities by the Louvain method

Groups nodes to raise the modularity (see
[`ntmod()`](https://rootcoder007.github.io/crimlab/reference/ntmod.md))
of the undirected, weighted graph (Blondel et al. 2008). Each level
visits the nodes in node order when `seed` is `NULL`, or in a random
order drawn from `seed` (one stream for the whole run), and moves each
to the neighbouring community with the largest gain in modularity, the
first met among equal gains, when that gain beats staying by more than
1e-12 of the node's strength; passes repeat until none moves a node. The
communities then become the nodes of the next level, and a level that
moves nothing ends the run. The result is the same on every platform for
the same input and seed.

## Usage

``` r
ntlvn(from, to, w = NULL, nodes = NULL, res = 1, seed = NULL)
```

## Arguments

- from, to:

  Labels of the two ends of each edge: vectors of the same length (zero
  when `nodes` gives every node) with no missing values.

- w:

  `NULL`, or one positive weight per edge (a tie's strength, such as the
  number of co-offences).

- nodes:

  `NULL`, or every node's label once: the node order of the results and
  the way to include isolated nodes.

- res:

  The resolution: one positive number; above 1 favours smaller
  communities.

- seed:

  `NULL` for node order, or one whole number from 0.

## Value

A list: `memb`, a data frame with `node` and `comm` (communities
numbered 1, 2, ... in the order of their first node); `q`, its
modularity; `levels`, the levels that moved a node.

## Details

The method is a heuristic: different orders reach different local optima
(on Zachary's karate club, modularity from about 0.39 to the optimum
0.4198), so run several seeds and keep the best. Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md);
isolated nodes stay alone.

## References

Blondel, V. D., Guillaume, J.-L., Lambiotte, R. and Lefebvre, E. (2008).
Fast unfolding of communities in large networks. *Journal of Statistical
Mechanics*, P10008.
[doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

## Examples

``` r
e <- c(1, 1, 2, 4, 4, 5, 3)
f <- c(2, 3, 3, 5, 6, 6, 4)
ntlvn(e, f)
#> $memb
#>   node comm
#> 1    1    1
#> 2    2    1
#> 3    3    1
#> 4    4    2
#> 5    5    2
#> 6    6    2
#> 
#> $q
#> [1] 0.3571429
#> 
#> $levels
#> [1] 1
#> 
```
