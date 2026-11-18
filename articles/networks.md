# Networks of co-offending

Building a co-offending network and measuring who is central in it.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Co-offending graph from an incident table

[`coofg()`](https://rootcoder007.github.io/crimlab/reference/coofg.md):
The one-mode projection of offenders onto incidents: two offenders are
joined when they are recorded in the same incident, and the weight of
the tie is the number of incidents they share. An offender recorded
twice in one incident counts once there. Offenders who never offend with
anyone appear in `nodes` but in no edge; give `nodes` to the network
functions to keep them as isolated nodes.

``` r

g <- coofg(inc = c(1, 1, 1, 2, 2, 3), off = c("a", "b", "c", "a", "b", "d"))
g$edges
#>   from to w
#> 1    a  b 2
#> 2    a  c 1
#> 3    b  c 1
ntdeg(g$edges$from, g$edges$to, w = g$edges$w, nodes = g$nodes)
#>   node deg wdeg
#> 1    a   2    3
#> 2    b   2    3
#> 3    c   2    2
#> 4    d   0    0
```

## Degree and strength of each node

[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md):
A graph is an edge list: `from[k]` is tied to `to[k]`. The rules below
hold for every network function in crimlab.

``` r

ntdeg(c("a", "a", "b", "a"), c("b", "c", "c", "b"))
#>   node deg wdeg
#> 1    a   2    3
#> 2    b   2    3
#> 3    c   2    2
ntdeg(c("a", "a"), c("b", "c"), nodes = c("c", "b", "a", "z"), dir = TRUE)
#>   node in out deg win wout wdeg
#> 1    c  1   0   1   1    0    1
#> 2    b  1   0   1   1    0    1
#> 3    a  0   2   2   0    2    2
#> 4    z  0   0   0   0    0    0
```

## Betweenness of each node

[`ntbtw()`](https://rootcoder007.github.io/crimlab/reference/ntbtw.md):
The sum, over pairs of other nodes, of the share of their shortest paths
that pass through the node (Freeman 1977), by Brandes’ algorithm in
O(nodes x edges) time. Paths are counted in steps: weights are tie
strengths, not lengths, and are not taken. An undirected pair is counted
once. Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

``` r

ntbtw(c(1, 2, 3), c(2, 3, 4))
#>   node btw
#> 1    1   0
#> 2    2   2
#> 3    3   2
#> 4    4   0
```

## Closeness of each node

[`ntcls()`](https://rootcoder007.github.io/crimlab/reference/ntcls.md):
Distances are counted in steps (weights are not taken), along edges’
directions when `dir = TRUE` (out-distances). A disconnected graph is
handled in the two standard ways, both returned: `cls` and `ncls` count
only the nodes a node reaches (as `igraph::closeness()` does), so a node
close to the few it reaches scores high; `hrm`, the harmonic centrality,
sums `1 / distance` over every other node, an unreachable one adding 0,
so reaching more nodes always scores more. Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

``` r

ntcls(c(1, 2, 3), c(2, 3, 4), nodes = 1:5)
#>   node reach       cls ncls      hrm
#> 1    1     3 0.1666667 0.50 1.833333
#> 2    2     3 0.2500000 0.75 2.500000
#> 3    3     3 0.2500000 0.75 2.500000
#> 4    4     3 0.1666667 0.50 1.833333
#> 5    5     0       NaN  NaN 0.000000
```

## Eigenvector centrality of each node

[`ntevc()`](https://rootcoder007.github.io/crimlab/reference/ntevc.md):
The leading eigenvector of the (weighted) adjacency matrix of the
undirected graph, scaled to a maximum of 1 (Bonacich 1972). It is found
by power iteration on A + c I, `c` the largest weight, which has the
same eigenvectors and whose leading eigenvalue is strictly dominant even
on a bipartite graph (a co-offending network between two roles, a star),
where plain power iteration oscillates. The iteration starts from 1 for
every node and stops when no value changes by `tol` or more between two
steps, or after `maxit` steps (then `conv` is `FALSE`).

``` r

ntevc(c("hub", "hub", "hub", "a"), c("a", "b", "c", "b"))
#> $cent
#>   node       evc
#> 1  hub 1.0000000
#> 2    a 0.8546377
#> 3    b 0.8546377
#> 4    c 0.4608111
#> 
#> $value
#> [1] 2.170086
#> 
#> $iter
#> [1] 33
#> 
#> $conv
#> [1] TRUE
```

## k-core of each node

[`ntcore()`](https://rootcoder007.github.io/crimlab/reference/ntcore.md):
The coreness of a node is the largest k such that it belongs to a
subgraph in which every node has at least k neighbours (Seidman 1983),
computed in O(nodes + edges) by the algorithm of Batagelj and Zaversnik.
The cores are nested: the (k + 1)-core lies inside the k-core
(`proof("ntcore")`). Edges are taken as undirected, without weights.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

``` r

ntcore(c(1, 1, 1, 2, 2, 3, 4), c(2, 3, 4, 3, 4, 4, 5))
#>   node core
#> 1    1    3
#> 2    2    3
#> 3    3    3
#> 4    4    3
#> 5    5    1
```

## Connected components

[`ntcmp()`](https://rootcoder007.github.io/crimlab/reference/ntcmp.md):
Weak components join nodes linked by edges in either direction; strong
components (Tarjan’s algorithm, without recursion) join nodes that reach
each other along the edges’ directions. Undirected, the two coincide.
Components are numbered 1, 2, … in the order of their first node. Graph
rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

``` r

ntcmp(c("a", "b", "c", "x"), c("b", "c", "a", "y"), dir = TRUE)
#>   node weak strong
#> 1    a    1      1
#> 2    b    1      1
#> 3    c    1      1
#> 4    x    2      2
#> 5    y    2      3
```

## Communities by the Louvain method

[`ntlvn()`](https://rootcoder007.github.io/crimlab/reference/ntlvn.md):
Groups nodes to raise the modularity (see
[`ntmod()`](https://rootcoder007.github.io/crimlab/reference/ntmod.md))
of the undirected, weighted graph (Blondel et al. 2008). Each level
visits the nodes in node order when `seed` is `NULL`, or in a random
order drawn from `seed` (one stream for the whole run), and moves each
to the neighbouring community with the largest gain in modularity, the
first met among equal gains, when that gain beats staying by more than
1e-12 of the node’s strength; passes repeat until none moves a node. The
communities then become the nodes of the next level, and a level that
moves nothing ends the run. The result is the same on every platform for
the same input and seed.

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
```

## Modularity of a partition

[`ntmod()`](https://rootcoder007.github.io/crimlab/reference/ntmod.md):
Q = sum over communities of in / 2m - res (tot / 2m)^2, where in is the
weight of edges inside the community (counted from both ends), `tot` the
summed strength of its nodes and 2m the total strength (Newman and
Girvan 2004), on the undirected, weighted graph. At `res = 1` it lies
between -1/2 and 1 (`proof("ntmod")`). `NaN` for a graph without edges.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

``` r

ntmod(c(1, 1, 2, 4, 4, 5, 3), c(2, 3, 3, 5, 6, 6, 4),
      memb = c(1, 1, 1, 2, 2, 2))
#> [1] 0.3571429
```

## Brokerage roles

[`ntbrk()`](https://rootcoder007.github.io/crimlab/reference/ntbrk.md):
Counts, for each node v, the two-paths `i -> v -> j` with `i != j` and
no tie `i -> j` (v brokers between them), by the roles of Gould and
Fernandez (1989), given each node’s group: `coord` (coordinator: i, v
and j in one group), `itin` (itinerant: i and j in one group, v in
another), `repr` (representative: v in i’s group, j outside), `gate`
(gatekeeper: v in j’s group, i outside) and `liais` (liaison: three
groups). Undirected, each tie runs both ways, so each unordered pair of
neighbours counts twice and `repr` equals `gate`. Weights are not taken.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).
Raw counts only: the expected counts and z-scores of Gould and Fernandez
are not given.

``` r

ntbrk(c("a", "b", "b"), c("b", "c", "d"),
      grp = c(a = "x", b = "y", c = "y", d = "z"), dir = TRUE)
#>   node grp coord itin repr gate liais total
#> 1    a   x     0    0    0    0     0     0
#> 2    b   y     0    0    0    1     1     2
#> 3    c   y     0    0    0    0     0     0
#> 4    d   z     0    0    0    0     0     0
```

## Cheeger bounds on a network’s bottleneck

[`ntchg()`](https://rootcoder007.github.io/crimlab/reference/ntchg.md):
The conductance of a set of nodes S is the weight of the edges leaving
it over the smaller of its volume and that of the rest (volume: summed
strength); the graph’s conductance h is the least over all sets, a
measure of how cleanly the network splits in two. Finding h is hard, but
the second-smallest eigenvalue `lambda2` of the normalised Laplacian I -
D^-1/2 A D^-1/2 brackets it: lambda2 / 2 \<= h \<= sqrt(2 lambda2)
(Cheeger’s inequality; the lower half is `proof("ntchg")`). The sweep
over the eigenvector finds a set whose conductance (`sweep`) is at most
sqrt(2 lambda2); h is at most `sweep`.

``` r

# two triangles joined by one edge
ntchg(c(1, 1, 2, 4, 4, 5, 3), c(2, 3, 3, 5, 6, 6, 4))
#> $bound
#>   lambda2     lower     upper     sweep 
#> 0.2046664 0.1023332 0.6397911 0.1428571 
#> 
#> $set
#> [1] "4" "5" "6"
#> 
#> $iter
#> [1] 31
#> 
#> $conv
#> [1] TRUE
```
