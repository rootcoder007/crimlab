# The spatial encoding of graph attention

Turns a matrix of hop distances between the nodes of a network (areas,
places, people) into the bias that Graphormer (Ying et al. 2021) adds to
the attention logits: a learnable value per distance, `b[d + 1]` (the
last for every larger distance), and `ub` for pairs that cannot reach
each other (distance `-1`). A bias inside the softmax keeps distant
nodes reachable but discouraged.

## Usage

``` r
spbias(d, b, ub = -10)
```

## Arguments

- d:

  The hop distances: a matrix of whole numbers from 0, `-1` for an
  unreachable pair.

- b:

  The bias table, indexed by distance from 0.

- ub:

  The bias of an unreachable pair: one number.

## Value

A list: `bias` (a matrix the shape of `d`) and `ub`.

## References

Ying, C., Cai, T., Luo, S., Zheng, S., Ke, G., He, D., Shen, Y. and Liu,
T.-Y. (2021). Do transformers really perform badly for graph
representation? *Advances in Neural Information Processing Systems*, 34,
28877-28888.
[doi:10.48550/arXiv.2106.05234](https://doi.org/10.48550/arXiv.2106.05234)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
d <- rbind(c(0, 1, 2), c(1, 0, -1), c(2, -1, 0))
spbias(d, b = c(1, 0.5, 0.25))
#> $bias
#>      [,1]  [,2]   [,3]
#> [1,] 1.00   0.5   0.25
#> [2,] 0.50   1.0 -10.00
#> [3,] 0.25 -10.0   1.00
#> 
#> $ub
#> [1] -10
#> 
```
