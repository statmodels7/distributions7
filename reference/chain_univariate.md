# One Univariate Chain Rule Over a Multivariate Map

Evaluates \\\partial^{S} f(u(v))\\ for a scalar function \\f\\ of one
coordinate, by the partition form of Faa di Bruno \$\$\partial^{S} f(u)
= \sum\_{\pi} f^{(\lvert\pi\rvert)}(u) \prod\_{B \in \pi} \partial^{B}
u,\$\$ the sum running over the set partitions of the multiset \\S\\.
Everything is vectorized over the coordinates, so one call serves every
\\j\\ at once.

## Usage

``` r
chain_univariate(tuple, fd, ud)
```

## Arguments

- tuple:

  An integer vector of free-value indices, with repeats.

- fd:

  A list whose \\m\\-th element is \\f^{(m)}\\ evaluated at every
  coordinate.

- ud:

  A named list of the map's derivative tensors, keyed by sorted tuple as
  [`simplex_map_tensors`](https://statmodels7.github.io/distributions7/reference/simplex_map_tensors.md)
  returns them.

## Value

A numeric vector over the coordinates.
