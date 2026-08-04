# Geometric Quantile Function

The generalized inverse \\Q(p) = \min\\y \in \mathbb{N}\_0 : F(y) \ge
p\\\\.

## Arguments

- distrib:

  A `GeometricDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), `p` is \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, `p` is given as its logarithm.

## Value

A numeric vector of quantiles.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
