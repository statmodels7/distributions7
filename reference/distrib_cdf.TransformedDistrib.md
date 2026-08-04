# Transformed Cumulative Distribution Function

\\F_Y(q) = F_X(g^{-1}(q))\\, with tails swapped for decreasing
transformations.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list of the parent's parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logs.

## Value

A numeric vector of cumulative probabilities.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
