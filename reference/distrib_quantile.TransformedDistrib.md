# Transformed Quantile Function

\\Q_Y(p) = g(Q_X(p))\\, with tails swapped for decreasing
transformations.

## Arguments

- distrib:

  A `TransformedDistrib` object.

- p:

  A numeric vector of probabilities.

- theta:

  A list of the parent's parameters.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le p)\\.

- log.p:

  Logical; if `TRUE`, probabilities are given as logs.

## See also

[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)
