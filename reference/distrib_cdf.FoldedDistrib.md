# Folded Distribution Function

\$\$P(\|Y\| \le q) = F(q) - F(-q)\$\$

## Arguments

- distrib:

  A `FoldedDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical; if `TRUE` (default), \\P(\|Y\| \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
