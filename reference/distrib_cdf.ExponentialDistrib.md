# Exponential Distribution Function

\$\$F(q; \mu) = 1 - e^{-q/\mu}\$\$

## Arguments

- distrib:

  An `ExponentialDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing the parameter `mu`.

- lower.tail:

  Logical; if `TRUE` (default), \\P(Y \le q)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned as logarithms.

## Value

A numeric vector of cumulative probabilities.

## See also

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
