# Exponential Quantile Function

\$\$Q(p; \mu) = -\mu \log(1 - p)\$\$

## Arguments

- distrib:

  An `ExponentialDistrib` object.

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

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
