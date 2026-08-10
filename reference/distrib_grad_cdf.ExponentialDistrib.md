# Exponential Log-CDF Derivatives

Closed form at every order from the survival function \\S =
\exp(-q/\mu)\\, whose logarithm has the partial derivatives
\\\partial^{j}L/\partial\mu^{j} = -q(-1)^{j}j!/\mu^{j+1}\\.

## Arguments

- distrib:

  An `ExponentialDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
