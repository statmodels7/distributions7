# Poisson Log-CDF Gradient

Closed form: the sum defining \\F\\ telescopes, leaving \\\partial
F(k)/\partial\mu = -f(k)\\.

## Arguments

- distrib:

  A `PoissonDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list with one element.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
