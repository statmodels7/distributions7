# Exponential Log-CDF Derivatives

Closed form. With \\F = 1 - e^{-q/\mu}\\, \\\partial F/\partial\mu =
-(q/\mu) f\\ and \\\partial^2 F/\partial\mu^2 = (q/\mu^2) f (2 -
q/\mu)\\.

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

## Value

A named list with one element.

## See also

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
