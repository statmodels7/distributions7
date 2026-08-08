# Laplace Log-CDF Derivatives in Location and Rate

Closed form: \\\partial F/\partial\mu = -f(q)\\ and \\\partial
F/\partial\lambda = (q-\mu)\\f(q)/\lambda\\. The second derivatives
inherit the kink at \\q = \mu\\ exactly as those of
[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
do.

## Arguments

- distrib:

  A `Laplace2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `lambda`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md)
