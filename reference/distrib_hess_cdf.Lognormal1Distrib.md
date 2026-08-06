# Lognormal Log-CDF Hessian

Closed form. On the log scale the family is location-scale, so with \\z
= (\log q - \mu)/\sigma\\ and \\\varphi\\ the standard normal density,
\\\partial^2 F/\partial\mu^2 = -z\varphi/\sigma^2\\, \\\partial^2
F/\partial\mu\partial\sigma^2 = \varphi(1-z^2)/(2\sigma^3)\\ and
\\\partial^2 F/\partial(\sigma^2)^2 = \varphi z(3-z^2)/(4\sigma^4)\\.

## Arguments

- distrib:

  A `Lognormal1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## See also

[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
