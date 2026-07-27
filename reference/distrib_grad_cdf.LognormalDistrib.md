# Lognormal Log-CDF Gradient

Closed form. On the log scale the lognormal is a location-scale family,
so \\\partial F/\partial\mu = -y f(y)\\ and \\\partial
F/\partial\sigma^{2} = -y f(y) z/(2\sigma)\\ with \\z = (\log y -
\mu)/\sigma\\.

## Arguments

- distrib:

  A `LognormalDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma2`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`lognormal_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
