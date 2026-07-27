# Gaussian Log-CDF Derivatives

Closed form, from the location-scale structure: \\\partial F/\partial\mu
= -f(y)\\ and \\\partial F/\partial\sigma = -z f(y)\\ with \\z =
(y-\mu)/\sigma\\.

## Arguments

- distrib:

  A `GaussianDistrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `sigma`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

## Value

A named list, one vector per parameter.

## See also

[`gaussian_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian_distrib.md)
