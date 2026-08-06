# Gumbel Log-CDF Derivatives

Closed form from the location-scale structure, as for the Gaussian:
\\\partial F/\partial\mu = -f\\ and \\\partial F/\partial\sigma = -zf\\.

## Arguments

- distrib:

  A `GumbelDistrib` object.

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

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
