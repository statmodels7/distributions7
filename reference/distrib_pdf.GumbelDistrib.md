# Gumbel Probability Density Function

Computes the probability density function for the Gumbel distribution,
with \\z = (y - \mu)/\sigma\\: \$\$f(y; \mu, \sigma) = \dfrac{1}{\sigma}
\exp\left\\-z - e^{-z}\right\\\$\$

## Arguments

- distrib:

  A `GumbelDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
