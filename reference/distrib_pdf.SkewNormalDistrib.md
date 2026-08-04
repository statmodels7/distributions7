# Skew Normal Probability Density Function

Computes the probability density function, with \\z = (y-\mu)/\sigma\\:
\$\$f(y; \mu, \sigma, \alpha) = \dfrac{2}{\sigma}\\\phi(z)\\\Phi(\alpha
z)\$\$

## Arguments

- distrib:

  A `SkewNormalDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`skewnormal_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md)
