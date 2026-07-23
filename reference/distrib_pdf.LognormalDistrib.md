# Lognormal Probability Density Function

Computes the probability density function for the Lognormal
distribution: \$\$f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}}
\exp\left\\-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\\\$\$

## Arguments

- distrib:

  A `LognormalDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`lognormal_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal_distrib.md)
