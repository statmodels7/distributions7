# Inverse-Gaussian Probability Density Function

Computes the probability density function for the Inverse-Gaussian
distribution: \$\$f(y; \mu, \phi) = \sqrt{\dfrac{1}{2\pi\phi y^3}}
\exp\left\\-\dfrac{(y-\mu)^2}{2\phi\mu^2 y}\right\\\$\$

## Arguments

- distrib:

  An `InvGaussDistrib` object.

- y:

  A numeric vector of observations (\\y \> 0\\).

- theta:

  A list containing the parameters `mu` and `phi`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`invgauss_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss_distrib.md)
