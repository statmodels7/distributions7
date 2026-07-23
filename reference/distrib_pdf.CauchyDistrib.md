# Cauchy Probability Density Function

Computes the probability density function for the Cauchy distribution:
\$\$f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left\[1 +
\left(\dfrac{y-\mu}{\sigma}\right)^2\right\]}\$\$

## Arguments

- distrib:

  A `CauchyDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
