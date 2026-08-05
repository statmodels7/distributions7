# Gaussian Probability Density Function

Computes the probability density function for the Gaussian distribution:
\$\$f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\sigma}
\exp\left\\-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\\\$\$

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
