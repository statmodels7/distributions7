# Weibull Probability Density Function

Computes the probability density function for the Weibull distribution:
\$\$f(y; \mu, \sigma) = \dfrac{\sigma}{\mu}
\left(\dfrac{y}{\mu}\right)^{\sigma - 1}
\exp\left\\-\left(\dfrac{y}{\mu}\right)^{\sigma}\right\\\$\$

## Arguments

- distrib:

  A `WeibullDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`weibull_distrib`](https://statmodels7.github.io/distributions7/reference/weibull_distrib.md)
