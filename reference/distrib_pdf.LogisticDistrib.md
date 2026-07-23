# Logistic Probability Density Function

Computes the probability density function for the Logistic distribution:
\$\$f(y; \mu, \sigma) =
\dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left\[1 +
\exp\left(-\dfrac{y-\mu}{\sigma}\right)\right\]^2}\$\$

## Arguments

- distrib:

  A `LogisticDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`logistic_distrib`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md)
