# Gamma Probability Density Function

Computes the probability density function for the Gamma distribution:
\$\$f(y; \mu, \sigma^2) =
\dfrac{1}{\Gamma\left(\dfrac{\mu^2}{\sigma^2}\right)}
\left(\dfrac{\mu}{\sigma^2}\right)^{\dfrac{\mu^2}{\sigma^2}}
y^{\dfrac{\mu^2}{\sigma^2}-1} \exp\left(-\dfrac{\mu}{\sigma^2}
y\right)\$\$

## Arguments

- distrib:

  A `Gamma2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma2`.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector of density values.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
