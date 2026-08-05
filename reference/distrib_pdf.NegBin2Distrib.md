# Negative Binomial Probability Mass Function

Computes the probability mass function for the Negative Binomial
distribution (NB2): \$\$P(Y=y; \mu, \theta) =
\dfrac{\Gamma(y+\theta)}{y!\\\Gamma(\theta)}
\left(\dfrac{\theta}{\theta+\mu}\right)^\theta
\left(\dfrac{\mu}{\theta+\mu}\right)^y\$\$

## Arguments

- distrib:

  A `NegBin2Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `theta`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
