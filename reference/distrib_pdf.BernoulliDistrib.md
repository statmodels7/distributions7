# Bernoulli Probability Mass Function

Computes the probability mass function for the Bernoulli distribution:
\$\$P(Y=y; \mu) = \mu^y (1-\mu)^{1-y}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- y:

  A numeric vector of observations (`0` or `1`).

- theta:

  A list containing the parameter `mu`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
