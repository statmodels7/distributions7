# Binomial Probability Mass Function

Computes the probability mass function for the Binomial distribution:
\$\$P(Y=y; \mu, n) = \dbinom{n}{y} \mu^y (1-\mu)^{n-y}\$\$

## Arguments

- distrib:

  A `BinomialDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- log:

  Logical; if `TRUE`, returns the log-probability.

## Value

A numeric vector of probability values.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
