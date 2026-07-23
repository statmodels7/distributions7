# Binomial Analytical Gradient

Computes the analytical gradient (first derivative) of the Binomial
log-probability with respect to the parameter \\\mu\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
n\mu}{\mu(1-\mu)}\$\$

## Arguments

- distrib:

  A `BinomialDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

## Value

A list containing the vector of first derivatives.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
