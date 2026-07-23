# Binomial Analytical Observed Hessian

Computes the analytical observed Hessian (second derivative) of the
Binomial log-probability with respect to the parameter \\\mu\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} -
\dfrac{n-y}{(1-\mu)^2}\$\$

## Arguments

- distrib:

  A `BinomialDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

## Value

A list containing the vector of second derivatives.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
