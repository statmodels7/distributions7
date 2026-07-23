# Bernoulli Analytical Observed Hessian

Computes the analytical observed Hessian (second derivative) of the
Bernoulli log-probability with respect to the parameter \\\mu\\.

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} -
\dfrac{1-y}{(1-\mu)^2}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

## Value

A list containing the vector of second derivatives.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
