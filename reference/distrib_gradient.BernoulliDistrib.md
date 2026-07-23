# Bernoulli Analytical Gradient

Computes the analytical gradient (first derivative) of the Bernoulli
log-probability with respect to the parameter \\\mu\\.

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
\mu}{\mu(1-\mu)}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

## Value

A list containing the vector of first derivatives.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
