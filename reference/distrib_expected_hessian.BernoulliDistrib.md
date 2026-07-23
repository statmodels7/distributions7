# Bernoulli Analytical Expected Hessian

Computes the analytical expected Hessian of the Bernoulli
log-probability with respect to the parameter \\\mu\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu(1-\mu)}\$\$

## Arguments

- distrib:

  A `BernoulliDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

## Value

A list containing the vector of expected second derivatives.

## See also

[`bernoulli_distrib`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md)
