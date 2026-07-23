# Poisson Analytical Expected Hessian

Computes the analytical expected Hessian of the Poisson log-probability
with respect to the parameter \\\mu\\.

\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu}\$\$

## Arguments

- distrib:

  A `PoissonDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

## Value

A list containing the vector of expected second derivatives.

## See also

[`poisson_distrib`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
