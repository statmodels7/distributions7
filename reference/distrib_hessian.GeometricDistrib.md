# Geometric Analytical Observed Hessian

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} +
\dfrac{y+1}{(1+\mu)^2}\$\$

## Arguments

- distrib:

  A `GeometricDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameter `mu`.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- ...:

  Unused.

## Value

A named list with the `mu_mu` component.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
