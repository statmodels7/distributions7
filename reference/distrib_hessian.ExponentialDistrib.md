# Exponential Analytical Observed Hessian

\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\mu^2} -
\dfrac{2y}{\mu^3}\$\$

## Arguments

- distrib:

  An `ExponentialDistrib` object.

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

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
