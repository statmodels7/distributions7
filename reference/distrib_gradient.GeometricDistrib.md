# Geometric Analytical Gradient

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y -
\mu}{\mu(1+\mu)}\$\$ the deviation from the mean over the variance.

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

A named list with the `mu` component.

## See also

[`geometric_distrib`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
