# Exponential Analytical Gradient

\$\$\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\mu^2}\$\$ the
score of a one-parameter family written as the deviation from the mean
over the variance.

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

A named list with the `mu` component.

## See also

[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
