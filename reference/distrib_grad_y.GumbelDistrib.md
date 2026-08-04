# Gumbel Response Derivative

Closed form: \\\partial \ell / \partial y = (w - 1)/\sigma\\, which is
minus the derivative in \\\mu\\, as it must be for a location family.

## Arguments

- distrib:

  A `GumbelDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector.

## See also

[`gumbel_distrib`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
