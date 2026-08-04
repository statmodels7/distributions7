# Skew Normal Response Derivative

Closed form: \\\partial \ell / \partial y = (\alpha R - z)/\sigma\\,
which is minus the derivative in \\\mu\\, as it must be for a location
family.

## Arguments

- distrib:

  A `SkewNormalDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma` and `alpha`.

## Value

A numeric vector.

## See also

[`skewnormal_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal_distrib.md)
