# Skew t Response Derivative

Closed form: \\\partial \ell / \partial y = D/\sigma\\, which is minus
the derivative in \\\mu\\, as it must be for a location family.

## Arguments

- distrib:

  A `SkewTDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `sigma`, `alpha` and `nu`.

## Value

A numeric vector.

## See also

[`skewt_distrib`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
