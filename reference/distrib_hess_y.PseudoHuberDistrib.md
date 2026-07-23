# Pseudo-Huber Response Second Derivative

Closed-form \\\partial^2 \ell / \partial y^2 = -\nu/(\sigma^2 D^3)\\,
\\D = \sqrt{\nu + ((y-\mu)/\sigma)^2}\\.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu`, `sigma` and `nu`.

## Value

A numeric vector.

## See also

[`pseudohuber_distrib`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
