# Pseudo-Huber Response Derivatives

Closed-form derivatives of the Pseudo-Huber log-density with respect to
the response. Let \\r = y - \mu\\ and \\D = \sqrt{\nu + (r/\sigma)^2}\\:
\\\partial \ell / \partial y = -r/(\sigma^2 D)\\ and \\\partial^2 \ell /
\partial y^2 = -\nu/(\sigma^2 D^3)\\.

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
