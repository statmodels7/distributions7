# Laplace Response Derivatives

Closed-form derivative of the Laplace log-density with respect to the
response, \\\partial \ell / \partial y = -\mathrm{sign}(y-\mu)/\sigma\\
(the second derivative is 0 almost everywhere). The analytic form is
provided because finite differences would be inaccurate across the kink
at \\y = \mu\\.

## Arguments

- distrib:

  A `LaplaceDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector.

## See also

[`laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
