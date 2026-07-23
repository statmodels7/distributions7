# Cauchy Response Derivatives

Closed-form derivatives of the Cauchy log-density with respect to the
response. Let \\r = y - \mu\\ and \\d = \sigma^2 + r^2\\: \\\partial
\ell / \partial y = -2r/d\\ and \\\partial^2 \ell / \partial y^2 =
2(r^2 - \sigma^2)/d^2\\.

## Arguments

- distrib:

  A `CauchyDistrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector.

## See also

[`cauchy_distrib`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
