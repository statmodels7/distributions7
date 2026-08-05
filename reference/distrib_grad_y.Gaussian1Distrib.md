# Gaussian Response Derivatives

Closed-form derivatives of the Gaussian log-density with respect to the
response: \\\partial \ell / \partial y = -(y-\mu)/\sigma^2\\ and
\\\partial^2 \ell / \partial y^2 = -1/\sigma^2\\.

## Arguments

- distrib:

  A `Gaussian1Distrib` object.

- y:

  A numeric vector of observations.

- theta:

  A list containing the parameters `mu` and `sigma`.

## Value

A numeric vector.

## See also

[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
